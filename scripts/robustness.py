"""Pre-registered robustness grid (research design, section 8).

One month loop covers the variants that need no new inputs: hedge moneyness,
butterfly reading, contributor, stale butterflies treated as missing, event
exclusions, implementable returns with option half-spreads, dollar-carry and
ten-currency portfolios, log returns and previous-day spot for payoffs. Each
variant reports E1 (regime difference in φ), E2 (bias-corrected slope of
HML^U on φ and the Clark–West statistic) and E3 (mean returns, skew term and
θ_UB). The base variant must reproduce E1 and Stage 4; the script stops if it
does not.

Definitions fixed here, before any variant was estimated:

- Stale butterflies: the audit rule. A currency-month is dropped when the
  selected 25Δ butterfly quote lies in a run of five or more identical
  two-sided business-day quotes.
- Exclusions: month-ends whose quote date or return window falls in the
  excluded period. March 2020 removes the February and March 2020
  month-ends; the CHF exclusion removes CHF from the ranking at the
  November 2014 to March 2015 month-ends.
- Implementable returns: forwards traded at bid and ask, closed at the spot
  bid or ask on the payoff date; options bought at the smile volatility plus
  k times the ATM half-spread, converted to premium by vega.
- Dollar carry: long all available currencies when their average forward
  discount is positive, short all otherwise, equal weights.
- Ten currencies: USD enters the ranking with fd = 0; a USD leg has zero
  return and no option.
- Log returns: each leg's forward return is replaced by ±(s_{t+1} − f_t);
  option payoffs and premia are unchanged.
- Previous-day spot: option payoffs use the spot of the New York business
  day before expiry; forward returns keep the expiry-date spot.
- Clark–West: expanding window from the primary start, first forecast for
  January 2017, the same rule for every variant including the base.

The three-month tenor and vanna–volga smiles need new inputs and are run
separately. Outputs are written to data/private/results/.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from estimate_stage4 import BDAY, clark_west, daily_spot_mid, spot_on, stambaugh_bootstrap  # noqa: E402

from qef.data.panel import on_business_days, read_raw, sample_month_ends, stale_mask, two_sided  # noqa: E402
from qef.data.smile_inputs import option_dates  # noqa: E402
from qef.fx.conventions import G10  # noqa: E402
from qef.fx.crash import forward_discount, forward_return, leg_skew_cost, option_payoff, protective_option  # noqa: E402
from qef.fx.gk import atm_dns_strike, forward_premium, vega  # noqa: E402
from qef.fx.sabr import sabr_vol  # noqa: E402
from qef.stats.bootstrap import bootstrap_distribution  # noqa: E402
from qef.stats.hac import mean_and_se, ols_hac  # noqa: E402

REGIME_BREAK = pd.Timestamp("2022-01-01")
PRIMARY_START = pd.Timestamp("2013-05-31")
STALE_RUN = 5


@dataclass
class Variant:
    name: str
    contributor: str = ""  # "" composite, "fn" Fenics
    reading: str = "market"
    hedge: str = "10d"  # "10d", "25d", "atm"
    portfolio: str = "hml"  # "hml", "dollar", "ten"
    exclude_months: set = field(default_factory=set)
    exclude_pairs: set = field(default_factory=set)  # (currency, month-end)
    stale_missing: bool = False
    cost_k: float | None = None
    log_returns: bool = False
    payoff_lag: int = 0  # New York business days before expiry for the payoff spot


def stale_flags(raw: Path, month_ends) -> set:
    """(currency, month-end) whose selected 25Δ butterfly quote is stale under the audit rule."""
    out = set()
    for c in G10:
        f = on_business_days(read_raw(raw / "vol_bf25" / f"{c}1MBF=.csv"))
        q = f[two_sided(f)]
        mask = stale_mask(q, STALE_RUN)
        me = sample_month_ends(f, month_ends)
        me = me[me.status != "missing"]
        for m, src in me.source_date.items():
            if bool(mask.get(src, False)):
                out.add((c, m))
    return out


def quote_sides(raw: Path, month_ends):
    """Month-end forward outright bid/ask and ATM half-spread (decimal vol); daily spot bid/ask."""
    sides, daily = {}, {}
    for c, conv in G10.items():
        sp = sample_month_ends(read_raw(raw / "spot" / f"{c}=.csv"), month_ends)
        fw = sample_month_ends(read_raw(raw / "forward" / f"{c}1M=.csv"), month_ends)
        at = sample_month_ends(read_raw(raw / "vol_atm" / f"{c}1MO=.csv"), month_ends)
        for m in month_ends:
            sides[(c, m)] = {"F_bid": sp.loc[m, "bid"] + fw.loc[m, "bid"] * conv.pip_factor,
                             "F_ask": sp.loc[m, "ask"] + fw.loc[m, "ask"] * conv.pip_factor,
                             "atm_hs": (at.loc[m, "ask"] - at.loc[m, "bid"]) / 200.0}
        f = read_raw(raw / "spot" / f"{c}=.csv")
        f = f[~f.index.duplicated(keep="last")].sort_index()
        daily[c] = f.loc[two_sided(f), ["BID", "ASK"]].astype(float)
    return sides, daily


def select_legs(fd: dict, portfolio: str):
    """(currency, long?, weight) for each leg."""
    if portfolio == "dollar":
        long_all = np.mean(list(fd.values())) > 0
        return [(c, long_all, 1 / len(fd)) for c in fd]
    if portfolio == "ten":
        fd = {**fd, "USD": 0.0}
    order = sorted(fd, key=fd.get)
    return [(c, True, 1 / 3) for c in order[-3:]] + [(c, False, 1 / 3) for c in order[:3]]


def implementable_return(c, is_long, sd, spot_row):
    """Forward return with the forward at bid/ask and the closing spot at bid/ask."""
    if G10[c].usd_base:  # long j = sell USD forward at the bid, buy USD spot at the ask
        return sd["F_bid"] / spot_row.ASK - 1 if is_long else 1 - sd["F_ask"] / spot_row.BID
    return spot_row.BID / sd["F_ask"] - 1 if is_long else 1 - spot_row.ASK / sd["F_bid"]


def run_variant(v: Variant, inputs, panel, spots, sides, daily, stale):
    inp = inputs.set_index(["currency", "date"])
    fits = panel[(panel.reading == v.reading) & panel.status.isin(["ok", "not_converged"])].set_index(["currency", "date"])
    delta = {"10d": 0.10, "25d": 0.25}.get(v.hedge)
    rows = []
    for t in sorted(pd.Timestamp(x) for x in inputs.date.unique()):
        if t < PRIMARY_START or t in v.exclude_months:
            continue
        avail = {}
        for c in G10:
            if (c, t) not in inp.index or (c, t) not in fits.index or not inp.loc[(c, t)].complete_calib:
                continue
            if (c, t) in v.exclude_pairs or (v.stale_missing and (c, t) in stale):
                continue
            avail[c] = (inp.loc[(c, t)], fits.loc[(c, t)])
        if len(avail) < 6:
            rows.append({"date": t, "status": "too_few_currencies"})
            continue
        fd = {c: forward_discount(r.S, r.F, G10[c].usd_base) for c, (r, _) in avail.items()}
        acc = {"C": 0.0, "FD": 0.0, "U": 0.0, "H": 0.0, "c3": 0.0}
        status = "ok"
        for c, is_long, w in select_legs(fd, v.portfolio):
            if c == "USD":
                continue
            r, f = avail[c]
            vol = lambda K, r=r, f=f: sabr_vol(K, r.F, r.tau, f.alpha, f.rho, f.nu, 1.0)
            if delta is None:
                K = atm_dns_strike(r.F, r.atm, r.tau, G10[c].delta)
                vs = vf = float(forward_premium(r.F, K, r.atm, r.tau, protective_option(c, is_long))) / r.F
            else:
                try:
                    K, vs, vf = leg_skew_cost(vol, r.F, r.tau, r.df_base, c, r.atm, is_long, delta)
                except ValueError:
                    status = f"inversion_failed:{c}"
                    break
            acc["C"] += w * (vs - vf)
            acc["FD"] += w * (fd[c] if is_long else -fd[c])
            acc["c3"] -= w * (vs - vf)
            if status != "ok":
                continue
            expiry = option_dates(t, c)[2]
            S_end, _ = spot_on(spots[c], expiry)
            S_pay, _ = spot_on(spots[c], expiry - v.payoff_lag * BDAY) if v.payoff_lag else (S_end, 0)
            if not (np.isfinite(S_end) and np.isfinite(S_pay)):
                status = "return_unrealised"
                continue
            prem = vs
            if v.cost_k is not None:
                sd = sides[(c, t)]
                rx = implementable_return(c, is_long, sd, daily[c].loc[:expiry].iloc[-1])
                prem += v.cost_k * sd["atm_hs"] * float(vega(r.F, K, float(vol(K)), r.tau, 1.0)) / r.F
            elif v.log_returns:
                rl = np.log(forward_return(c, True, S_end, r.F) + 1.0)
                rx = rl if is_long else -rl
            else:
                rx = forward_return(c, is_long, S_end, r.F)
            acc["U"] += w * rx
            acc["H"] += w * (rx + option_payoff(c, is_long, K, r.F, S_pay) - prem)
        if status.startswith("inversion_failed"):
            rows.append({"date": t, "status": status})
            continue
        if status == "return_unrealised":
            acc["U"] = acc["H"] = np.nan
        phi = acc["C"] / acc["FD"] if acc["FD"] > 0 and delta is not None else np.nan
        rows.append({"date": t, "status": "ok", **acc, "phi": phi})
    return pd.DataFrame(rows)


def summarise(df: pd.DataFrame, B: int) -> dict:
    ok = df[df.status == "ok"]
    ret = ok.dropna(subset=["U"])
    out = {"n_phi": int(ok.phi.notna().sum()), "n_returns": len(ret)}
    s = ok.dropna(subset=["phi"])
    if len(s):
        D = (s.date >= REGIME_BREAK).astype(float).to_numpy()
        reg = ols_hac(s.phi.to_numpy(), np.column_stack([np.ones(len(s)), D]))
        a, b = s.phi[D == 0].to_numpy(), s.phi[D == 1].to_numpy()
        boot = bootstrap_distribution(lambda x, y: y.mean() - x.mean(), [a, b], B=B)
        rp = ret.dropna(subset=["phi"])
        e2 = stambaugh_bootstrap(rp.U.to_numpy(), rp.phi.to_numpy(), B)
        cw = clark_west(rp.date.to_numpy(), rp.U.to_numpy(), rp.phi.to_numpy())
        out.update(phi_zero_rate=a.mean(), phi_hiking=b.mean(), phi_diff=reg["beta"][1], phi_diff_se=reg["se"][1],
                   phi_diff_boot_lo=np.percentile(boot, 2.5), phi_diff_boot_hi=np.percentile(boot, 97.5),
                   E2_b_bc=e2["b_bias_corrected"], E2_p=e2["p_one_sided_b_gt_0"],
                   CW_t=cw["cw_t"], CW_p=cw["cw_p_one_sided"], CW_n=cw["n_forecasts"])
    u, h, c3 = (mean_and_se(ret[k].to_numpy()) for k in ("U", "H", "c3"))
    out.update(C_skew_bp=1e4 * ok.C.mean(), U_bp=1e4 * u["mean"], U_t=u["mean"] / u["se"], H_bp=1e4 * h["mean"],
               H_t=h["mean"] / h["se"], skew_term_bp=1e4 * c3["mean"], skew_term_se_bp=1e4 * c3["se"],
               theta_UB=1 - h["mean"] / u["mean"])
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--bootstrap", type=int, default=1999)
    args = p.parse_args()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    raw = ROOT / "data" / "private" / "lseg" / args.retrieval_date / "raw"
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    data = {"": (load("smile_inputs.csv"), load("smile_panel.csv")),
            "fn": (load("smile_inputs_fn.csv"), load("smile_panel_fn.csv"))}
    month_ends = sorted(pd.Timestamp(x) for x in data[""][0].date.unique() if pd.Timestamp(x) >= PRIMARY_START)
    spots = {c: daily_spot_mid(raw, c) for c in G10}
    sides, daily = quote_sides(raw, month_ends)
    stale = stale_flags(raw, month_ends)
    in_period = lambda lo, hi: {m for m in month_ends if pd.Timestamp(lo) <= m <= pd.Timestamp(hi)}
    variants = [
        Variant("base"),
        Variant("hedge 25d", hedge="25d"),
        Variant("hedge ATM", hedge="atm"),
        Variant("smile-strangle reading", reading="smile"),
        Variant("Fenics quotes", contributor="fn"),
        Variant("stale butterflies missing", stale_missing=True),
        Variant("excluding March 2020", exclude_months=in_period("2020-02-01", "2020-03-31")),
        Variant("excluding CHF Dec 2014-Mar 2015", exclude_pairs={("CHF", m) for m in in_period("2014-11-01", "2015-03-31")}),
        Variant("implementable, k = 1", cost_k=1.0),
        Variant("implementable, k = 1.5", cost_k=1.5),
        Variant("implementable, k = 2", cost_k=2.0),
        Variant("dollar carry", portfolio="dollar"),
        Variant("ten currencies", portfolio="ten"),
        Variant("log returns", log_returns=True),
        Variant("previous-day spot for payoffs", payoff_lag=1),
    ]
    e1 = pd.read_csv(d / "e1_summary.csv")
    e1 = e1[(e1["sample"] == "primary") & (e1.reading == "market") & (e1.delta == 0.10) & (e1.variable == "phi")].iloc[0]
    s4 = load("stage4_months.csv")
    s4 = s4[(s4["sample"] == "primary") & (s4.status == "ok")]
    results, series = {}, []
    for v in variants:
        df = run_variant(v, *data[v.contributor], spots, sides, daily, stale)
        df["variant"] = v.name
        series.append(df)
        results[v.name] = summarise(df, args.bootstrap)
        print(f"done: {v.name}", flush=True)
        if v.name == "base":
            base = results["base"]
            checks = {"E1 difference": (base["phi_diff"], e1.difference), "mean HML^U": (base["U_bp"], 1e4 * s4.U.mean()),
                      "mean HML^H": (base["H_bp"], 1e4 * s4.H10.mean()),
                      "skew term": (base["skew_term_bp"], 1e4 * s4.c3.mean())}
            for k, (x, y) in checks.items():
                if not abs(x - y) < 1e-8:
                    raise SystemExit(f"base variant does not reproduce {k}: {x} vs {y}")

    pd.concat(series, ignore_index=True).to_csv(d / "robustness_series.csv", index=False)
    table = pd.DataFrame(results).T
    table.to_csv(d / "robustness_summary.csv")
    lines = ["# Robustness grid (restricted)", "",
             f"Base reproduces E1 and Stage 4 ({', '.join(checks)}). Bootstrap draws: {args.bootstrap}.",
             f"Stale 25-delta butterfly currency-months in the primary sample: {len(stale)}", "", "```",
             table.to_string(float_format=lambda x: f"{x:.4g}"), "```", ""]
    (d / "robustness_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
