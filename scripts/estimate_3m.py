"""Robustness variant 2: three-month tenor with quarterly rebalancing (research design, section 8).

The rules below were fixed before estimation (research log, 24 September 2026,
choices for the remaining Stage 5 variants); where the log is silent they are
those of the one-month scripts (``estimate_e1.py``, ``estimate_stage4.py``,
``robustness.py``).

- Inputs: composite three-month smiles calibrated by ``calibrate_smiles.py
  --tenor 3M`` (market reading, beta = 1), whose forwards and rates come from
  the retrieval of 24 September 2026 (``build_month_end_inputs``).
- Rebalancing: quarterly and non-overlapping, one position per quarter. The
  primary phase uses calendar-quarter month-ends (the last New York business
  day of March, June, September and December) from 28 June 2013. The phases
  of January, April, July and October and of February, May, August and
  November month-ends are supplementary. A date is used if at least six
  currencies have the three-month calibration set; non-converged smiles are
  used and counted.
- Legs: 3L/3S on the three-month forward discount fd = log(X/F_X).
- Hedge: each leg's protective option at the 10Δ strike of the three-month
  market-reading SABR smile, bought at the smile premium carried to delivery.
- Returns: forward returns and option payoffs at the spot on the three-month
  option expiry date (``spot_on``: substitution within five New York business
  days, and windows ending after the last spot observation excluded). The
  option dates follow the spot and delivery conventions, so an expiry can fall
  a few days after the next rebalancing date, as at one month; such quarters
  are counted.

Estimands, per quarter and not annualised:

- E1: φ = C_skew/FD with C_skew = (1/3) Σ_legs (V_smile − V_flat) and
  FD = (1/3)(Σ_long fd − Σ_short fd); regime means (zero-rate to December
  2021, hiking from January 2022), their difference with the Newey–West
  standard error, and the stationary-bootstrap 95% percentile interval
  resampling each regime separately.
- E2 (in sample): HML^U = a + b φ + e with Newey–West errors and the residual
  bootstrap under the null (``estimate_stage4.stambaugh_bootstrap``).
- E3: means of HML^U, HML^H and the skew term −(1/3) Σ_legs (V_smile − V_flat)
  with Newey–West t-statistics, θ_UB = 1 − mean(HML^H)/mean(HML^U), and its
  test-inversion set (design, section 7).

Check: run at one month over the primary month-ends, the same loop reproduces
φ of E1 and HML^U, HML^H and the skew term of Stage 4 within 1e-12; the script
stops if it does not.

Data diagnostics: the first month-end at which all nine currencies have
complete three-month inputs including the 10Δ quotes, and the counts of
substituted and missing quotes. Outputs are LSEG-derived and are written to
data/private/results/2026-09-24/. Reproduction, after the three-month
calibration (``calibrate_smiles.py --tenor 3M``):

    .venv/bin/python scripts/estimate_3m.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from estimate_stage4 import daily_spot_mid, spot_on, stambaugh_bootstrap, theta_confidence_set  # noqa: E402
from robustness import PRIMARY_START, REGIME_BREAK  # noqa: E402

from qef.data.panel import ny_month_ends  # noqa: E402
from qef.data.smile_inputs import QUOTE_FIELDS, build_month_end_inputs, option_dates  # noqa: E402
from qef.fx.conventions import G10  # noqa: E402
from qef.fx.crash import forward_discount, forward_return, leg_skew_cost, option_payoff, rank_legs  # noqa: E402
from qef.fx.sabr import sabr_vol  # noqa: E402
from qef.stats.bootstrap import bootstrap_distribution  # noqa: E402
from qef.stats.hac import mean_and_se, ols_hac  # noqa: E402

MONTHS = 3
N_LEGS = 3
DELTA = 0.10
PRIMARY_FIRST = pd.Timestamp("2013-06-28")
EXTRA_RAW_ROOT = "data/private/lseg/2026-09-24/raw"
PHASES = {"Mar/Jun/Sep/Dec": (3, 6, 9, 12), "Jan/Apr/Jul/Oct": (1, 4, 7, 10), "Feb/May/Aug/Nov": (2, 5, 8, 11)}


def quarter_rows(dates, inputs, panel, spots, months=MONTHS):
    """One row per rebalancing date: legs, E1 terms and the returns realised at the ``months`` expiry."""
    inp = inputs.set_index(["currency", "date"])
    fits = panel[(panel.reading == "market") & panel.status.isin(["ok", "not_converged"])].set_index(["currency", "date"])
    rows = []
    for t in dates:
        avail = {}
        for c in G10:
            if (c, t) in inp.index and (c, t) in fits.index and inp.loc[(c, t)].complete_calib:
                avail[c] = (inp.loc[(c, t)], fits.loc[(c, t)])
        if len(avail) < 2 * N_LEGS:
            rows.append({"date": t, "status": "too_few_currencies", "n_ccy": len(avail)})
            continue
        fd = {c: forward_discount(r.S, r.F, G10[c].usd_base) for c, (r, _) in avail.items()}
        longs, shorts = rank_legs(fd, N_LEGS)
        acc = {"C_skew": 0.0, "FD": 0.0, "U": 0.0, "H": 0.0, "c3": 0.0}
        status, spot_sub, n_nc, last_expiry = "ok", 0, 0, t
        for c, is_long in [(c, True) for c in longs] + [(c, False) for c in shorts]:
            r, f = avail[c]
            vol = lambda K, r=r, f=f: sabr_vol(K, r.F, r.tau, f.alpha, f.rho, f.nu, 1.0)
            try:
                K, vs, vf = leg_skew_cost(vol, r.F, r.tau, r.df_base, c, r.atm, is_long, DELTA)
            except ValueError:
                status = f"inversion_failed:{c}"
                break
            n_nc += int(f.status == "not_converged")
            w = 1.0 / N_LEGS
            acc["C_skew"] += w * (vs - vf)
            acc["FD"] += w * (fd[c] if is_long else -fd[c])
            acc["c3"] -= w * (vs - vf)
            expiry = option_dates(t, c, months)[2]
            last_expiry = max(last_expiry, expiry)
            S_end, sub = spot_on(spots[c], expiry)
            if not np.isfinite(S_end):
                status = "return_unrealised"
                continue
            spot_sub += sub
            rx = forward_return(c, is_long, S_end, r.F)
            acc["U"] += w * rx
            acc["H"] += w * (rx + option_payoff(c, is_long, K, r.F, S_end) - vs)
        if status.startswith("inversion_failed"):
            rows.append({"date": t, "status": status, "n_ccy": len(avail)})
            continue
        if status == "return_unrealised":
            acc["U"] = acc["H"] = np.nan
        phi = acc["C_skew"] / acc["FD"] if acc["FD"] > 0 else np.nan
        rows.append({"date": t, "status": "ok", "n_ccy": len(avail), "n_not_converged": n_nc,
                     "longs": " ".join(longs), "shorts": " ".join(shorts), **acc, "phi": phi,
                     "last_expiry": last_expiry, "spot_sub": spot_sub, "returns_realised": status == "ok"})
    return pd.DataFrame(rows)


def summarise(q: pd.DataFrame, B: int) -> dict:
    ok = q[q.status == "ok"]
    s = ok.dropna(subset=["phi"])
    ret = ok[ok.returns_realised.astype(bool)]
    out = {"n_dates": len(q), "n_too_few_currencies": int((q.status == "too_few_currencies").sum()),
           "n_inversion_failed": int(q.status.str.startswith("inversion_failed").sum()),
           "n_phi": len(s), "n_returns": len(ret), "first": s.date.min(), "last_phi": s.date.max(),
           "last_return": ret.date.max(), "legs_not_converged": int(ok.n_not_converged.sum()),
           "payoff_spot_substitutions": int(ret.spot_sub.sum())}
    nxt = ok.date.shift(-1)
    over = (ok.last_expiry > nxt) & nxt.notna()
    out["n_expiry_after_next_date"] = int(over.sum())
    out["max_overlap_days"] = int((ok.last_expiry - nxt)[over].dt.days.max()) if over.any() else 0

    # E1
    D = (s.date >= REGIME_BREAK).to_numpy()
    a, b = s.phi[~D].to_numpy(), s.phi[D].to_numpy()
    reg = ols_hac(s.phi.to_numpy(), np.column_stack([np.ones(len(s)), D.astype(float)]))
    boot = bootstrap_distribution(lambda x, y: y.mean() - x.mean(), [a, b], B=B)
    ra, rb = mean_and_se(a), mean_and_se(b)
    out.update(n_zero_rate=len(a), n_hiking=len(b), phi_zero_rate=ra["mean"], phi_zero_rate_se=ra["se"],
               phi_hiking=rb["mean"], phi_hiking_se=rb["se"], phi_diff=reg["beta"][1], phi_diff_se=reg["se"][1],
               phi_diff_t=reg["beta"][1] / reg["se"][1], phi_diff_hac_lag=reg["lag"],
               phi_diff_boot_lo=np.percentile(boot, 2.5), phi_diff_boot_hi=np.percentile(boot, 97.5))
    for name, x in (("zero_rate", a), ("hiking", b)):
        bm = bootstrap_distribution(np.mean, [x], B=B)
        out[f"phi_{name}_boot_lo"], out[f"phi_{name}_boot_hi"] = np.percentile(bm, [2.5, 97.5])
    for col in ("C_skew", "FD"):
        out[f"{col}_zero_rate_bp"] = 1e4 * s.loc[~D, col].mean()
        out[f"{col}_hiking_bp"] = 1e4 * s.loc[D, col].mean()

    # E2, in sample
    rp = ret.dropna(subset=["phi"])
    y, x = rp.U.to_numpy(), rp.phi.to_numpy()
    reg = ols_hac(y, np.column_stack([np.ones(len(y)), x]))
    e2 = stambaugh_bootstrap(y, x, B)
    out.update(E2_n=len(y), E2_b=reg["beta"][1], E2_se_hac=reg["se"][1], E2_t_hac=reg["beta"][1] / reg["se"][1],
               E2_bias_bootstrap=e2["bias_bootstrap"], E2_bias_analytic=e2["bias_analytic"],
               E2_b_bias_corrected=e2["b_bias_corrected"], E2_p_one_sided=e2["p_one_sided_b_gt_0"],
               E2_rho_phi=e2["rho_predictor"], E2_corr_u_v=e2["corr_u_v"])

    # E3
    for name in ("U", "H", "c3"):
        r = mean_and_se(ret[name].to_numpy())
        out[f"{name}_bp"], out[f"{name}_se_bp"], out[f"{name}_t"] = 1e4 * r["mean"], 1e4 * r["se"], r["mean"] / r["se"]
    out["theta_UB"] = 1 - ret.H.mean() / ret.U.mean()
    out["theta_UB_meaningful"] = abs(out["U_t"]) >= 0.5
    kind, lo, hi = theta_confidence_set(ret.H.to_numpy(), ret.U.to_numpy())
    out.update(theta_UB_set=kind, theta_UB_set_lo=lo, theta_UB_set_hi=hi)
    return out


def one_month_check(d: Path, spots) -> str:
    """Run the loop at one month over the primary month-ends; it must reproduce E1 and Stage 4."""
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    inputs = load("smile_inputs.csv")
    dates = [t for t in pd.DatetimeIndex(sorted(inputs.date.unique())) if t >= PRIMARY_START]
    q = quarter_rows(dates, inputs, load("smile_panel.csv"), spots, months=1).set_index("date")
    e1 = load("e1_series.csv")
    e1 = e1[(e1["sample"] == "primary") & (e1.reading == "market") & (e1.delta == DELTA) & (e1.status == "ok")].set_index("date")
    s4 = load("stage4_months.csv")
    s4 = s4[(s4["sample"] == "primary") & (s4.status == "ok")].set_index("date")
    if not (len(q.dropna(subset=["phi"])) == len(e1) and int(q.returns_realised.eq(True).sum()) == len(s4)):
        raise SystemExit("one-month check: sample sizes differ from E1 and Stage 4")
    gaps = {"phi": (q.loc[e1.index, "phi"] - e1.phi).abs().max(), "HML^U": (q.loc[s4.index, "U"] - s4.U).abs().max(),
            "HML^H": (q.loc[s4.index, "H"] - s4.H10).abs().max(), "skew term": (q.loc[s4.index, "c3"] - s4.c3).abs().max()}
    for k, v in gaps.items():
        if not v < 1e-12:
            raise SystemExit(f"one-month check: {k} differs from E1 and Stage 4 by {v}")
    return f"At one month the loop reproduces E1 ({len(e1)} month-ends) and Stage 4 ({len(s4)} returns) within 1e-12."


def quote_diagnostics(raw, extra, inputs_file, calib_month_ends, primary_dates):
    """Completeness of the three-month inputs; checks the rebuilt inputs against the calibration file."""
    st_cols = [f"status_{f}" for f in QUOTE_FIELDS]
    cal = build_month_end_inputs(raw, calib_month_ends, tenor="3M", extra_roots=extra, quote_status=True)
    if cal.drop(columns=st_cols).to_csv(index=False) != inputs_file.read_text():
        raise SystemExit("rebuilt three-month inputs differ from the calibration file")
    counts = {}
    for label, frame in (("calibration month-ends", cal), ("primary quarters", cal[cal.date.isin(primary_dates)])):
        tab = pd.DataFrame({f: frame[f"status_{f}"].value_counts() for f in QUOTE_FIELDS}).reindex(
            ["observed", "substituted", "missing"]).fillna(0).astype(int)
        tab["all fields"] = tab.sum(axis=1)
        counts[label] = (len(frame), tab)
    full = build_month_end_inputs(raw, ny_month_ends("1995-01-01", calib_month_ends.max()), tenor="3M", extra_roots=extra)
    n_complete = full.groupby("date").complete.sum()
    all_nine = n_complete[n_complete == len(G10)]
    first = all_nine.index.min() if len(all_nine) else pd.NaT
    later = n_complete[n_complete.index >= first] if len(all_nine) else n_complete.iloc[:0]
    return {"first_all_nine_complete": first, "month_ends_after_first_not_all_nine": int((later < len(G10)).sum()),
            "month_ends_from_first": len(later)}, counts


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23", help="retrieval of spot, volatilities and the calibration")
    p.add_argument("--extra-raw-root", action="append", default=None,
                   help="raw directory of the three-month forwards and rates (default: the 24 September retrieval)")
    p.add_argument("--out-date", default="2026-09-24")
    p.add_argument("--bootstrap", type=int, default=1999)
    args = p.parse_args()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    raw = ROOT / "data" / "private" / "lseg" / args.retrieval_date / "raw"
    extra = [Path(r) if Path(r).is_absolute() else ROOT / r for r in args.extra_raw_root or [EXTRA_RAW_ROOT]]
    out_dir = ROOT / "data" / "private" / "results" / args.out_date
    out_dir.mkdir(parents=True, exist_ok=True)
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    inputs, panel = load("smile_inputs_3m.csv"), load("smile_panel_3m.csv")
    spots = {c: daily_spot_mid(raw, c) for c in G10}
    calib_dates = pd.DatetimeIndex(sorted(inputs.date.unique()))
    check = one_month_check(d, spots)

    series, results = [], {}
    for phase, months in PHASES.items():
        dates = [t for t in calib_dates if t.month in months]
        if phase == "Mar/Jun/Sep/Dec" and dates[0] != PRIMARY_FIRST:
            raise SystemExit(f"primary phase starts {dates[0]:%Y-%m-%d}, not {PRIMARY_FIRST:%Y-%m-%d}")
        q = quarter_rows(dates, inputs, panel, spots)
        q["phase"] = phase
        series.append(q)
        results[phase] = summarise(q, args.bootstrap)
        print(f"done: {phase}", flush=True)
    primary_dates = series[0].date
    data_check, counts = quote_diagnostics(raw, extra, d / "smile_inputs_3m.csv", calib_dates, primary_dates)

    pd.concat(series, ignore_index=True).to_csv(out_dir / "tenor3m_quarters.csv", index=False)
    table = pd.DataFrame(results).T
    table.to_csv(out_dir / "tenor3m_summary.csv")
    fmt = lambda v: f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
    lines = ["# Three-month tenor, quarterly rebalancing (restricted)", "",
             f"Primary phase Mar/Jun/Sep/Dec; the other phases are supplementary. Bootstrap draws: {args.bootstrap}.",
             check,
             "Returns, C_skew, FD and the skew term are in basis points per quarter.", "", "```",
             table.T.to_string(float_format=fmt), "```", "", "## Three-month inputs", "",
             *[f"- {k}: {fmt(v)}" for k, v in data_check.items()], ""]
    for label, (n, tab) in counts.items():
        lines += [f"Quote status, {label} ({n} currency-dates):", "", "```", tab.to_string(), "```", ""]
    (out_dir / "tenor3m_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
