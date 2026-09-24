"""Estimand E1 (and E4): the ex-ante skew price of crash insurance per unit of carry.

Research design, section 6. At each month-end the available currencies are
ranked by forward discount; the portfolio is long the top n and short the
bottom n (n = 3 in the primary sample, 2 in the extended sample). A currency
without the calibration set is excluded that month; a month is used only if at
least 2n currencies remain. Non-converged smiles are used and counted. For
each leg the protective option at the smile delta strike is priced from the
calibrated smile (V_smile) and at the ATM volatility (V_flat). Then

    C_skew = (1/n) Σ_legs (V_smile − V_flat),   FD = (1/n)(Σ_long fd − Σ_short fd),
    φ = C_skew / FD.

E1 uses the market-strangle reading; E4 repeats it under the smile-strangle
reading. 10-delta is primary and 25-delta a comparison. The ATM comparison is
identically zero because the smile passes through the ATM volatility.

Pre-registered inference: regime means (zero-rate: primary start to December
2021; hiking: January 2022 onwards) with Newey–West standard errors and a
stationary-bootstrap percentile interval for the difference. Supplementary
analyses added after the first results (research log, 24 September 2026) are
reported separately and labelled as such.

Splice: leg-level comparison of Fenics and composite skew costs on every
currency-month where both are available. Outputs are LSEG-derived and are
written to data/private/results/.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from qef.fx.conventions import G10
from qef.fx.crash import forward_discount, leg_skew_cost, rank_legs
from qef.fx.sabr import sabr_vol
from qef.stats.bootstrap import bootstrap_distribution, optimal_block_length
from qef.stats.hac import mean_and_se, ols_hac

ROOT = Path(__file__).resolve().parents[1]
REGIME_BREAK = pd.Timestamp("2022-01-01")
PRIMARY_START = pd.Timestamp("2013-05-31")  # private audit, primary-sample rule


def _vol(r, f):
    return lambda K: sabr_vol(K, r.F, r.tau, f.alpha, f.rho, f.nu, 1.0)


def e1_series(inputs, panel, reading, delta, n_legs):
    fits = panel[(panel.reading == reading) & panel.status.isin(["ok", "not_converged"])].set_index(["currency", "date"])
    rows = []
    for date, month in inputs[inputs.complete_calib].groupby("date"):
        avail = {r.currency: (r, fits.loc[(r.currency, date)]) for _, r in month.iterrows() if (r.currency, date) in fits.index}
        if len(avail) < 2 * n_legs:
            rows.append({"date": date, "status": "too_few_currencies", "n_ccy": len(avail)})
            continue
        fd = {c: forward_discount(r.S, r.F, G10[c].usd_base) for c, (r, _) in avail.items()}
        longs, shorts = rank_legs(fd, n_legs)
        c_skew, failed = 0.0, None
        for c, is_long in [(c, True) for c in longs] + [(c, False) for c in shorts]:
            r, f = avail[c]
            try:
                _, v_s, v_f = leg_skew_cost(_vol(r, f), r.F, r.tau, r.df_base, c, r.atm, is_long, delta)
            except ValueError:
                failed = c
                break
            c_skew += (v_s - v_f) / n_legs
        if failed:
            rows.append({"date": date, "status": f"inversion_failed:{failed}", "n_ccy": len(avail)})
            continue
        FD = (sum(fd[c] for c in longs) - sum(fd[c] for c in shorts)) / n_legs
        n_nc = sum(avail[c][1].status == "not_converged" for c in longs + shorts)
        rows.append({"date": date, "status": "ok", "n_ccy": len(avail), "n_not_converged": int(n_nc),
                     "longs": " ".join(longs), "shorts": " ".join(shorts),
                     "C_skew": c_skew, "FD": FD, "phi": c_skew / FD if FD > 0 else np.nan})
    out = pd.DataFrame(rows)
    out["reading"], out["delta"] = reading, delta
    return out


def regime_split(s, col):
    s = s[(s.status == "ok")].dropna(subset=[col])
    return s, s.loc[s.date < REGIME_BREAK, col].to_numpy(), s.loc[s.date >= REGIME_BREAK, col].to_numpy()


def regime_summary(s, col, B):
    s, a, b = regime_split(s, col)
    ra, rb = mean_and_se(a), mean_and_se(b)
    X = np.column_stack([np.ones(len(s)), (s.date >= REGIME_BREAK).astype(float)])
    reg = ols_hac(s[col].to_numpy(), X)
    boot = bootstrap_distribution(lambda x, y: y.mean() - x.mean(), [a, b], B=B)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return {"variable": col, "n_zero_rate": ra["n"], "mean_zero_rate": ra["mean"], "se_zero_rate": ra["se"],
            "n_hiking": rb["n"], "mean_hiking": rb["mean"], "se_hiking": rb["se"], "difference": reg["beta"][1],
            "se_difference_hac": reg["se"][1], "hac_lag": reg["lag"], "boot_lo": lo, "boot_hi": hi}


def _ar1(x):
    x = x - x.mean()
    return float(x[1:] @ x[:-1] / (x @ x))


def supplementary(s, B):
    """Post-hoc analyses (research log, 24 September 2026); not pre-registered."""
    s, a, b = regime_split(s, "phi")
    diff = b.mean() - a.mean()
    out = {}
    X = np.column_stack([np.ones(len(s)), (s.date >= REGIME_BREAK).astype(float)])
    y = s["phi"].to_numpy()
    e = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    XtX_inv = np.linalg.inv(X.T @ X)
    for L in (0, 3, 6, 9, 12, 15, 20):
        scores = X * e[:, None]
        S = scores.T @ scores / len(y)
        for j in range(1, L + 1):
            G = scores[j:].T @ scores[: len(y) - j] / len(y)
            S += (1 - j / (L + 1)) * (G + G.T)
        out[f"se_difference_lag_{L}"] = float(np.sqrt((len(y) * XtX_inv @ S @ XtX_inv)[1, 1]))
    se_ar1 = []
    for name, x in (("zero_rate", a), ("hiking", b)):
        rho = _ar1(x)
        out[f"ar1_{name}"] = rho
        out[f"effective_n_{name}"] = len(x) * (1 - rho) / (1 + rho)
        v = np.var(x) * (1 + rho) / (1 - rho)  # long-run variance of an AR(1)
        se_ar1.append(np.sqrt(v / len(x)))
        out[f"median_{name}"] = float(np.median(x))
        k = int(0.1 * len(x))
        xs = np.sort(x)
        out[f"trimmed10_{name}"] = float(xs[k: len(x) - k].mean())
        out[f"mean_log_phi_{name}"] = float(np.log(x).mean()) if np.all(x > 0) else np.nan
        out[f"block_length_{name}"] = optimal_block_length(x)
    out["se_difference_ar1_plugin"] = float(np.hypot(*se_ar1))
    out["t_difference_ar1_plugin"] = diff / out["se_difference_ar1_plugin"]
    out["log_phi_difference"] = out["mean_log_phi_hiking"] - out["mean_log_phi_zero_rate"]
    zr, hk = s[s.date < REGIME_BREAK], s[s.date >= REGIME_BREAK]
    out["ratio_of_means_zero_rate"] = zr.C_skew.mean() / zr.FD.mean()
    out["ratio_of_means_hiking"] = hk.C_skew.mean() / hk.FD.mean()
    keep = zr[(zr.date < "2020-01-01")]
    out["difference_excluding_2020_2021"] = b.mean() - keep.phi.mean()
    out["n_zero_rate_excluding_2020_2021"] = len(keep)
    for k in (1, 3, 5, 10):
        out[f"difference_drop_top{k}_zero_rate"] = b.mean() - np.sort(a)[: len(a) - k].mean()
    boot = bootstrap_distribution(lambda x, y: y.mean() - x.mean(), [a, b], B=B)
    out["bootstrap_p_two_sided_percentile"] = float(2 * min((boot >= 0).mean(), (boot <= 0).mean()))
    return out


def leg_splice(comp_in, comp_panel, fn_in, fn_panel, delta=0.10):
    fc = comp_panel[comp_panel.reading == "market"].set_index(["currency", "date"])
    ff = fn_panel[fn_panel.reading == "market"].set_index(["currency", "date"])
    ci = comp_in[comp_in.complete_calib].set_index(["currency", "date"])
    fi = fn_in[fn_in.complete_calib].set_index(["currency", "date"])
    common = ci.index.intersection(fi.index).intersection(fc.index).intersection(ff.index)
    rows = []
    for key in common:
        c, date = key
        for long_leg in (True, False):
            vals = []
            for inp, fit in ((ci, fc), (fi, ff)):
                r, f = inp.loc[key], fit.loc[key]
                try:
                    _, v_s, v_f = leg_skew_cost(_vol(r, f), r.F, r.tau, r.df_base, c, r.atm, long_leg, delta)
                    vals.append(v_s - v_f)
                except ValueError:
                    vals.append(np.nan)
            rows.append({"currency": c, "date": date, "side": "long" if long_leg else "short",
                         "composite": vals[0], "fenics": vals[1]})
    d = pd.DataFrame(rows).dropna()
    d["window"] = np.where(d.date < PRIMARY_START, "pre_primary", "primary")
    g = d.assign(diff=d.fenics - d.composite).groupby("window")
    return pd.DataFrame({"n": g.size(), "mean_diff": g["diff"].mean(), "sd_diff": g["diff"].std(),
                         "corr": g.apply(lambda x: np.corrcoef(x.fenics, x.composite)[0, 1])})


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--bootstrap", type=int, default=9999)
    args = p.parse_args()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    comp_in, comp_panel = load("smile_inputs.csv"), load("smile_panel.csv")
    fn_in, fn_panel = load("smile_inputs_fn.csv"), load("smile_panel_fn.csv")
    prim_in = comp_in[comp_in.date >= PRIMARY_START]

    series, summaries = [], []
    for reading in ("market", "smile"):
        for delta in (0.10, 0.25):
            s = e1_series(prim_in, comp_panel, reading, delta, n_legs=3)
            s["sample"] = "primary"
            series.append(s)
            for col in ("phi", "C_skew", "FD"):
                summaries.append({"sample": "primary", "reading": reading, "delta": delta, **regime_summary(s, col, args.bootstrap)})
    primary = series[0]
    supp = supplementary(primary, args.bootstrap)
    splice = leg_splice(comp_in, comp_panel, fn_in, fn_panel)
    ext = e1_series(fn_in[fn_in.date < PRIMARY_START], fn_panel, "market", 0.10, n_legs=2)
    ext["sample"] = "extended"
    series.append(ext)
    ext_ok = ext[ext.status == "ok"]

    out = pd.concat(series, ignore_index=True)
    out.to_csv(d / "e1_series.csv", index=False)
    summ = pd.DataFrame(summaries)
    summ.to_csv(d / "e1_summary.csv", index=False)
    fmt = lambda v: f"{v:.6g}"
    status_counts = out.groupby(["sample", "reading", "delta"])["status"].value_counts().unstack(fill_value=0)
    lines = ["# E1 and E4 (restricted)", "",
             "## Month status", "", "```", status_counts.to_string(), "```", "",
             "## Pre-registered: regime means and differences", "", "```", summ.to_string(index=False, float_format=fmt), "```", "",
             "## Supplementary (post hoc): phi, primary sample, market reading, 10-delta", "", "```",
             pd.Series(supp).to_string(float_format=fmt), "```", "",
             "## Leg-level splice: Fenics minus composite skew cost (market reading, 10-delta)", "", "```",
             splice.to_string(float_format=fmt), "```", "",
             "## Extended sample (Fenics, 2L/2S, market reading, 10-delta)", "", "```",
             ext_ok[["C_skew", "FD", "phi"]].describe().to_string(float_format=fmt), "```",
             f"Months used: {len(ext_ok)}; first {ext_ok.date.min():%Y-%m}, last {ext_ok.date.max():%Y-%m}", ""]
    (d / "e1_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
