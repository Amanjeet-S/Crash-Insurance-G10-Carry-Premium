"""Option-implied variance and skewness intervals (R1, R3, R4) and the secondary E2 regressions.

Research design, sections 3, 6 and 7. For every currency-month of the primary
sample (from 2013-05-31) with a market-reading one-month SABR smile (status ok
or not_converged), ``qef.fx.moments`` gives intervals for the raw moments of
y = ln(X_T/F_X) under Q^USD, with P and C observed between the smile 10Δ put
and call strikes in the pair's convention, and the exact ranges of the variance
and skewness over the box of those intervals. Two exponent settings are
reported:

- (i) boundary elasticities of the smile at the 10Δ strikes. This is the
  assumption that the tails beyond the last quotes are no heavier than the
  smile implies at those quotes. As a diagnostic, the script records whether
  the calibrated smile's own elasticities satisfy it up to two ATM standard
  deviations beyond each boundary strike;
- (ii) heavy tails, γ = η = 2.

Portfolio predictors (E2, secondary). At each E1 month-end (primary sample,
market reading, 10Δ, status ok) the six legs give: variance, the mean of v over
the legs; oriented skewness, the mean of +skewness over long legs and
−skewness over short legs, so that a short leg's interval is
[−skew_hi, −skew_lo]. Interval endpoints are averaged endpoint by endpoint.

Regressions. Next-month HML^U (stage4_months.csv, primary, status ok; the
return from t to the expiry of the option dated t, aligned by date as in
estimate_stage4) on each predictor at its lower endpoint, midpoint and upper
endpoint under setting (i), with Newey–West t statistics (qef.stats.hac) and
the Stambaugh residual bootstrap of estimate_stage4 (one-sided p for b > 0).
The bootstrap p for b < 0 applies the same function to the negated predictor,
which negates every bootstrap slope under the same seed. The one-sided
alternatives are b > 0 for variance and b < 0 for oriented skewness (more
negative oriented skewness is more crash risk borne by the carry portfolio),
the counterparts of the b > 0 alternative used for φ and the oriented risk
reversal. Outputs are LSEG-derived and are written to data/private/results/.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from estimate_stage4 import stambaugh_bootstrap  # noqa: E402

from qef.fx.conventions import G10  # noqa: E402
from qef.fx.moments import HEAVY_TAIL_EXPONENT, implied_moment_intervals, smile_distribution  # noqa: E402
from qef.fx.sabr import sabr_vol  # noqa: E402
from qef.stats.hac import ols_hac  # noqa: E402

PRIMARY_START = pd.Timestamp("2013-05-31")  # private audit, primary-sample rule
DELTA = 0.10
SETTINGS = ("i", "ii")
PREDICTORS = (("var", "> 0"), ("oskew", "< 0"))


def smile_tail_consistent(F, tau, vol, K_min, K_max, gamma, eta, sigma_atm, span=2.0, n=8):
    """Whether the smile's own elasticities stay at or above γ below K_min and η above K_max.

    Checked on n strikes spaced evenly in log strike up to ``span`` ATM standard
    deviations beyond each boundary strike. Also returns whether the smile
    density g is positive at all of those strikes.
    """
    steps = np.linspace(0.0, span * sigma_atm * np.sqrt(tau), n + 1)[1:]
    lower = [(K, smile_distribution(K, F, tau, vol)) for K in K_min * np.exp(-steps)]
    upper = [(K, smile_distribution(K, F, tau, vol)) for K in K_max * np.exp(steps)]
    ok_lo = min(K * d["g"] / d["G"] for K, d in lower) >= gamma * (1 - 1e-9)
    ok_hi = min(K * d["g"] / d["Gbar"] for K, d in upper) >= eta * (1 - 1e-9)
    return bool(ok_lo), bool(ok_hi), bool(min(d["g"] for _, d in lower + upper) > 0)


def currency_months(inputs: pd.DataFrame, panel: pd.DataFrame) -> pd.DataFrame:
    fits = panel[(panel.reading == "market") & panel.status.isin(["ok", "not_converged"]) & (panel.date >= PRIMARY_START)]
    inp = inputs.set_index(["currency", "date"])
    rows = []
    for _, f in fits.sort_values(["date", "currency"]).iterrows():
        c, t = f.currency, f.date
        r, conv = inp.loc[(c, t)], G10[f.currency]
        vol = lambda K, r=r, f=f: sabr_vol(K, r.F, r.tau, f.alpha, f.rho, f.nu, 1.0)
        base = {"date": t, "currency": c, "pair": conv.pair, "usd_base": conv.usd_base, "smile_status": f.status,
                "F": r.F, "tau": r.tau, "atm": r.atm}
        try:
            res = implied_moment_intervals(r.F, r.tau, vol, conv.usd_base, conv.delta, r.df_base, DELTA)
        except ValueError:
            rows.append({**base, "status": "inversion_failed"})
            continue
        if res["status_i"] != "exponent_invalid":
            res["smile_tail_ok_lower_i"], res["smile_tail_ok_upper_i"], res["smile_tail_density_positive"] = smile_tail_consistent(
                r.F, r.tau, vol, res["K_min"], res["K_max"], res["gamma_i"], res["eta_i"], r.atm)
        rows.append({**base, "status": "ok", **res})
    return pd.DataFrame(rows)


def portfolio_predictors(e1: pd.DataFrame, mom: pd.DataFrame) -> pd.DataFrame:
    m = mom.set_index(["currency", "date"])
    rows = []
    for _, e in e1.iterrows():
        legs = [(c, 1.0) for c in e.longs.split()] + [(c, -1.0) for c in e.shorts.split()]
        row = {"date": e.date, "n_legs": len(legs)}
        if not all((c, e.date) in m.index and m.loc[(c, e.date), "status"] == "ok" for c, _ in legs):
            rows.append({**row, "status": "missing_leg"})
            continue
        L = [(m.loc[(c, e.date)], o) for c, o in legs]
        for s in SETTINGS:
            if any(x[f"status_{s}"] != "ok" for x, _ in L):
                row[f"status_{s}"] = "leg_not_ok"
                continue
            row[f"var_lo_{s}"] = np.mean([x[f"var_lo_{s}"] for x, _ in L])
            row[f"var_hi_{s}"] = np.mean([x[f"var_hi_{s}"] for x, _ in L])
            row[f"oskew_lo_{s}"] = np.mean([x[f"skew_lo_{s}"] if o > 0 else -x[f"skew_hi_{s}"] for x, o in L])
            row[f"oskew_hi_{s}"] = np.mean([x[f"skew_hi_{s}"] if o > 0 else -x[f"skew_lo_{s}"] for x, o in L])
            for p in ("var", "oskew"):
                row[f"{p}_mid_{s}"] = 0.5 * (row[f"{p}_lo_{s}"] + row[f"{p}_hi_{s}"])
            row[f"status_{s}"] = "ok"
        rows.append({**row, "status": "ok"})
    return pd.DataFrame(rows)


def e2_secondary(pm: pd.DataFrame, B: int) -> pd.DataFrame:
    y = pm.U.to_numpy()
    rows = []
    for pred, alt in PREDICTORS:
        for end in ("lo", "mid", "hi"):
            x = pm[f"{pred}_{end}_i"].to_numpy()
            reg = ols_hac(y, np.column_stack([np.ones(len(x)), x]))
            up, down = stambaugh_bootstrap(y, x, B), stambaugh_bootstrap(y, -x, B)
            p_alt = up["p_one_sided_b_gt_0"] if alt == "> 0" else down["p_one_sided_b_gt_0"]
            t = reg["beta"][1] / reg["se"][1]
            rows.append({"predictor": pred, "endpoint": end, "n": reg["n"], "b": reg["beta"][1], "se_hac": reg["se"][1],
                         "t_hac": t, "hac_lag": reg["lag"], "b_bias_corrected": up["b_bias_corrected"],
                         "rho_predictor": up["rho_predictor"], "p_boot_b_gt_0": up["p_one_sided_b_gt_0"],
                         "p_boot_b_lt_0": down["p_one_sided_b_gt_0"], "alternative": f"b {alt}",
                         "sig5_hac_two_sided": bool(abs(t) >= 1.96), "sig5_boot_one_sided": bool(p_alt < 0.05)})
    return pd.DataFrame(rows)


def agreement(e2: pd.DataFrame) -> list[str]:
    out = []
    for pred, g in e2.groupby("predictor", sort=False):
        same = lambda col: g[col].nunique() == 1
        out.append(f"- {pred}: sign {'agrees' if np.sign(g.b).nunique() == 1 else 'differs'} across endpoints "
                   f"({', '.join('+' if b > 0 else '-' for b in g.b)}); Newey-West 5% two-sided significance "
                   f"{'agrees' if same('sig5_hac_two_sided') else 'differs'} ({', '.join(map(str, g.sig5_hac_two_sided))}); "
                   f"bootstrap 5% one-sided ({g.alternative.iloc[0]}) {'agrees' if same('sig5_boot_one_sided') else 'differs'} "
                   f"({', '.join(map(str, g.sig5_boot_one_sided))})")
    return out


def excludes_zero(lo, hi):
    return (lo > 0) | (hi < 0)


def summary_lines(mom, pm_all, pm, e2, B):
    ok = mom[mom.status == "ok"]
    q = lambda s, ps=(0.05, 0.5, 0.95): ", ".join(f"{s.quantile(p):.4g}" for p in ps)
    fmt = lambda v: f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
    lines = ["# Option-implied variance and skewness intervals; secondary E2 (restricted)", "",
             "Aggregates only. Moments of y = ln(X_T/F_X) under Q^USD (R1, R3, R4); prices between the smile 10-delta "
             "strikes; setting (i) boundary elasticities (tails beyond the last quotes no heavier than the smile "
             f"implies there), setting (ii) gamma = eta = {HEAVY_TAIL_EXPONENT:g}.", "",
             "## Currency-months", "",
             f"- Attempted (primary, market reading, smile ok or not_converged): {len(mom)}; ok: {len(ok)}; "
             f"inversion failed: {int((mom.status == 'inversion_failed').sum())}; non-converged smiles used: "
             f"{int((ok.smile_status == 'not_converged').sum())}",
             f"- Status by setting: " + "; ".join(f"({s}) " + ", ".join(f"{k} {v}" for k, v in ok[f'status_{s}'].value_counts().items())
                                                  for s in SETTINGS),
             f"- Simpson halving check: converged {int(ok.simpson_converged.sum())}/{len(ok)}; subintervals per side "
             f"median {ok.n_simpson.median():.0f}, max {ok.n_simpson.max():.0f}; largest relative change {ok.simpson_rel_change.max():.2e}",
             "- Tail quadrature: largest relative error estimate "
             + ", ".join(f"({s}) {ok[f'tail_quad_err_{s}'].max():.2e}" for s in SETTINGS)
             + "; largest relative difference from the R4 closed forms (K^-2 and logarithmic weights) "
             + ", ".join(f"({s}) {ok[f'closed_form_rel_diff_{s}'].max():.2e}" for s in SETTINGS),
             "- Range search (vertices and edge critical points, then multi-start L-BFGS-B): largest gain of the search "
             + ", ".join(f"({s}) {ok[f'range_gap_{s}'].max():.2e}" for s in SETTINGS),
             f"- Setting (i) exponents, 5%/50%/95% quantiles: gamma {q(ok.gamma_i)}; eta {q(ok.eta_i)}",
             f"- Setting (i) diagnostic, calibrated smile's own elasticities up to 2 ATM s.d. beyond the boundary "
             f"at least the boundary value: lower tail {ok.smile_tail_ok_lower_i.mean():.1%}, upper tail {ok.smile_tail_ok_upper_i.mean():.1%} "
             f"of currency-months (smile density positive on the checked strikes: {ok.smile_tail_density_positive.mean():.1%})",
             ""]
    lines += ["Interval summaries (medians over currency-months; annualised volatility bounds are sqrt(v/tau)):", "",
              "| setting | pairs | n | vol lower | vol upper | median v width | median skew lower | median skew upper | median skew width | skew interval excludes 0 |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for s in SETTINGS:
        for label, sub in (("all", ok), ("EURUSD-type", ok[~ok.usd_base]), ("USD-base", ok[ok.usd_base])):
            sub = sub[sub[f"status_{s}"] == "ok"]
            lines.append(f"| ({s}) | {label} | {len(sub)} | {np.sqrt(sub[f'var_lo_{s}'] / sub.tau).median():.4f} | "
                         f"{np.sqrt(sub[f'var_hi_{s}'] / sub.tau).median():.4f} | {(sub[f'var_hi_{s}'] - sub[f'var_lo_{s}']).median():.3e} | "
                         f"{sub[f'skew_lo_{s}'].median():.3g} | {sub[f'skew_hi_{s}'].median():.3g} | "
                         f"{(sub[f'skew_hi_{s}'] - sub[f'skew_lo_{s}']).median():.3g} | "
                         f"{excludes_zero(sub[f'skew_lo_{s}'], sub[f'skew_hi_{s}']).mean():.1%} |")
    lines += ["", "## Portfolio predictors (E1 legs, primary sample)", "",
              f"- E1 months: {len(pm_all)}; with all six legs ok: {int((pm_all.status == 'ok').sum())}; "
              f"used in E2 (matched to HML^U, status ok): {len(pm)}", "",
              "| setting | n | median var width | median oriented skew lower | median oriented skew upper | median width | oriented skew interval excludes 0 | share with upper < 0 |",
              "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    ready = pm_all[pm_all.status == "ok"]
    for s in SETTINGS:
        sub = ready[ready[f"status_{s}"] == "ok"]
        lo, hi = sub[f"oskew_lo_{s}"], sub[f"oskew_hi_{s}"]
        lines.append(f"| ({s}) | {len(sub)} | {(sub[f'var_hi_{s}'] - sub[f'var_lo_{s}']).median():.3e} | {lo.median():.3g} | "
                     f"{hi.median():.3g} | {(hi - lo).median():.3g} | {excludes_zero(lo, hi).mean():.1%} | {(hi < 0).mean():.1%} |")
    lines += ["", f"## E2 secondary: HML^U(t+1) on the setting (i) predictor at each endpoint (bootstrap draws: {B})", "",
              "Variance is that of the one-month log return (not annualised). Newey-West t with the automatic bandwidth; "
              "bootstrap p values from the Stambaugh residual bootstrap under the null of no predictability.", "", "```",
              e2.to_string(index=False, float_format=fmt), "```", "", "Agreement across endpoints:", "", *agreement(e2), ""]
    return lines


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--bootstrap", type=int, default=1999)
    args = p.parse_args()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])

    mom = currency_months(load("smile_inputs.csv"), load("smile_panel.csv"))
    mom.to_csv(d / "moments_1m.csv", index=False)

    e1 = load("e1_series.csv")
    e1 = e1[(e1["sample"] == "primary") & (e1.reading == "market") & (e1.delta == DELTA) & (e1.status == "ok")].sort_values("date")
    pm_all = portfolio_predictors(e1, mom)
    s4 = load("stage4_months.csv")
    s4 = s4[(s4["sample"] == "primary") & (s4.status == "ok")][["date", "U"]]
    pm = pm_all[(pm_all.status == "ok") & (pm_all.status_i == "ok")].merge(s4, on="date", how="inner").sort_values("date")
    e2 = e2_secondary(pm, args.bootstrap)

    lines = summary_lines(mom, pm_all, pm, e2, args.bootstrap)
    (d / "moments_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
