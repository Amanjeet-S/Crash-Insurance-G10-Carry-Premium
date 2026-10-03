"""Design items completed after the main results (research log, 2 October 2026).

A read of the whole project against the research design on 2 October 2026
found three items that the design or the log had committed to and that no
earlier script computed. This script computes them, under rules fixed here
before any of them was run:

1. E3 by regime. The design's intended contribution (section 2) decomposes
   the cost of crash insurance across the two rate regimes, but Stage 4
   reported the three E3 terms for the pooled sample only. For the primary
   sample (159 return months) this reports the regime means of HML^U, HML^H
   (10Δ) and the terms (i) payoff, (ii) volatility level and (iii) skew, each
   with its Newey–West standard error, and the difference (hiking minus
   zero-rate) with the Newey–West error of the coefficient on a regime dummy
   and the stationary-bootstrap 95% percentile interval with 9,999 draws,
   resampling each regime separately, exactly as for E1. Regimes are assigned
   by the month-end at which the position is opened, as in E1.

2. Stale butterflies in the extended sample. The research log of
   23 September 2026 states that the stale-butterfly robustness run applies to
   the extended sample, but robustness.py flags primary-sample month-ends only.
   The same audit rule is applied to the Fenics 25Δ butterflies of the
   extended sample (a currency-month is dropped when the selected quote lies
   in a run of five or more identical two-sided business-day quotes), and the
   extended E1 and E3 quantities are recomputed with those currency-months
   treated as missing, together with the Clark–West statistic of E2, whose
   training window starts in the extended sample. The unflagged run must
   reproduce the stored extended series and the Clark–West statistic, or the
   script stops.

3. Bootstrap draws. The design fixes 9,999 stationary-bootstrap draws
   (section 7); the robustness grid, the three-month tenor and the secondary
   moment predictors were first run with 1,999. Their scripts now use 9,999.
   Here each is computed with both numbers of draws, under the same seeds, and
   a conclusion changes if a bootstrap interval's exclusion of zero, or a
   one-sided bootstrap test at 5%, differs between the two.

Outputs, all aggregate, are written to data/private/results/2026-09-23/:
outstanding_items_summary.md and the tables behind it.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import estimate_3m  # noqa: E402
import estimate_moments  # noqa: E402
import robustness  # noqa: E402
from estimate_e1 import PRIMARY_START, REGIME_BREAK, e1_series  # noqa: E402
from estimate_stage4 import clark_west, daily_spot_mid, month_rows  # noqa: E402

from qef.data.panel import on_business_days, read_raw, sample_month_ends, stale_mask, two_sided  # noqa: E402
from qef.fx.conventions import G10  # noqa: E402
from qef.stats.bootstrap import bootstrap_distribution  # noqa: E402
from qef.stats.hac import mean_and_se, ols_hac  # noqa: E402

DRAWS = 9999
E3_COLUMNS = {"U": "HML^U", "H10": "HML^H, 10-delta", "c1": "(i) payoff", "c2": "(ii) volatility level", "c3": "(iii) skew"}


def e3_by_regime(pm: pd.DataFrame, B: int) -> pd.DataFrame:
    D = (pm.date >= REGIME_BREAK).to_numpy()
    rows = []
    for col, label in E3_COLUMNS.items():
        x = pm[col].to_numpy()
        a, b = x[~D], x[D]
        ra, rb = mean_and_se(a), mean_and_se(b)
        reg = ols_hac(x, np.column_stack([np.ones(len(x)), D.astype(float)]))
        boot = bootstrap_distribution(lambda u, v: v.mean() - u.mean(), [a, b], B=B)
        lo, hi = np.percentile(boot, [2.5, 97.5])
        rows.append({"quantity": label, "n_zero_rate": ra["n"], "zero_rate_bp": 1e4 * ra["mean"], "zero_rate_se_bp": 1e4 * ra["se"],
                     "n_hiking": rb["n"], "hiking_bp": 1e4 * rb["mean"], "hiking_se_bp": 1e4 * rb["se"],
                     "difference_bp": 1e4 * reg["beta"][1], "difference_se_bp": 1e4 * reg["se"][1],
                     "boot_lo_bp": 1e4 * lo, "boot_hi_bp": 1e4 * hi})
    return pd.DataFrame(rows)


def fenics_stale_flags(raw: Path, month_ends) -> set:
    """(currency, month-end) whose selected Fenics 25Δ butterfly is stale under the audit rule."""
    out = set()
    for c in G10:
        f = on_business_days(read_raw(raw / "vol_bf25" / f"{c}1MBF=FN.csv"))
        mask = stale_mask(f[two_sided(f)], robustness.STALE_RUN)
        me = sample_month_ends(f, month_ends)
        me = me[me.status != "missing"]
        for m, src in me.source_date.items():
            if bool(mask.get(src, False)):
                out.add((c, m))
    return out


def extended_quantities(e1, months) -> dict:
    ok = months[months.status == "ok"]
    out = {"months_phi": int((e1.status == "ok").sum()), "months_returns": len(ok),
           "C_skew_bp": 1e4 * e1[e1.status == "ok"].C_skew.mean(), "phi": e1[e1.status == "ok"].phi.mean()}
    for col in ("U", "H10", "c1", "c2", "c3"):
        r = mean_and_se(ok[col].to_numpy())
        out[f"{col}_bp"], out[f"{col}_se_bp"] = 1e4 * r["mean"], 1e4 * r["se"]
    return out


def stale_extended(d: Path, raw: Path, spots) -> tuple[pd.DataFrame, int, dict]:
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    fn_in, fn_panel = load("smile_inputs_fn.csv"), load("smile_panel_fn.csv")
    ext_in = fn_in[fn_in.date < PRIMARY_START].copy()
    stored = load("stage4_months.csv")
    prim = stored[(stored["sample"] == "primary") & (stored.status == "ok")]
    month_ends = sorted(pd.Timestamp(x) for x in ext_in.date.unique())
    stale = fenics_stale_flags(raw, month_ends)

    results, cw = {}, {}
    for label, drop in (("as reported", set()), ("stale butterflies missing", stale)):
        inp = ext_in.copy()
        if drop:
            flagged = [(c, t) in drop for c, t in zip(inp.currency, inp.date)]
            inp.loc[flagged, "complete_calib"] = False
        e1 = e1_series(inp, fn_panel, "market", 0.10, n_legs=2)
        months = month_rows(e1[e1.status == "ok"], inp, fn_panel, spots, 2)
        results[label] = extended_quantities(e1, months)
        allm = pd.concat([months[months.status == "ok"], prim], ignore_index=True).sort_values("date")
        r = clark_west(allm.date.to_numpy(), allm.U.to_numpy(), allm.phi.to_numpy())
        cw[label] = r
        results[label].update(CW_t=r["cw_t"], CW_p=r["cw_p_one_sided"], CW_n=r["n_forecasts"])
        if not drop:
            ext_stored = stored[(stored["sample"] == "extended") & (stored.status == "ok")].set_index("date")
            mine = months[months.status == "ok"].set_index("date")
            gap = float((mine[["U", "H10", "c1", "c2", "c3"]] - ext_stored.loc[mine.index, ["U", "H10", "c1", "c2", "c3"]]).abs().max().max())
            if len(mine) != len(ext_stored) or not gap < 1e-12:
                raise SystemExit(f"unflagged extended run does not reproduce stage4_months.csv (gap {gap})")
            if abs(r["cw_t"] - 1.639) > 0.0006:
                raise SystemExit(f"unflagged run does not reproduce the Clark-West statistic: {r['cw_t']}")
    stale_in_window = sorted({m for _, m in stale})
    return pd.DataFrame(results).T, len(stale), {"months_with_flag": len(stale_in_window)}


def robustness_rerun(d: Path, raw: Path, spots, B: int, B_old: int = 1999) -> pd.DataFrame:
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    data = {"": (load("smile_inputs.csv"), load("smile_panel.csv")),
            "fn": (load("smile_inputs_fn.csv"), load("smile_panel_fn.csv"))}
    vv_panel = load("smile_panel_vv.csv")
    month_ends = sorted(pd.Timestamp(x) for x in data[""][0].date.unique() if pd.Timestamp(x) >= PRIMARY_START)
    sides, daily = robustness.quote_sides(raw, month_ends)
    stale = robustness.stale_flags(raw, month_ends)
    in_period = lambda lo, hi: {m for m in month_ends if pd.Timestamp(lo) <= m <= pd.Timestamp(hi)}
    V = robustness.Variant
    variants = [
        V("base"), V("hedge 25d", hedge="25d"), V("hedge ATM", hedge="atm"), V("smile-strangle reading", reading="smile"),
        V("Fenics quotes", contributor="fn"), V("stale butterflies missing", stale_missing=True),
        V("excluding March 2020", exclude_months=in_period("2020-02-01", "2020-03-31")),
        V("excluding CHF Dec 2014-Mar 2015", exclude_pairs={("CHF", m) for m in in_period("2014-11-01", "2015-03-31")}),
        V("implementable, k = 1", cost_k=1.0), V("implementable, k = 1.5", cost_k=1.5), V("implementable, k = 2", cost_k=2.0),
        V("dollar carry", portfolio="dollar"), V("ten currencies", portfolio="ten"), V("log returns", log_returns=True),
        V("previous-day spot for payoffs", payoff_lag=1), V("vanna-volga smiles", smile="vv"),
    ]
    stored = pd.read_csv(d / "robustness_summary.csv", index_col=0)
    rows = []
    for v in variants:
        df = robustness.run_variant(v, *data[v.contributor], spots, sides, daily, stale, vv_panel)
        s, old = robustness.summarise(df, B), robustness.summarise(df, B_old)
        ref = stored.loc[v.name]
        if abs(s["U_bp"] - ref["U_bp"]) > 1e-9 or (np.isfinite(ref.get("phi_diff", np.nan)) and abs(s["phi_diff"] - ref["phi_diff"]) > 1e-12):
            raise SystemExit(f"rerun of {v.name} does not reproduce the stored point estimates")
        row = {"variant": v.name}
        if "phi_diff_boot_lo" in s:
            row.update(boot_lo_1999=old["phi_diff_boot_lo"], boot_hi_1999=old["phi_diff_boot_hi"],
                       boot_lo_9999=s["phi_diff_boot_lo"], boot_hi_9999=s["phi_diff_boot_hi"],
                       E2_p_1999=old["E2_p"], E2_p_9999=s["E2_p"], E2_b_1999=old["E2_b_bc"], E2_b_9999=s["E2_b_bc"])
        rows.append(row)
        print(f"robustness {v.name} done", flush=True)
    return pd.DataFrame(rows)


def three_month_rerun(d: Path, raw: Path, spots, B: int, B_old: int = 1999) -> pd.DataFrame:
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    inputs, panel = load("smile_inputs_3m.csv"), load("smile_panel_3m.csv")
    stored = pd.read_csv(ROOT / "data" / "private" / "results" / "2026-09-24" / "tenor3m_summary.csv", index_col=0)
    calib_dates = pd.DatetimeIndex(sorted(inputs.date.unique()))
    rows = []
    for phase, months in estimate_3m.PHASES.items():
        dates = [t for t in calib_dates if t.month in months]
        q = estimate_3m.quarter_rows(dates, inputs, panel, spots)
        s, old = estimate_3m.summarise(q, B), estimate_3m.summarise(q, B_old)
        if abs(s["phi_diff"] - stored.loc[phase, "phi_diff"]) > 1e-12:
            raise SystemExit(f"rerun of the {phase} phase does not reproduce the stored point estimate")
        rows.append({"phase": phase, "boot_lo_1999": old["phi_diff_boot_lo"], "boot_hi_1999": old["phi_diff_boot_hi"],
                     "boot_lo_9999": s["phi_diff_boot_lo"], "boot_hi_9999": s["phi_diff_boot_hi"],
                     "E2_p_1999": old["E2_p_one_sided"], "E2_p_9999": s["E2_p_one_sided"]})
        print(f"three-month {phase} done", flush=True)
    return pd.DataFrame(rows)


def moments_rerun(d: Path, B: int) -> pd.DataFrame:
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    mom = load("moments_1m.csv")
    e1 = load("e1_series.csv")
    e1 = e1[(e1["sample"] == "primary") & (e1.reading == "market") & (e1.delta == estimate_moments.DELTA) & (e1.status == "ok")].sort_values("date")
    pm_all = estimate_moments.portfolio_predictors(e1, mom)
    s4 = load("stage4_months.csv")
    s4 = s4[(s4["sample"] == "primary") & (s4.status == "ok")][["date", "U"]]
    pm = pm_all[(pm_all.status == "ok") & (pm_all.status_i == "ok")].merge(s4, on="date", how="inner").sort_values("date")
    old, new = estimate_moments.e2_secondary(pm, 1999), estimate_moments.e2_secondary(pm, B)
    out = old[["predictor", "endpoint", "alternative"]].copy()
    for col in ("p_boot_b_gt_0", "p_boot_b_lt_0", "b_bias_corrected"):
        out[f"{col}_1999"], out[f"{col}_9999"] = old[col].to_numpy(), new[col].to_numpy()
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--bootstrap", type=int, default=DRAWS)
    args = p.parse_args()
    t0 = time.time()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    raw = ROOT / "data" / "private" / "lseg" / args.retrieval_date / "raw"
    spots = {c: daily_spot_mid(raw, c) for c in G10}

    s4 = pd.read_csv(d / "stage4_months.csv", parse_dates=["date"])
    pm = s4[(s4["sample"] == "primary") & (s4.status == "ok")].sort_values("date")
    e3r = e3_by_regime(pm, args.bootstrap)
    e3r.to_csv(d / "e3_by_regime.csv", index=False)
    print("E3 by regime done", flush=True)

    stale_tab, n_stale, stale_info = stale_extended(d, raw, spots)
    stale_tab.to_csv(d / "extended_stale_check.csv")
    print("extended stale check done", flush=True)

    rob = robustness_rerun(d, raw, spots, args.bootstrap)
    rob.to_csv(d / "robustness_draws_check.csv", index=False)
    tm = three_month_rerun(d, raw, spots, args.bootstrap)
    tm.to_csv(d / "tenor3m_draws_check.csv", index=False)
    mo = moments_rerun(d, args.bootstrap)
    mo.to_csv(d / "moments_draws_check.csv", index=False)

    excl = lambda lo, hi: bool(lo > 0 or hi < 0)
    changes = []
    for _, r in rob.dropna(subset=["boot_lo_1999"]).iterrows():
        if excl(r.boot_lo_1999, r.boot_hi_1999) != excl(r.boot_lo_9999, r.boot_hi_9999):
            changes.append(f"robustness {r.variant}: E1 bootstrap interval")
        if (r.E2_p_1999 < 0.05) != (r.E2_p_9999 < 0.05):
            changes.append(f"robustness {r.variant}: E2 one-sided test")
    for _, r in tm.iterrows():
        if excl(r.boot_lo_1999, r.boot_hi_1999) != excl(r.boot_lo_9999, r.boot_hi_9999):
            changes.append(f"three-month {r.phase}: E1 bootstrap interval")
        if (r.E2_p_1999 < 0.05) != (r.E2_p_9999 < 0.05):
            changes.append(f"three-month {r.phase}: E2 one-sided test")
    for _, r in mo.iterrows():
        col = "p_boot_b_gt_0" if r.alternative == "b > 0" else "p_boot_b_lt_0"
        if (r[f"{col}_1999"] < 0.05) != (r[f"{col}_9999"] < 0.05):
            changes.append(f"moments {r.predictor} {r.endpoint}: one-sided test")

    fmt = lambda v: f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
    lines = ["# Design items completed after the main results (restricted)", "",
             f"Bootstrap draws: {args.bootstrap}. Run time {time.time() - t0:.0f} s.", "",
             "## 1. E3 by regime (primary sample, bp per month)", "", "```", e3r.to_string(index=False, float_format=fmt), "```", "",
             "## 2. Extended sample with stale Fenics butterflies treated as missing", "",
             f"Flagged currency-months: {n_stale}, in {stale_info['months_with_flag']} month-ends.", "", "```",
             stale_tab.to_string(float_format=fmt), "```", "",
             "## 3. 1,999 against 9,999 bootstrap draws", "", "Robustness grid (E1 interval for the regime difference in phi; E2 one-sided p):", "",
             "```", rob.to_string(index=False, float_format=fmt), "```", "", "Three-month tenor:", "", "```",
             tm.to_string(index=False, float_format=fmt), "```", "", "Secondary moment predictors:", "", "```",
             mo.to_string(index=False, float_format=fmt), "```", "",
             "Conclusions that change at the 5% level: " + ("; ".join(changes) if changes else "none"), ""]
    (d / "outstanding_items_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
