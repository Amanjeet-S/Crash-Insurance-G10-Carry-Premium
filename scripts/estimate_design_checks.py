"""Four design checks that the paper lists as not reported (paper, Section 4.7).

Added on 27 September 2026, after every pre-registered result had been
recorded. Each check is named in the research design or the data plan, but
the design leaves details open; I fixed them as rules R1 to R12 below before
computing any of the results, except for the items that R10 and R11 mark as
added after a first run with 99 bootstrap draws and the amendments that follow
R12, made after an independent review. Citations are inherited from the design
and the paper and were not consulted for this script.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/estimate_design_checks.py

1. SOFR OIS check (design, section 3: USD discounting uses one-month Fed Funds
   OIS, USD1MOIS=, "with SOFR OIS as a check from its start", USDSROIS1M=).
   The USD rate r enters E1 only for the four pairs quoted against USD
   (EURUSD, GBPUSD, AUDUSD, NZDUSD). There it gives the quote-currency factor
   D_q = 1/(1 + r d/360), with d the forward's days to maturity (as in
   qef.data.smile_inputs), and the base-currency factor D_b = F D_q/S, which
   scales the pips spot delta Δ = ω D_b Φ(ω d₊) (R5(a) of the paper;
   Reiswich and Wystup, 2012). D_b therefore moves the 25Δ strikes of the
   calibration, hence the smile, and the 10Δ strike of the hedge. For the
   dollar-base pairs D_q comes from the foreign deposit rate and D_b = F D_q/S
   is the USD factor implied by the forward, so no USD rate enters. Forward
   discounts depend only on S and F, so the legs and FD are unchanged. I
   recalibrate the four affected smiles with SOFR through the call path of
   calibrate_smiles.py (its function _calibrate_currency, imported and not
   modified) over the whole history of the stored panel, and recompute E1
   with e1_series and regime_summary of estimate_e1.py (also imported).

2. Named-broker check (data plan: TIFO, EUR1MO=TIFO, EUR1MRR=TIFO and
   EUR1MBF=TIFO, "only as a check over its available period"). The pooled
   daily comparison is read from the private audit's splice.csv (check 8 of
   audit_fx_panel.py). I add the comparison at month-ends under the design's
   month-end rule, from the raw tables, and the effect on the skew cost of the
   EUR legs: a SABR smile calibrated to the TIFO calibration set, with spot,
   forward and rate from the composite inputs, and the 10Δ leg skew cost
   V_smile − V_flat (leg_skew_cost of qef.fx.crash) compared with the
   composite and Fenics ones at the same month-ends.

3. Delta-method intervals for θ_UB = 1 − μ̂_H/μ̂_U (design, section 7: "delta-
   method intervals are secondary"), for the 10Δ, 25Δ and ATM hedges of the
   primary sample (stage4_months.csv). With g = ∂θ/∂(μ_H, μ_U) =
   (−1/μ_U, μ_H/μ_U²) and Ω the long-run covariance of (H_t, U_t), the delta-
   method variance of θ̂ is g'Ωg/T. Because μ̂_H = (1 − θ̂)μ̂_U,
   g'(H_t − μ̂_H, U_t − μ̂_U)' = −(H_t − (1 − θ̂)U_t)/μ̂_U, so

       se_δ(θ̂) = ŝe(θ̂)/|μ̂_U|,

   where ŝe(ϑ) is the Newey–West standard error of the mean of H − (1 − ϑ)U,
   the statistic inverted for the test-inversion set (Fieller, 1954;
   Newey and West, 1987, 1994, through qef.stats.hac). Since
   μ̂_H − (1 − ϑ)μ̂_U = μ̂_U(ϑ − θ̂), the delta-method interval
   θ̂ ± 1.96 se_δ(θ̂) is {ϑ : |μ̂_H − (1 − ϑ)μ̂_U| ≤ 1.96 ŝe(θ̂)}, and the
   test-inversion set is the same inequality with ŝe(ϑ) in place of ŝe(θ̂).
   For a fixed lag, ŝe(ϑ)² is a quadratic in ϑ whose leading coefficient is
   ŝe_U², the squared standard error of μ̂_U; the automatic lag of
   H − (1 − ϑ)U tends to that of U as |ϑ| grows, because the bandwidth rule is
   scale invariant. For large |ϑ| the test-inversion inequality therefore
   tends to |μ̂_U| ≤ 1.96 ŝe_U: when μ̂_U is not significant every ϑ far
   enough from θ̂ is retained and the set is unbounded. The delta-method
   interval freezes the standard error at θ̂, so its half-width
   1.96 ŝe(θ̂)/|μ̂_U| is finite whenever μ̂_U ≠ 0, however small μ̂_U is
   relative to its standard error. A procedure whose sets are bounded with
   probability one has zero worst-case coverage for a ratio of this kind
   (Gleser and Hwang, 1987; Dufour, 1997; R9(c) of the paper), so the
   delta-method interval is misleading in particular when the test-inversion
   set is unbounded. Away from θ̂ the two procedures can also disagree about a
   single value such as ϑ = 0, because ŝe(ϑ) differs from ŝe(θ̂) whether or
   not the set is bounded.

4. Closed-form checks of the moment code (moments_1m.csv, written by
   estimate_moments.py through qef.fx.moments): closed_form_rel_diff, the
   largest relative difference between the adaptive tail quadrature and the
   closed-form R4 bounds for the weights K⁻² and the logarithmic weights at
   each currency-month's own boundary prices and exponents; tail_quad_err, the
   largest error estimate of the adaptive quadrature relative to its integral,
   which is quad's own estimate and so confirms only that no tail integral
   stopped at quad's subdivision limit before meeting its relative tolerance
   of 1e-12 (closed_form_rel_diff is the independent check of accuracy); and
   simpson_converged, whether halving the Simpson step moved every middle
   contract by less than 1e-10 relative.

Rules fixed before computing the results, except where R10 and R11 say
otherwise (the design leaves these open; the amendments after R12 were made
after an independent review):

R1  SOFR window: the primary-sample month-ends from the first at which
    USDSROIS1M= is available under the design's month-end rule (two-sided on
    the month-end, or the last two-sided quote within the five preceding New
    York business days) to August 2026. A month-end in the window at which
    SOFR is missing keeps Fed Funds OIS in the spliced series and is excluded
    from the overlap comparison; such month-ends are counted.
R2  SOFR discounting: mid = (bid + ask)/2, and the same money-market formula,
    ACT/360 basis and day count d as for Fed Funds OIS. Only the four
    USD-quoted pairs are recalibrated; the stored smiles of the dollar-base
    pairs are reused unchanged.
R3  SOFR recalibration runs _calibrate_currency per currency over every
    month-end of smile_inputs.csv (July 2010 to August 2026, the range of the
    stored panel), with inputs identical to the stored ones before the window.
    Each month-end is then warm-started from the previous month's solution
    exactly as in the published calibration, so the SOFR smiles differ from
    the stored ones only through the rate. The same call on the unchanged Fed
    Funds inputs is compared with smile_panel.csv as a check of the call path
    (largest parameter gap and whether the statuses agree), and the SOFR
    smiles before the window are compared with the stored ones. The Fed Funds
    E1 series is the published one, recomputed from the stored inputs and
    panel and checked against e1_series.csv.
R4  SOFR inference: E1 as specified (primary sample, market-strangle reading,
    10Δ), regimes as in E1, regime means with Newey–West standard errors and
    the difference with the Newey–West standard error of the regime-dummy
    coefficient and the stationary-bootstrap 95% percentile interval (9,999
    draws, seed and block lengths as in estimate_e1.regime_summary). Reported
    (a) over the overlap under both discountings and (b) over the whole
    primary sample with SOFR spliced in from the start of the window.
R5  A conclusion changes if, for the regime difference of φ (and, as a
    supplementary check, of C_skew and FD), the sign, or whether the
    Newey–West interval (difference ± 1.96 s.e.) excludes zero, or whether the
    bootstrap interval excludes zero, differs between the two discountings.
    The largest absolute changes in monthly φ and C_skew are taken over the
    month-ends at which both series have status ok.
R6  TIFO month-end comparison: mids = (bid + ask)/2 of two-sided quotes,
    sampled at New York month-ends from the primary start to August 2026 under
    the same five-business-day rule; TIFO minus the reference (composite or
    Fenics), in volatility points, at month-ends where both are available.
    Equal mids: absolute difference below 1e-12, as in the audit.
R7  TIFO smile: EUR month-ends at which all three TIFO calibration quotes are
    available under the rule and the composite EUR row has the calibration
    set. Spot, forward, rate and dates come from the composite inputs. TIFO
    quotes no 10Δ, so no held-out prediction is made. The market-strangle
    smile is used, as in E1.
R8  TIFO leg skew cost: 10Δ, both sides (a long EUR leg is protected by a put
    on EURUSD, a short leg by a call), in bp of notional per month, compared
    with the composite at the same month-ends and with Fenics where a Fenics
    market-reading smile exists (status ok or not_converged). Effect on E1:
    at the E1 month-ends (primary, market, 10Δ, status ok) where EUR is a leg
    and a TIFO smile exists, the EUR leg's contribution to C_skew is replaced
    by the TIFO one, with FD and the other legs unchanged. No regime inference
    is made on this subsample.
R9  Delta method: the lag is chosen by the automatic Newey–West (1994)
    procedure applied to the linearisation z_t = −(H_t − (1 − θ̂)U_t)/μ̂_U,
    which makes se_δ equal to the Newey–West standard error of the mean of z
    and to ŝe(θ̂)/|μ̂_U|; the explicit 2 × 2 Bartlett long-run covariance of
    (H, U) at that lag is computed as a check. Interval θ̂ ± 1.96 se_δ.
R10 Test-inversion sets for the three hedges use theta_confidence_set of
    estimate_stage4.py (grid −10 to 10 in steps of 0.001), the rule behind the
    published 10Δ set; for the 25Δ and ATM hedges this rule is applied here
    for the first time. An unbounded set is described by the grid points it
    rejects: none (the whole grid) or an interval [a, b] of rejected points,
    in which case the set is the grid minus [a, b] (first written as the
    complement of (a, b); amendment A1). Whether the delta-method interval
    contains θ_0 is reported for the 10Δ hedge only (θ_0 is the mean of the
    theta0 column of stage4_months.csv, the reference of the published E3;
    the file holds none for the other hedges). Added after a first run
    with 99 bootstrap draws had shown the sets, and therefore not fixed in
    advance: whether each test-inversion set contains 0 and θ_0, whether the
    delta-method interval contains 0 (missing from this list as first
    written; amendment A1), and the t statistics of the test of ϑ = 0 under
    both procedures.
R11 Closed-form checks: over the currency-months of moments_1m.csv with status
    ok, the maximum of each diagnostic where it is finite, per tail setting,
    with its count; Simpson convergence as a share of those currency-months.
    The number of tail error estimates above the quadrature tolerance of 1e-12
    was added after the first run.
R12 Only pooled statistics (means, standard deviations, counts, extremes,
    test statistics, intervals) and overlap windows go to the summary. The
    month-level series and calibrated parameters go to CSV files next to it.

Amendments of 27 September 2026, made after an independent review of the full
run. They correct the record and add post hoc items; I tuned no tolerance,
grid, lag rule or estimate to the results, and every number reported before
them is unchanged.

A1  Corrections. In R10, a and b are the smallest and largest rejected grid
    points, so both are rejected and the set is the grid minus the closed
    interval [a, b]; "the complement of (a, b)", as first written, put them
    in the set. R10's list of items added after the first run left out
    whether the delta-method interval contains 0, which no rule covers; it
    is now listed. The opening paragraph and the heading of the rules said
    that every rule was fixed before any result was computed; both now
    except the items of R10 and R11, and the summary marks those items post
    hoc. The delta-method t statistic at ϑ = 0 is computed as
    −sign(μ̂_U) θ̂/se_δ, the mean of H − U over |μ̂_U| se_δ, so that it is
    oriented as the test-inversion statistic, the mean of H − U over its own
    Newey–West standard error; with μ̂_U > 0 its value is unchanged. The
    summary's account of why the delta-method interval misleads, and of what
    tail_quad_err checks, is also corrected; neither changes a number.
A2  Post hoc, added after the review. For each hedge (10Δ, 25Δ and ATM), a
    studentised (bootstrap-t) stationary bootstrap of the mean of H − U:
    9,999 draws (the --bootstrap default), seed 20260924, indices from
    stationary_indices with the expected block length from
    optimal_block_length (qef.stats.bootstrap, through bootstrap_distribution
    as in estimate_e1.regime_summary), and in each draw
    t* = (mean* − mean)/se*, with se* the Newey–West standard error at the
    automatic lag (mean_and_se) recomputed on the draw. Reported: the
    equal-tailed 95% interval [mean − q(0.975) se, mean − q(0.025) se], with
    q the quantiles of t* and se the Newey–West standard error of the sample,
    and the share of t* at or below the observed t, a one-sided bootstrap p
    value. The percentile interval of the mean is not reported, because it
    does not correct for the skewness of H − U. Reported as diagnostics: the
    skewness of H − U; the t statistic after dropping the month with the
    largest H − U; the t statistic at ϑ = 0, and for the 10Δ hedge at θ_0, at
    each fixed Bartlett lag from 0 to 12; the t statistic of the mean of
    H10 − U over the extended sample of stage4_months.csv; and, as a check of
    the sign reading, the largest gap between H10 − U and the E3 total
    c1 + c2 + c3. The conclusion I draw from them in the summary was written
    after seeing them, and delta_check prints it only while the numbers still
    support each of its clauses.

Outputs are LSEG-derived and are written to data/private/results/<retrieval
date>/: design_checks_summary.md and design_checks_sofr_summary.csv (pooled
statistics only), and the private month-level files design_checks_sofr_e1.csv,
design_checks_sofr_smiles.csv, design_checks_tifo_smiles.csv and
design_checks_tifo_legs.csv. The calibration progress lines on standard error
name currencies, months and statuses only. The whole run takes about two
minutes on eight cores.
"""

from __future__ import annotations

import os

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse  # noqa: E402
import sys  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from calibrate_smiles import _calibrate_currency  # noqa: E402
from estimate_e1 import PRIMARY_START, REGIME_BREAK, _vol, e1_series, regime_summary  # noqa: E402
from estimate_stage4 import theta_confidence_set  # noqa: E402

from qef.data.panel import ny_month_ends, read_raw, sample_month_ends  # noqa: E402
from qef.data.smile_inputs import MM_BASIS  # noqa: E402
from qef.fx.conventions import G10  # noqa: E402
from qef.fx.crash import leg_skew_cost  # noqa: E402
from qef.stats.bootstrap import bootstrap_distribution, optimal_block_length  # noqa: E402
from qef.stats.hac import long_run_variance, mean_and_se  # noqa: E402

SAMPLE_END = pd.Timestamp("2026-08-31")
DELTA = 0.10
N_LEGS = 3
USD_QUOTED = tuple(c for c, conv in G10.items() if not conv.usd_base)
SOFR = ("rates", "USDSROIS1M=")
TIFO_QUOTES = ("atm", "rr25", "bf25")
VOL_FILES = {  # quote -> (block, RIC stem without contributor suffix)
    "atm": ("vol_atm", "EUR1MO="), "rr25": ("vol_rr25", "EUR1MRR="), "bf25": ("vol_bf25", "EUR1MBF="),
}
CONTRIBUTORS = {"tifo": "TIFO", "composite": "", "fenics": "FN"}
EQUAL_TOL = 1e-12
BP = 1e4
BOOT_SEED = 20260924  # the seed of estimate_e1.regime_summary (A2)
FIXED_LAGS = range(13)  # Bartlett lags 0 to 12 (A2)


# ---------------------------------------------------------------------------
# Shared helpers


def discount_factors(r_quote, days, S, F, basis=MM_BASIS["USD"]):
    """D_q = 1/(1 + r d/basis) and D_b = F D_q/S, as in qef.data.smile_inputs."""
    df_q = 1.0 / (1.0 + r_quote * days / basis)
    return df_q, F * df_q / S


def calibrate_groups(groups: dict, workers: int) -> dict:
    """Run calibrate_smiles._calibrate_currency on each labelled group of input rows."""
    labels = list(groups)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        parts = list(ex.map(_calibrate_currency, [groups[k] for k in labels]))
    return {k: pd.DataFrame(p) for k, p in zip(labels, parts)}


def ci_flags(row) -> dict:
    """Sign and whether each 95% interval of a regime difference excludes zero."""
    d, se = row["difference"], row["se_difference_hac"]
    return {"sign": int(np.sign(d)), "hac_excludes_0": bool(abs(d) > 1.96 * se),
            "boot_excludes_0": bool(row["boot_lo"] > 0 or row["boot_hi"] < 0)}


def fmt_summary(r, scale=1.0) -> str:
    return (f"zero-rate {scale * r['mean_zero_rate']:.4f} ({scale * r['se_zero_rate']:.4f}, n={r['n_zero_rate']}); "
            f"hiking {scale * r['mean_hiking']:.4f} ({scale * r['se_hiking']:.4f}, n={r['n_hiking']}); "
            f"difference {scale * r['difference']:.4f} (HAC s.e. {scale * r['se_difference_hac']:.4f}, lag {r['hac_lag']}; "
            f"bootstrap 95% [{scale * r['boot_lo']:.4f}, {scale * r['boot_hi']:.4f}])")


# ---------------------------------------------------------------------------
# 1. SOFR OIS


def sofr_inputs(inputs: pd.DataFrame, raw: Path, month_ends):
    """Inputs with SOFR OIS in place of Fed Funds OIS for USD-quoted pairs from the SOFR start (R1, R2)."""
    sofr = sample_month_ends(read_raw(raw / SOFR[0] / f"{SOFR[1]}.csv"), month_ends)
    avail = sofr["status"] != "missing"
    start = sofr.index[avail].min()
    window = sofr.index[sofr.index >= start]
    use = window[avail.reindex(window).to_numpy()]
    out = inputs.copy()
    m = out.currency.isin(USD_QUOTED) & out.date.isin(use)
    out.loc[m, "r_quote"] = out.loc[m, "date"].map(sofr["mid"]).to_numpy() / 100.0
    df_q, df_b = discount_factors(out.loc[m, "r_quote"], out.loc[m, "days"], out.loc[m, "S"], out.loc[m, "F"])
    out.loc[m, "df_quote"], out.loc[m, "df_base"] = df_q, df_b
    info = {"start": start, "end": window.max(), "n_window": len(window), "n_missing": int(len(window) - len(use)),
            "status_counts": sofr.loc[window, "status"].value_counts().to_dict()}
    return out, use, m, info


def sofr_check(inputs, panel, e1_pub, raw, month_ends, B, workers, out_dir):
    L = ["## 1. SOFR OIS in place of Fed Funds OIS (E1: primary sample, market reading, 10-delta)", ""]
    # Reproduction of the discount-factor formula on the stored inputs.
    usdq = inputs.currency.isin(USD_QUOTED)
    q_chk, b_chk = discount_factors(inputs.loc[usdq, "r_quote"], inputs.loc[usdq, "days"], inputs.loc[usdq, "S"], inputs.loc[usdq, "F"])
    formula_gap = max(float(np.nanmax(np.abs(q_chk - inputs.loc[usdq, "df_quote"]))),
                      float(np.nanmax(np.abs(b_chk - inputs.loc[usdq, "df_base"]))))

    s_in, use, changed, info = sofr_inputs(inputs, raw, month_ends)
    ff_rate = inputs.loc[changed, "r_quote"].to_numpy()
    so_rate = s_in.loc[changed, "r_quote"].to_numpy()
    by_date = pd.DataFrame({"date": inputs.loc[changed, "date"].to_numpy(), "diff": so_rate - ff_rate}).groupby("date")["diff"].first()
    rel_db = pd.Series(np.abs(s_in.loc[changed, "df_base"].to_numpy() / inputs.loc[changed, "df_base"].to_numpy() - 1.0)).dropna()
    # complete_calib is kept: it changes under SOFR only where the Fed Funds rate is missing in the window.
    n_ff_missing = int(np.isnan(ff_rate).sum())

    # Recalibrate the four USD-quoted smiles over the whole stored history (R3), with Fed Funds
    # (call-path check) and with SOFR spliced in from the start of the window.
    groups = {}
    for c in USD_QUOTED:
        groups[("ff", c)] = inputs[inputs.currency == c]
        groups[("sofr", c)] = s_in[s_in.currency == c]
    fits = calibrate_groups(groups, workers)
    ff_fit = pd.concat([fits[("ff", c)] for c in USD_QUOTED], ignore_index=True)
    so_fit = pd.concat([fits[("sofr", c)] for c in USD_QUOTED], ignore_index=True)
    ff_fit["rate"], so_fit["rate"] = "fed_funds_ois", "sofr_ois"
    pd.concat([ff_fit, so_fit], ignore_index=True).to_csv(out_dir / "design_checks_sofr_smiles.csv", index=False)

    key = ["currency", "date", "reading"]
    par = ["alpha", "rho", "nu"]
    stored = panel.set_index(key)

    def gap_to_stored(fit):
        j = fit.set_index(key).join(stored[par + ["status"]], rsuffix="_stored", how="inner")
        g = np.abs(j[par].to_numpy() - j[[f"{x}_stored" for x in par]].to_numpy())
        return (float(np.nanmax(g)) if np.isfinite(g).any() else 0.0), bool((j["status"] == j["status_stored"]).all()), len(j)

    replica_gap, replica_status_same, n_replica = gap_to_stored(ff_fit)
    pre_gap, pre_status_same, n_pre = gap_to_stored(so_fit[so_fit.date < info["start"]])
    in_win = so_fit.date.isin(use)
    so_ok = so_fit[(so_fit.reading == "market") & in_win]
    so_status = so_ok["status"].value_counts().to_dict()
    so_done = so_ok[so_ok.status.isin(["ok", "not_converged"])]
    so_resid = float(so_done["max_abs_residual"].max())
    so_arb = float(so_done[["arb_ok_inner", "arb_ok_wide"]].all(axis=1).mean())
    win_gap = so_fit[in_win & (so_fit.reading == "market")].set_index(key).join(
        stored[par], rsuffix="_stored", how="inner")
    par_change = {x: float(np.nanmax(np.abs(win_gap[x] - win_gap[f"{x}_stored"]))) for x in par}

    # SOFR panel: stored panel with the recalibrated USD-quoted smiles substituted.
    drop = panel.currency.isin(USD_QUOTED)
    panel_s = pd.concat([panel[~drop], so_fit.drop(columns="rate")], ignore_index=True)

    prim_ff = inputs[inputs.date >= PRIMARY_START]
    prim_so = s_in[s_in.date >= PRIMARY_START]
    e_ff = e1_series(prim_ff, panel, "market", DELTA, N_LEGS)
    e_so = e1_series(prim_so, panel_s, "market", DELTA, N_LEGS)
    pub = e1_pub.set_index("date")
    chk = e_ff[e_ff.status == "ok"].set_index("date")
    pub_gap = float(np.max(np.abs(chk["phi"] - pub.loc[chk.index, "phi"]))) if len(chk) else np.nan
    pub_same_months = set(chk.index) == set(pub[pub.status == "ok"].index)

    both = e_ff[["date", "status", "longs", "shorts", "phi", "C_skew", "FD"]].merge(
        e_so[["date", "status", "phi", "C_skew", "FD"]], on="date", suffixes=("_ff", "_sofr"))
    both["in_overlap"] = both.date.isin(use)
    both.to_csv(out_dir / "design_checks_sofr_e1.csv", index=False)
    ov = both[both.in_overlap & (both.status_ff == "ok") & (both.status_sofr == "ok")]
    dphi, dc, dfd = ov.phi_sofr - ov.phi_ff, BP * (ov.C_skew_sofr - ov.C_skew_ff), BP * (ov.FD_sofr - ov.FD_ff)
    n_legs_aff = int(sum(sum(c in USD_QUOTED for c in (a + " " + b).split()) for a, b in zip(ov.longs, ov.shorts)))
    n_status_diff = int((both[both.in_overlap].status_ff != both[both.in_overlap].status_sofr).sum())

    # Inference over the overlap (both discountings) and over the whole primary sample (SOFR spliced).
    rows = []
    for label, s in (("fed_funds", e_ff), ("sofr", e_so)):
        for col in ("phi", "C_skew", "FD"):
            ovs = s[s.date.isin(ov.date)]
            rows.append({"window": "overlap", "rate": label, **regime_summary(ovs, col, B)})
            rows.append({"window": "primary", "rate": label, **regime_summary(s, col, B)})
    summ = pd.DataFrame(rows)
    summ.to_csv(out_dir / "design_checks_sofr_summary.csv", index=False)

    L += [f"- SOFR window under the month-end rule: {info['start']:%Y-%m-%d} to {info['end']:%Y-%m-%d}, "
          f"{info['n_window']} month-ends ({', '.join(f'{k} {v}' for k, v in info['status_counts'].items())}); "
          f"month-ends missing and excluded: {info['n_missing']}.",
          f"- Overlap used: {len(ov)} month-ends with status ok under both rates "
          f"({int((ov.date < REGIME_BREAK).sum())} zero-rate, {int((ov.date >= REGIME_BREAK).sum())} hiking); "
          f"month-ends whose status differs between the rates: {n_status_diff}. "
          f"Legs in USD-quoted pairs over the overlap: {n_legs_aff} of {2 * N_LEGS * len(ov)}.",
          f"- SOFR minus Fed Funds OIS, one-month, over the window (bp a year): mean {BP * by_date.mean():.2f}, "
          f"s.d. {BP * by_date.std():.2f}, min {BP * by_date.min():.2f}, max {BP * by_date.max():.2f}, "
          f"largest absolute {BP * by_date.abs().max():.2f} ({len(by_date)} month-ends).",
          f"- Largest relative change in D_b: {rel_db.max():.2e} (mean {rel_db.mean():.2e}, {len(rel_db)} currency-months); "
          f"currency-months in the window with Fed Funds OIS missing: {n_ff_missing}.",
          f"- Checks: discount-factor formula against the stored inputs, largest gap {formula_gap:.1e}; "
          f"Fed Funds recalibration through the same call path against the stored panel ({n_replica} calibrations, "
          f"both readings), largest parameter gap {replica_gap:.1e}, statuses identical: {replica_status_same}; "
          f"SOFR run before the window against the stored panel ({n_pre} calibrations), largest parameter gap "
          f"{pre_gap:.1e}, statuses identical: {pre_status_same}; Fed Funds E1 recomputed against "
          f"e1_series.csv, largest phi gap {pub_gap:.1e}, same month-ends: {pub_same_months}.",
          f"- SOFR recalibrations in the window (market reading, {len(so_ok)} currency-months): statuses {so_status}; "
          f"largest absolute residual {so_resid:.1e}; share passing both arbitrage checks {so_arb:.1%}; "
          f"largest absolute change in the SABR parameters against the stored smiles: "
          + ", ".join(f"{x} {v:.2e}" for x, v in par_change.items()) + ".",
          "",
          "Change in monthly values, SOFR minus Fed Funds, over the overlap:", "",
          f"- phi: mean {dphi.mean():.2e}, largest absolute {dphi.abs().max():.2e}, "
          f"largest absolute relative {(dphi / ov.phi_ff).abs().max():.2e}",
          f"- C_skew (bp): mean {dc.mean():.2e}, largest absolute {dc.abs().max():.2e}, "
          f"largest absolute relative {(dc / (BP * ov.C_skew_ff)).abs().max():.2e}",
          f"- FD (bp): largest absolute {dfd.abs().max():.1e} (FD does not depend on the USD rate)",
          "", "Regime means (Newey-West s.e.) and differences, hiking minus zero-rate. phi dimensionless; C_skew and FD in bp:", ""]
    flags = {}
    for window, title in (("overlap", "Overlap"), ("primary", "Whole primary sample, SOFR spliced in from the start of the window")):
        L += [f"{title}:", ""]
        for col in ("phi", "C_skew", "FD"):
            scale = 1.0 if col == "phi" else BP
            for rate in ("fed_funds", "sofr"):
                r = summ[(summ.window == window) & (summ.rate == rate) & (summ.variable == col)].iloc[0]
                flags[(window, col, rate)] = ci_flags(r)
                L.append(f"- {col}, {rate}: {fmt_summary(r, scale)}")
        L.append("")
    changes = [(w, c) for (w, c, rate) in flags if rate == "sofr" and flags[(w, c, "sofr")] != flags[(w, c, "fed_funds")]]
    L += ["Conclusion rule (R5): sign, Newey-West interval excluding zero, bootstrap interval excluding zero.", ""]
    for (w, c, rate), f in flags.items():
        if rate == "sofr":
            L.append(f"- {w}, {c}: Fed Funds {flags[(w, c, 'fed_funds')]}; SOFR {f}")
    L += ["", f"- Conclusions that change: {'none' if not changes else changes}", ""]
    return L


# ---------------------------------------------------------------------------
# 2. Named broker (TIFO)


def month_end_quotes(raw: Path, month_ends) -> dict:
    """Month-end samples of the EUR one-month calibration quotes per contributor (R6)."""
    out = {}
    for who, suffix in CONTRIBUTORS.items():
        for q, (block, stem) in VOL_FILES.items():
            out[(who, q)] = sample_month_ends(read_raw(raw / block / f"{stem}{suffix}.csv"), month_ends)
    return out


def tifo_check(inputs, panel, fn_in, fn_panel, e1_pub, raw, month_ends, audit_dir, workers, out_dir):
    L = ["## 2. Named-broker check: TIFO against the composite and Fenics (EUR, one month)", ""]
    spl = pd.read_csv(audit_dir / "splice.csv")
    ti = spl[spl.comparison.str.startswith("tifo") & (spl.n > 0)]
    L += ["Daily mids on the overlap, from the private audit (splice.csv; volatility points, TIFO minus reference):", "",
          "| comparison | quote | n days | overlap | mean | s.d. | mean abs | share equal | mean spread TIFO / ref. | n month-ends | month-end mean | month-end s.d. |",
          "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in ti.itertuples():
        L.append(f"| {r.comparison} | {r.quote} | {r.n} | {r.first_overlap} to {r.last_overlap} | {r.mean:+.4f} | {r.sd:.4f} | "
                 f"{r.mean_abs:.4f} | {r.share_equal_mid:.3f} | {r.mean_spread:.3f} / {r.mean_spread_reference:.3f} | "
                 f"{r.n_month_ends} | {r.mean_month_ends:+.4f} | {r.sd_month_ends:.4f} |")
    L += ["", "The audit's month-end columns use quotes dated on the month-end only. Under the month-end rule of the design "
          "(five-business-day substitution), primary-sample month-ends:", "",
          "| comparison | quote | n month-ends | overlap | TIFO substituted | mean | s.d. | mean abs | share equal |",
          "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    me = month_end_quotes(raw, month_ends)
    for ref in ("composite", "fenics"):
        for q in TIFO_QUOTES:
            a, b = me[("tifo", q)], me[(ref, q)]
            ok = (a.status != "missing") & (b.status != "missing")
            d = (a["mid"] - b["mid"])[ok]
            L.append(f"| tifo-{ref} | {q} | {len(d)} | {d.index.min():%Y-%m} to {d.index.max():%Y-%m} | "
                     f"{int((a.status[ok] == 'substituted').sum())} | {d.mean():+.4f} | {d.std():.4f} | {d.abs().mean():.4f} | "
                     f"{(d.abs() < EQUAL_TOL).mean():.3f} |")

    # TIFO smiles (R7).
    tifo_ok = np.logical_and.reduce([me[("tifo", q)].status != "missing" for q in TIFO_QUOTES])
    tifo_dates = me[("tifo", "atm")].index[tifo_ok]
    eur = inputs[(inputs.currency == "EUR") & inputs.date.isin(tifo_dates) & inputs.complete_calib].copy()
    for q in TIFO_QUOTES:
        eur[q] = eur.date.map(me[("tifo", q)]["mid"]).to_numpy() / 100.0
    eur["has_10d"] = False
    eur[["rr10", "bf10"]] = np.nan
    t_fit = calibrate_groups({"tifo": eur}, workers)["tifo"]
    t_fit.to_csv(out_dir / "design_checks_tifo_smiles.csv", index=False)
    t_mkt = t_fit[(t_fit.reading == "market") & t_fit.status.isin(["ok", "not_converged"])].set_index("date")

    # Leg skew costs (R8).
    ci = inputs[inputs.currency == "EUR"].set_index("date")
    fc = panel[(panel.currency == "EUR") & (panel.reading == "market") & panel.status.isin(["ok", "not_converged"])].set_index("date")
    fi = fn_in[(fn_in.currency == "EUR") & fn_in.complete_calib].set_index("date")
    ff = fn_panel[(fn_panel.currency == "EUR") & (fn_panel.reading == "market") & fn_panel.status.isin(["ok", "not_converged"])].set_index("date")
    ti_in = eur.set_index("date")

    def cost(r, f, long_leg):
        try:
            _, v_s, v_f = leg_skew_cost(_vol(r, f), r.F, r.tau, r.df_base, "EUR", r.atm, long_leg, DELTA)
            return BP * (v_s - v_f)
        except ValueError:
            return np.nan

    rows = []
    for t in t_mkt.index:
        if t not in fc.index:
            continue
        for long_leg in (True, False):
            rows.append({"date": t, "side": "long" if long_leg else "short",
                         "tifo": cost(ti_in.loc[t], t_mkt.loc[t], long_leg),
                         "composite": cost(ci.loc[t], fc.loc[t], long_leg),
                         "fenics": cost(fi.loc[t], ff.loc[t], long_leg) if (t in fi.index and t in ff.index) else np.nan})
    legs = pd.DataFrame(rows)
    legs.to_csv(out_dir / "design_checks_tifo_legs.csv", index=False)
    t_status = t_fit[t_fit.reading == "market"]["status"].value_counts().to_dict()
    t_ok = t_fit[(t_fit.reading == "market") & (t_fit.status == "ok")]
    L += ["", f"TIFO smiles (market reading): {len(eur)} EUR month-ends with the TIFO calibration set, "
          f"{eur.date.min():%Y-%m} to {eur.date.max():%Y-%m} ({int((eur.date < REGIME_BREAK).sum())} zero-rate, "
          f"{int((eur.date >= REGIME_BREAK).sum())} hiking); statuses {t_status}; largest absolute residual "
          f"{t_ok.max_abs_residual.max():.1e}; share passing both arbitrage checks "
          f"{t_ok[['arb_ok_inner', 'arb_ok_wide']].all(axis=1).mean():.1%}.", "",
          "10-delta leg skew cost V_smile - V_flat of the EUR protective option (bp of notional per month):", "",
          "| side | comparison | n | mean of reference | mean difference | s.d. | mean abs | largest abs | correlation |",
          "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for side in ("long", "short", "both"):
        g = legs if side == "both" else legs[legs.side == side]
        for a, b in (("tifo", "composite"), ("tifo", "fenics"), ("fenics", "composite")):
            x = g[[a, b]].dropna()
            d = x[a] - x[b]
            L.append(f"| {side} | {a} - {b} | {len(x)} | {x[b].mean():.3f} | {d.mean():+.3f} | {d.std():.3f} | "
                     f"{d.abs().mean():.3f} | {d.abs().max():.3f} | {np.corrcoef(x[a], x[b])[0, 1]:.4f} |")

    # Effect on E1 (R8).
    e1 = e1_pub[e1_pub.status == "ok"]
    eff = []
    lk = legs.set_index(["date", "side"])
    for r in e1.itertuples():
        side = "long" if "EUR" in r.longs.split() else ("short" if "EUR" in r.shorts.split() else None)
        if side is None or (r.date, side) not in lk.index:
            continue
        x = lk.loc[(r.date, side)]
        if not np.isfinite(x.tifo):
            continue
        dC = (x.tifo - x.composite) / N_LEGS / BP
        eff.append({"date": r.date, "side": side, "C_skew": r.C_skew, "FD": r.FD, "phi": r.phi,
                    "dC_bp": BP * dC, "phi_tifo": (r.C_skew + dC) / r.FD})
    eff = pd.DataFrame(eff)
    L += [""]
    if len(eff):
        dphi = eff.phi_tifo - eff.phi
        L += [f"Effect on E1 (primary, market reading, 10-delta) at the {len(eff)} E1 month-ends where EUR is a leg and a "
              f"TIFO smile exists ({int((eff.side == 'long').sum())} long, {int((eff.side == 'short').sum())} short; "
              f"{int((eff.date < REGIME_BREAK).sum())} zero-rate, {int((eff.date >= REGIME_BREAK).sum())} hiking):", "",
              f"- change in C_skew (bp): mean {eff.dC_bp.mean():+.4f}, s.d. {eff.dC_bp.std():.4f}, largest absolute {eff.dC_bp.abs().max():.4f}; "
              f"mean C_skew {BP * eff.C_skew.mean():.3f} (composite) against {BP * eff.C_skew.mean() + eff.dC_bp.mean():.3f} (TIFO for EUR)",
              f"- change in phi: mean {dphi.mean():+.5f}, s.d. {dphi.std():.5f}, largest absolute {dphi.abs().max():.5f}; "
              f"mean phi {eff.phi.mean():.4f} against {eff.phi_tifo.mean():.4f}", ""]
    else:
        L += ["No E1 month-end has EUR as a leg with a TIFO smile.", ""]
    return L


# ---------------------------------------------------------------------------
# 3. Delta-method intervals for theta_UB


def bartlett_lrcov(X: np.ndarray, lag: int) -> np.ndarray:
    """Bartlett-weighted long-run covariance of the columns of X (demeaned here)."""
    X = X - X.mean(axis=0)
    T = X.shape[0]
    S = X.T @ X / T
    for j in range(1, lag + 1):
        G = X[j:].T @ X[: T - j] / T
        S += (1.0 - j / (lag + 1.0)) * (G + G.T)
    return S


def delta_method(H, U) -> dict:
    """θ̂ = 1 − μ̂_H/μ̂_U with the delta-method standard error of R9."""
    H, U = np.asarray(H, float), np.asarray(U, float)
    T = len(H)
    mH, mU = H.mean(), U.mean()
    theta = 1.0 - mH / mU
    z = -(H - (1.0 - theta) * U) / mU
    lrv, lag = long_run_variance(z - z.mean())
    se = float(np.sqrt(lrv / T))
    g = np.array([-1.0 / mU, mH / mU**2])
    se_explicit = float(np.sqrt(g @ bartlett_lrcov(np.column_stack([H, U]), lag) @ g / T))
    se_inv = mean_and_se(H - (1.0 - theta) * U)
    u = mean_and_se(U)
    return {"theta": theta, "se_delta": se, "lag": lag, "se_explicit_2x2": se_explicit,
            "se_via_inversion_stat": se_inv["se"] / abs(mU), "lag_inversion_stat": se_inv["lag"],
            "lo": theta - 1.96 * se, "hi": theta + 1.96 * se, "mean_U": mU, "se_U": u["se"], "t_U": mU / u["se"]}


def rejected_range(H, U, grid=np.linspace(-10, 10, 20001)):
    """Smallest and largest grid points rejected by the test-inversion rule of theta_confidence_set."""
    rej = []
    for th in grid:
        r = mean_and_se(H - (1 - th) * U)
        rej.append(abs(r["mean"]) > 1.96 * r["se"])
    rej = np.array(rej)
    if not rej.any():
        return np.nan, np.nan, 0
    return float(grid[rej].min()), float(grid[rej].max()), int(np.all(np.diff(np.where(rej)[0]) == 1))


def inversion_accepts(H, U, v) -> tuple[bool, float]:
    """Whether the test-inversion rule of theta_confidence_set retains v, with the t statistic at v."""
    r = mean_and_se(H - (1 - v) * U)
    return bool(abs(r["mean"]) <= 1.96 * r["se"]), r["mean"] / r["se"]


def skewness(x) -> float:
    """Third central moment over the cubed standard deviation, both with divisor T."""
    e = np.asarray(x, float) - np.mean(x)
    return float(np.mean(e**3) / np.mean(e**2) ** 1.5)


def t_fixed_lag(x, lag: int) -> float:
    """t statistic of the mean of x with the Newey–West standard error at a fixed Bartlett lag."""
    x = np.asarray(x, float)
    lrv, _ = long_run_variance(x - x.mean(), lag)
    return float(x.mean() / np.sqrt(lrv / x.shape[0]))


def studentised_bootstrap(x, B: int, seed: int = BOOT_SEED) -> dict:
    """Studentised stationary bootstrap of the mean of x (A2, post hoc).

    Each draw gives t* = (mean* − mean)/se*, with se* the Newey–West standard
    error at the automatic lag of the draw. The equal-tailed 95% interval is
    [mean − q(0.975) se, mean − q(0.025) se], and share_below is the share of t*
    at or below the observed t = mean/se.
    """
    x = np.asarray(x, float)
    r = mean_and_se(x)
    m, se = r["mean"], r["se"]
    t_star = bootstrap_distribution(lambda xs: (xs.mean() - m) / mean_and_se(xs)["se"], [x], B=B, seed=seed)
    q_lo, q_hi = np.percentile(t_star, [2.5, 97.5])
    t = m / se
    return {"t": t, "block": optimal_block_length(x), "q_lo": float(q_lo), "q_hi": float(q_hi),
            "lo": m - q_hi * se, "hi": m - q_lo * se, "share_below": float(np.mean(t_star <= t))}


def delta_check(s4, B):
    L = ["## 3. Delta-method intervals for theta_UB (primary sample)", ""]
    pm = s4[(s4["sample"] == "primary") & (s4.status == "ok")].sort_values("date")
    theta0 = pm.theta0.mean()
    U = pm.U.to_numpy()
    u = mean_and_se(U)
    L += [f"Months: {len(pm)}. Mean HML^U {BP * u['mean']:.2f} bp, Newey-West s.e. {BP * u['se']:.2f} bp (lag {u['lag']}), "
          f"t = {u['mean'] / u['se']:.3f}; the test-inversion set is bounded only if |t| exceeds 1.96 "
          f"(asymptotically in |theta|). Reference theta_0 for the 10-delta hedge (mean over months of the average absolute forward "
          f"delta of the 10-delta options, from stage4_months.csv): {theta0:.4f}; stage4_months.csv holds no theta_0 for the 25-delta "
          f"and ATM hedges, so those cells are left empty.", "",
          "Columns and notes marked post hoc were added after a first run with 99 bootstrap draws (R10) or after an "
          "independent review (A2 of the script's docstring).", "",
          "| hedge | theta_UB | delta-method s.e. | lag | 95% delta-method interval | delta: contains 0 (post hoc) "
          "| delta: contains theta_0 | test-inversion 95% set | inversion: contains 0 (post hoc) "
          "| inversion: contains theta_0 (post hoc) |",
          "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    checks, notes, boot, diag, lag_rows, res = [], [], [], [], [], {}
    for hedge, label in (("H10", "10-delta"), ("H25", "25-delta"), ("Hatm", "ATM")):
        H = pm[hedge].to_numpy()
        r = delta_method(H, U)
        kind, lo, hi = theta_confidence_set(H, U)
        if kind == "interval":
            grid_note = f"[{lo:.3f}, {hi:.3f}]"
        else:
            a, b, contiguous = rejected_range(H, U)
            grid_note = (f"{kind}: the whole grid [-10, 10]" if not np.isfinite(a) else
                         f"{kind}: grid minus [{a:.3f}, {b:.3f}]" + ("" if contiguous else " (rejected points not contiguous)"))
        acc0, t0 = inversion_accepts(H, U, 0.0)
        acc_th0, t_th0 = inversion_accepts(H, U, theta0)
        is10 = hedge == "H10"
        L.append(f"| {label} | {r['theta']:.4f} | {r['se_delta']:.4f} | {r['lag']} | [{r['lo']:.4f}, {r['hi']:.4f}] | "
                 f"{r['lo'] <= 0 <= r['hi']} | {(r['lo'] <= theta0 <= r['hi']) if is10 else ''} | {grid_note} | "
                 f"{acc0} | {acc_th0 if is10 else ''} |")
        checks.append(max(abs(r["se_explicit_2x2"] / r["se_delta"] - 1), abs(r["se_via_inversion_stat"] / r["se_delta"] - 1)))

        # Test of v = 0 (post hoc). Both t statistics are the mean of H - U over a standard error: its own
        # Newey-West s.e. (test inversion) or |mean U| se_delta, the s.e. of the mean of H - (1 - theta_UB)U.
        cost = H - U
        c = mean_and_se(cost)
        ratio = c["se"] / (r["se_delta"] * abs(r["mean_U"]))
        t_delta = -np.sign(r["mean_U"]) * r["theta"] / r["se_delta"]
        note = (f"{label}: mean {BP * c['mean']:.2f} bp a month, t {t0:.2f} by test inversion and {t_delta:.2f} by the "
                f"delta method; s.e. of the mean of H - U over that of H - (1 - theta_UB)U {ratio:.3f}")
        if is10 and ratio < 1:
            note += (", below 1 because H - U, the option P&L alone, leaves out the share theta_UB U of the carry "
                     "return that H - (1 - theta_UB)U retains")
        notes.append(note)
        bt = studentised_bootstrap(cost, B)
        boot.append(f"{label}: block length {bt['block']:.2f}, 2.5% and 97.5% quantiles of t* {bt['q_lo']:.2f} and "
                    f"{bt['q_hi']:.2f}, 95% interval for the mean [{BP * bt['lo']:.1f}, {BP * bt['hi']:.1f}] bp, share "
                    f"of t* at or below the observed t of {t0:.2f}: {bt['share_below']:.3f}")
        drop = mean_and_se(np.delete(cost, np.argmax(cost)))
        t_drop, sk = drop["mean"] / drop["se"], skewness(cost)
        diag.append(f"{label}: skewness of H - U {sk:.2f}, t after dropping the month with the largest H - U {t_drop:.2f}")
        ts0 = [t_fixed_lag(cost, k) for k in FIXED_LAGS]
        lag_rows.append((f"{label}, v = 0", ts0, t0, c["lag"]))
        res[hedge] = {"t0": t0, "boot": bt, "skew": sk, "t_drop": t_drop, "n_rej0": sum(abs(t) > 1.96 for t in ts0)}
        if is10:
            x0 = H - (1 - theta0) * U
            ts_th0 = [t_fixed_lag(x0, k) for k in FIXED_LAGS]
            lag_x0 = mean_and_se(x0)["lag"]
            lag_rows.append((f"{label}, v = theta_0", ts_th0, t_th0, lag_x0))
            res[hedge].update(acc_th0=acc_th0, t_th0=t_th0, lag_th0=lag_x0,
                              lags_rej_th0=[k for k, t in zip(FIXED_LAGS, ts_th0) if abs(t) > 1.96])

    # Sign reading: H10 - U is the E3 total c1 + c2 + c3 (payoff less the premium carried to delivery).
    e3_gap = float(np.max(np.abs((pm.H10 - pm.U) - (pm.c1 + pm.c2 + pm.c3))))
    ext = s4[(s4["sample"] == "extended") & (s4.status == "ok")].sort_values("date")
    r_ext = mean_and_se((ext.H10 - ext.U).to_numpy())
    t_ext = r_ext["mean"] / r_ext["se"]

    r10 = res["H10"]
    b10 = r10["boot"]
    others_reject = any(abs(res[h]["t0"]) > 1.96 or not (res[h]["boot"]["lo"] <= 0 <= res[h]["boot"]["hi"])
                        for h in ("H25", "Hatm"))
    # The conclusion below was written for the run of 27 September 2026; it is printed only while the
    # numbers still support each of its clauses.
    as_written = (r10["t0"] < -1.96 and b10["lo"] <= 0 <= b10["hi"] and r10["skew"] > 0 and not others_reject
                  and r10["t_drop"] < r10["t0"] and abs(t_ext) < abs(r10["t0"])
                  and ext.date.min().year <= 2008 <= ext.date.max().year)
    if as_written:
        conclusion = (
            f"- Conclusion on v = 0 (post hoc): under the studentised bootstrap v = 0 is not rejected for the 10-delta "
            f"hedge at the two-sided 5% level, since its 95% interval for the mean of H - U contains 0 (the one-sided "
            f"share above is a one-sided p value), so the rejection at t = {r10['t0']:.2f} "
            f"depends on the normal reference. That reference is inaccurate here because H - U is strongly "
            f"right-skewed (skewness {r10['skew']:.2f}): for the mean of right-skewed data the t statistic has a "
            f"heavier left tail than the normal, and the 2.5% quantile of t* is {b10['q_lo']:.2f} rather than -1.96. "
            f"The 25-delta and ATM hedges reject v = 0 under neither reference. The rejection also reflects how few "
            f"large option payoffs the primary sample ({pm.date.min():%Y-%m} to {pm.date.max():%Y-%m}) contained, the "
            f"peso caveat of E3, and the direction is known: removing the month with the largest H - U strengthens it "
            f"(t {r10['t_drop']:.2f}), whereas more crash payoffs would weaken it, as over the extended sample "
            f"({ext.date.min():%Y-%m} to {ext.date.max():%Y-%m}, {len(ext)} months, which contains 2008; two legs per "
            f"side and Fenics smiles, so a different portfolio), where the t "
            f"of the mean of H10 - U is {t_ext:.2f}. Under the normal reference the rejection of v = 0 holds at "
            + ("every fixed lag from 0 to 12" if r10["n_rej0"] == len(FIXED_LAGS) else
               f"{r10['n_rej0']} of the {len(FIXED_LAGS)} fixed lags from 0 to 12")
            + " (table above): it depends on the reference distribution, not on the lag.")
    else:
        conclusion = ("- Conclusion on v = 0 (post hoc): the numbers no longer match the conclusion written for the run "
                      "of 27 September 2026 (see estimate_design_checks.delta_check); it must be revisited.")
    rej_th0 = ", ".join(map(str, r10["lags_rej_th0"])) or "none"
    L += ["", f"- Check: the explicit 2x2 Bartlett long-run covariance of (H, U) at the same lag and the inversion statistic "
          f"at theta_UB give the same standard error; largest relative gap {max(checks):.1e}.",
          "- Test of v = 0 (post hoc): H - U is the net P&L of the protective options, their payoff less the premium "
          "carried to delivery (for the 10-delta hedge the E3 total c1 + c2 + c3, which it matches to "
          f"{e3_gap:.1e}), so a negative mean means that the protection cost more than it paid. Both t statistics are the "
          "mean of H - U over a standard error: by test inversion its own Newey-West s.e., by the delta method "
          "|mean HML^U| se_delta, the s.e. of the mean of H - (1 - theta_UB)U, which gives "
          "-sign(mean HML^U) theta_UB/se_delta. " + "; ".join(notes) + ".",
          f"- Studentised (bootstrap-t) stationary bootstrap of the mean of H - U (post hoc, A2): {B:,} draws, seed "
          f"{BOOT_SEED}, expected block length from optimal_block_length, t* = (mean* - mean)/se* with the Newey-West "
          "s.e. at the automatic lag recomputed in each draw, interval [mean - q(97.5%) se, mean - q(2.5%) se]; the "
          "share of t* at or below the observed t is a one-sided bootstrap p value. " + "; ".join(boot) + ".",
          "- Skewness and influence (post hoc): " + "; ".join(diag) + ".",
          f"- theta_0 (post hoc): the 10-delta set {'retains' if r10['acc_th0'] else 'rejects'} theta_0 = {theta0:.4f} "
          f"at the automatic lag (t {r10['t_th0']:.2f}, lag {r10['lag_th0']}); the fixed lags from 0 to 12 at which the "
          f"test rejects it: {rej_th0}"
          + (". Its membership of the set is therefore fragile." if r10["acc_th0"] and r10["lags_rej_th0"] else "."),
          "",
          "Fixed-lag t statistics of the test inversion (post hoc): mean of H - (1 - v)U over its Newey-West s.e. at each "
          "Bartlett lag.", "",
          "| statistic | " + " | ".join(f"lag {k}" for k in FIXED_LAGS) + " | automatic lag |",
          "| --- | " + " | ".join("---" for _ in FIXED_LAGS) + " | --- |"]
    for name, ts, t_auto, lag_auto in lag_rows:
        L.append(f"| {name} | " + " | ".join(f"{t:.2f}" for t in ts) + f" | {t_auto:.2f} (lag {lag_auto}) |")
    L += ["", conclusion,
          "- Why the delta-method interval misleads here: write se(v) for the Newey-West s.e. of the mean of "
          "H - (1 - v)U. Since mu_H - (1 - v) mu_U = mu_U (v - theta_UB), the delta-method interval is the set where "
          "this is within 1.96 se(theta_UB) of zero (at v = theta_UB, se(v) is |mean HML^U| times the delta-method "
          "s.e. in the table, not that s.e.), while the test-inversion set uses se(v), the standard error under the "
          "value tested. For a fixed lag, se(v) is the square root of a quadratic in v with leading coefficient "
          "se(mu_U)^2 at that lag; since the automatic lag tends to that of U, se(v)/|v| tends to se(mu_U) as |v| grows, "
          "although near theta_UB se(v) need not be monotone. With |t| of mean HML^U below 1.96, every v far enough "
          "from theta_UB is retained and the test-inversion set is unbounded; the delta-method interval freezes the "
          "standard error at theta_UB and stays bounded however weakly mu_U is identified, so its finite length is an "
          "artefact of dividing by a mean that cannot be distinguished from zero. Freezing the standard error also "
          "moves the boundary between retained and rejected values: at v = 0 the tested statistic is the mean hedge cost "
          "H - U, whose standard error differs from that of H - (1 - theta_UB)U (ratio above), so the two procedures "
          "can disagree about v = 0. Procedures that always return bounded sets have zero worst-case coverage for such "
          "ratios (Gleser and Hwang, 1987; Dufour, 1997, as cited in R9(c) of the paper).", ""]
    return L


# ---------------------------------------------------------------------------
# 4. Closed-form checks of the moment code


def moment_checks(mom):
    ok = mom[mom.status == "ok"]
    L = ["## 4. Closed-form and convergence checks of the moment code (moments_1m.csv)", "",
         f"Currency-months with status ok: {len(ok)} of {len(mom)}.", ""]
    for s, label in (("i", "setting (i), boundary elasticities"), ("ii", "setting (ii), gamma = eta = 2")):
        cf = ok[f"closed_form_rel_diff_{s}"].astype(float)
        qe = ok[f"tail_quad_err_{s}"].astype(float)
        n_nonfinite = int((~np.isfinite(cf)).sum() + (~np.isfinite(qe)).sum() - 2 * ok[f"status_{s}"].eq("exponent_invalid").sum())
        cf, qe = cf[np.isfinite(cf)], qe[np.isfinite(qe)]
        L += [f"- {label}: largest closed_form_rel_diff {cf.max():.2e} (n = {len(cf)}, median {cf.median():.2e}); "
              f"largest tail quadrature error estimate {qe.max():.2e} (n = {len(qe)}, median {qe.median():.2e}); "
              f"above the tolerance of 1e-12 (post hoc): {int((qe > 1e-12).sum())}; statuses {ok[f'status_{s}'].value_counts().to_dict()}; "
              f"non-finite values outside exponent_invalid: {n_nonfinite}"]
    L += [f"- Simpson middle part: converged in {ok.simpson_converged.mean():.1%} of currency-months "
          f"({int(ok.simpson_converged.sum())} of {len(ok)}); subintervals per side median {ok.n_simpson.median():.0f}, "
          f"max {ok.n_simpson.max():.0f}; largest relative change at the last halving {ok.simpson_rel_change.max():.2e}", "",
          "What each checks:", "",
          "- closed_form_rel_diff: the adaptive tail quadrature of qef.fx.moments.tail_bounds against the closed-form R4 "
          "bounds for the weights K^-2 and the logarithmic weights, at each currency-month's own boundary prices and "
          "exponents (the upper logarithmic bound only where F < K_max <= eF, its stated domain).",
          "- tail_quad_err: the adaptive quadrature's own error estimate relative to its integral, the largest over the "
          "three contracts, both tails and each piece of constant sign. Being quad's own estimate, it confirms only that "
          "no tail integral stopped at quad's subdivision limit before meeting its relative tolerance of 1e-12; "
          "closed_form_rel_diff is the independent check of accuracy.",
          "- simpson_converged: whether halving the Simpson step between the 10-delta strikes changed every middle contract "
          "by less than 1e-10 relative, the stopping rule of the middle integral.", ""]
    return L


# ---------------------------------------------------------------------------


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--bootstrap", type=int, default=9999)
    p.add_argument("--workers", type=int, default=8)
    args = p.parse_args()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    raw = ROOT / "data" / "private" / "lseg" / args.retrieval_date / "raw"
    audit = ROOT / "data" / "private" / "audit" / args.retrieval_date
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    inputs, panel = load("smile_inputs.csv"), load("smile_panel.csv")
    fn_in, fn_panel = load("smile_inputs_fn.csv"), load("smile_panel_fn.csv")
    e1 = load("e1_series.csv")
    e1_pub = e1[(e1["sample"] == "primary") & (e1.reading == "market") & (e1.delta == DELTA)].sort_values("date")
    month_ends = ny_month_ends(PRIMARY_START, SAMPLE_END)

    lines = ["# Design checks listed as not reported (restricted)", "",
             "Written by scripts/estimate_design_checks.py. Pooled statistics only; the rules are in the script's "
             "docstring (R1 to R12, amended on 27 September 2026 after an independent review, A1 and A2); items added "
             "after a first run or after the review are marked post hoc.", ""]
    lines += sofr_check(inputs, panel, e1_pub, raw, month_ends, args.bootstrap, args.workers, d)
    lines += tifo_check(inputs, panel, fn_in, fn_panel, e1_pub, raw, month_ends, audit, args.workers, d)
    lines += delta_check(load("stage4_months.csv"), args.bootstrap)
    lines += moment_checks(load("moments_1m.csv"))
    (d / "design_checks_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
