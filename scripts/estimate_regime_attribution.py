"""Regime attribution of the change in the carry spread, post hoc.

Research log, 27 September 2026. Added after the E1 to E5 results were
recorded, and reported as post hoc. Primary sample, market reading, 10Δ
hedges. Regimes as in E1: zero-rate from the primary start to December 2021,
hiking from January 2022.

For the skew price C_skew, the forward-discount spread FD (both from
e1_series.csv, bp per month) and the oriented 10Δ risk reversal of the carry
portfolio (stage4_months.csv; the quoted volatility of the protective option
minus that of the opposite option, averaged over legs), it reports the regime
means and their difference, with the Newey–West standard error of the
coefficient on a regime dummy, an AR(1) plug-in standard error, and the
stationary-bootstrap 95% percentile interval (each regime resampled separately
with its own Politis–White block length; 9,999 draws), as for E1.

The ratio of regime differences r = ΔC_skew/ΔFD is the grouping estimator of
Wald (1940) with the regime as the grouping variable. Its 95% set is found by
test inversion, as for θ_UB: r is kept if the regime difference of
C_skew − r FD is not significant at 5%, with (a) the Newey–West standard error
and (b) the stationary bootstrap, resampling the pair (C_skew, FD) jointly
within each regime with the larger of the two block lengths. This is Fieller's
(1954) construction, which with one binary instrument coincides with the set of
Anderson and Rubin (1949). The one-for-one hypothesis r = 1 is tested in the
same way.

Writes regime_attribution.csv and regime_attribution_summary.md to
data/private/results/<retrieval date>/; they contain aggregate statistics only.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from estimate_e1 import REGIME_BREAK, regime_summary  # noqa: E402

from qef.stats.bootstrap import optimal_block_length, stationary_indices  # noqa: E402
from qef.stats.hac import ols_hac  # noqa: E402

GRID = np.linspace(-3.0, 3.0, 6001)
SEED = 20260926


def ar1_se(x):
    x = np.asarray(x, dtype=float)
    d = x - x.mean()
    rho = float(d[1:] @ d[:-1] / (d @ d))
    return np.sqrt(np.var(x) * (1 + rho) / (1 - rho) / len(x)), rho


def dummy_diff(y, hiking):
    reg = ols_hac(y, np.column_stack([np.ones(len(y)), hiking.astype(float)]))
    return reg["beta"][1], reg["se"][1]


def joint_bootstrap(c0, f0, c1, f1, B):
    """Draws of (ΔC, ΔFD), resampling the pairs jointly within each regime."""
    rng = np.random.default_rng(SEED)
    b0 = max(optimal_block_length(c0), optimal_block_length(f0))
    b1 = max(optimal_block_length(c1), optimal_block_length(f1))
    out = np.empty((B, 2))
    for r in range(B):
        i0, i1 = stationary_indices(len(c0), b0, rng), stationary_indices(len(c1), b1, rng)
        out[r] = c1[i1].mean() - c0[i0].mean(), f1[i1].mean() - f0[i0].mean()
    return out, (b0, b1)


def accepted_set(accept):
    inside = GRID[accept]
    if not accept.any():
        return "empty", np.nan, np.nan
    if accept[0] or accept[-1]:
        return "unbounded within the grid", float(inside.min()), float(inside.max())
    contiguous = np.all(np.diff(np.where(accept)[0]) == 1)
    return ("interval" if contiguous else "union of intervals"), float(inside.min()), float(inside.max())


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--bootstrap", type=int, default=9999)
    args = p.parse_args()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    e1 = pd.read_csv(d / "e1_series.csv", parse_dates=["date"])
    e1 = e1[(e1["sample"] == "primary") & (e1.reading == "market") & (e1.delta == 0.10) & (e1.status == "ok")]
    e1 = e1.sort_values("date").assign(C_skew_bp=lambda x: 1e4 * x.C_skew, FD_bp=lambda x: 1e4 * x.FD)
    s4 = pd.read_csv(d / "stage4_months.csv", parse_dates=["date"])
    s4 = s4[(s4["sample"] == "primary") & (s4.status == "ok")].sort_values("date").assign(rr_or_vol=lambda x: 100 * x.rr_or)

    rows = []
    for frame, col in ((e1, "C_skew_bp"), (e1, "FD_bp"), (s4, "rr_or_vol")):
        s = frame.assign(status="ok")
        r = regime_summary(s, col, args.bootstrap)
        a, b = s.loc[s.date < REGIME_BREAK, col].to_numpy(), s.loc[s.date >= REGIME_BREAK, col].to_numpy()
        (se_a, rho_a), (se_b, rho_b) = ar1_se(a), ar1_se(b)
        r.update(se_difference_ar1=float(np.hypot(se_a, se_b)), ar1_zero_rate=rho_a, ar1_hiking=rho_b)
        rows.append(r)
    table = pd.DataFrame(rows)

    hiking = (e1.date >= REGIME_BREAK).to_numpy()
    C, FD = e1.C_skew_bp.to_numpy(), e1.FD_bp.to_numpy()
    dC, dF = C[hiking].mean() - C[~hiking].mean(), FD[hiking].mean() - FD[~hiking].mean()
    ratio = dC / dF
    nw = np.array([abs(dummy_diff(C - g * FD, hiking)[0]) <= 1.96 * dummy_diff(C - g * FD, hiking)[1] for g in GRID])
    draws, blocks = joint_bootstrap(C[~hiking], FD[~hiking], C[hiking], FD[hiking], args.bootstrap)
    z = draws[:, [0]] - GRID[None, :] * draws[:, [1]]
    lo, hi = np.percentile(z, [2.5, 97.5], axis=0)
    boot = (lo <= 0) & (0 <= hi)
    one = dummy_diff(C - FD, hiking)
    z1 = draws[:, 0] - draws[:, 1]
    p_one = float(2 * min((z1 >= 0).mean(), (z1 <= 0).mean()))

    kind_nw, lo_nw, hi_nw = accepted_set(nw)
    kind_b, lo_b, hi_b = accepted_set(boot)
    extra = pd.DataFrame([{
        "ratio_dC_over_dFD": ratio, "set_nw": kind_nw, "set_nw_lo": lo_nw, "set_nw_hi": hi_nw,
        "set_boot": kind_b, "set_boot_lo": lo_b, "set_boot_hi": hi_b, "boot_blocks": blocks,
        "one_for_one_diff_bp": one[0], "one_for_one_se_nw": one[1], "one_for_one_t_nw": one[0] / one[1],
        "one_for_one_p_boot_two_sided": p_one, "corr_C_FD_months": float(np.corrcoef(C, FD)[0, 1]),
        "n_zero_rate": int((~hiking).sum()), "n_hiking": int(hiking.sum())}])
    table.to_csv(d / "regime_attribution.csv", index=False)
    extra.to_csv(d / "regime_attribution_ratio.csv", index=False)
    fmt = lambda v: f"{v:.6g}" if isinstance(v, (float, np.floating)) else str(v)
    lines = ["# Regime attribution of the change in the carry spread (restricted; post hoc)", "",
             f"Primary sample, market reading, 10Δ; bootstrap draws {args.bootstrap}.", "",
             "## Regime differences (hiking minus zero-rate)", "", "```",
             table.to_string(index=False, float_format=fmt), "```", "",
             "## Ratio of regime differences ΔC_skew/ΔFD (Wald grouping estimator) and test-inversion sets", "", "```",
             extra.T.to_string(header=False, float_format=fmt), "```", ""]
    (d / "regime_attribution_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
