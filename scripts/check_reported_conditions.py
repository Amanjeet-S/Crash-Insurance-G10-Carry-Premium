"""Checks behind statements in the paper that no estimation script produced.

A check of the project on 3 October 2026 (research log) found three groups of
reported statements without a script behind them. This script computes them
from the stored calibrations and results, under rules fixed before it was run.

1. Hypotheses of R5(c) and R8 on the calibrated one-month smiles of the
   primary sample (both butterfly readings, converged smiles):
   - the R5(c) slope statistic |∂σ/∂k| √τ |d−(K, σ(K))|, with k = ln(K/F),
     on 801 log-moneyness points over ±4 ATM standard deviations, the slope by
     central differences with step 1e-5 in k; the condition holds where the
     statistic is below 1;
   - for the five dollar-base pairs, whose convention is the premium-adjusted
     delta, whether the smile's premium-adjusted put delta is strictly monotone
     on the same grid and the call delta strictly decreasing to the right of its
     maximum on the grid;
   - the largest condition number of the calibration Jacobian, as stored by the
     calibration.

2. The grid of means in the identified sets (R10). SharpMomentSet._extreme
   searches the feasible means on a grid of 25 points (15 in the breakdown
   search), refined around the two best. Every mean it evaluates belongs to a
   feasible law, so the computed ranges are inner approximations, which can
   only err towards declaring the sign of the skewness identified. On the same
   40 currency-months as the step-point check of estimate_identification.py,
   the third-central-moment range is recomputed with 200 means refined around
   the 5 best (setting (i) where admissible, and θ = 0.5, 0.75 and 0.9), and θ*
   with 120 means; the script reports the largest widening of a range relative
   to its width, any change in whether the sign is identified, and the largest
   change in θ*.

3. The breakdown point with and without the 5Δ prices on the same
   currency-months: the quartiles of θ* without the 5Δ prices over the
   currency-months where θ* with them is defined.

Outputs, aggregate only, go to data/private/results/2026-09-23/:
reported_conditions_summary.md.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from estimate_identification import GRID_CHECK, THETAS  # noqa: E402

from qef.fx.conventions import G10  # noqa: E402
from qef.fx.gk import d_plus_minus, delta  # noqa: E402
from qef.fx.identification import SharpMomentSet, boundary, skewness_sign_breakdown  # noqa: E402
from qef.fx.moments import middle_contracts  # noqa: E402
from qef.fx.sabr import sabr_vol  # noqa: E402

PRIMARY_START = pd.Timestamp("2013-05-31")
N_GRID, WIDTH, H = 801, 4.0, 1e-5
N_T_FINE, N_REFINE_FINE, N_T_BREAK_FINE = 200, 5, 120


def smile_conditions(row) -> dict:
    vol = lambda K: sabr_vol(K, row.F, row.tau, row.alpha, row.rho, row.nu, 1.0)
    sd = row.atm * np.sqrt(row.tau)
    k = np.linspace(-WIDTH * sd, WIDTH * sd, N_GRID)
    K = row.F * np.exp(k)
    sig = np.asarray(vol(K), dtype=float)
    dsig = (np.asarray(vol(row.F * np.exp(k + H)), dtype=float) - np.asarray(vol(row.F * np.exp(k - H)), dtype=float)) / (2 * H)
    _, dm = d_plus_minus(row.F, K, sig, row.tau)
    out = {"r5c_max": float(np.max(np.abs(dsig) * np.sqrt(row.tau) * np.abs(dm)))}
    conv = G10[row.currency].delta
    if conv.premium_adjusted:
        put = np.asarray(delta(row.F, K, sig, row.tau, -1, conv, row.df_base), dtype=float)
        call = np.asarray(delta(row.F, K, sig, row.tau, 1, conv, row.df_base), dtype=float)
        dp = np.diff(put)
        out["pa_put_monotone"] = bool(np.all(dp < 0) or np.all(dp > 0))
        i = int(np.argmax(call))
        out["pa_call_decreasing_right"] = bool(np.all(np.diff(call[i:]) < 0))
    return out


def fine_identification(task):
    r, (alpha, rho, nu) = task
    vol = lambda K: sabr_vol(K, r["F"], r["tau"], alpha, rho, nu, 1.0)
    F, tau, usd = r["F"], r["tau"], bool(r["usd_base"])
    bd = boundary(F, tau, vol, r["K_min"], r["K_max"])
    mid = middle_contracts(F, bd.K_min, bd.K_max, tau, vol, usd)["values"]
    out = {"date": r["date"], "currency": r["currency"]}
    settings = {"i": (bd.gamma_i, bd.eta_i), **{f"t{int(round(100 * t)):02d}": bd.exponents(t) for t in THETAS}}
    for name, (g, e) in settings.items():
        if not (0 < g <= bd.gamma_bar and 1 < e <= bd.eta_bar):
            continue
        S = SharpMomentSet(F, tau, vol, usd, bd, g, e, n_c=1000, mid=mid)
        # _extreme with sign +1 returns the minimum over the means and with sign -1 the maximum (third_central_range)
        lo = S._extreme(S.third_central, 1.0, N_T_FINE, n_refine=N_REFINE_FINE)
        hi = S._extreme(S.third_central, -1.0, N_T_FINE, n_refine=N_REFINE_FINE)
        out.update({f"k3_lo_{name}": lo, f"k3_hi_{name}": hi})
    out["theta_star"], _ = skewness_sign_breakdown(F, tau, vol, usd, bd, n_c=1000, n_t=N_T_BREAK_FINE)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--workers", type=int, default=os.cpu_count())
    args = p.parse_args()
    t0 = time.time()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    lines = ["# Checks behind reported statements (restricted)", ""]

    # 1. Hypotheses of R5(c) and R8
    inputs, panel = load("smile_inputs.csv"), load("smile_panel.csv")
    fits = panel[(panel.date >= PRIMARY_START) & (panel.status == "ok")].merge(
        inputs[["date", "currency", "F", "tau", "atm", "df_base"]], on=["date", "currency"])
    cond = pd.DataFrame([{**smile_conditions(r), "reading": r.reading, "currency": r.currency} for r in fits.itertuples()])
    lines += ["## 1. Hypotheses of R5(c) and R8 (primary one-month smiles, converged)", ""]
    for reading, g in cond.groupby("reading"):
        pa = g.dropna(subset=["pa_put_monotone"])
        jc = fits[fits.reading == reading].jacobian_cond
        lines += [f"- {reading} reading: {len(g)} smiles; R5(c) statistic over ±4 ATM s.d.: maximum {g.r5c_max.max():.3f}, "
                  f"median of the per-smile maxima {g.r5c_max.median():.3f}, smiles with the condition failing somewhere "
                  f"{int((g.r5c_max >= 1).sum())}; premium-adjusted pairs: {len(pa)} smiles, put delta strictly monotone in "
                  f"{int(pa.pa_put_monotone.sum())}, call delta strictly decreasing right of its maximum in "
                  f"{int(pa.pa_call_decreasing_right.sum())}; Jacobian condition number: maximum {jc.max():.1f}, median {jc.median():.1f}."]
    lines.append("")

    # 2. Grid of means in the identified sets
    ident = load("identification_1m.csv")
    mom = load("moments_1m.csv")
    mom = mom[mom.status == "ok"].sort_values(["date", "currency"])
    fit_m = panel[panel.reading == "market"].set_index(["currency", "date"])
    cols = ["date", "currency", "F", "tau", "usd_base", "K_min", "K_max"]
    tasks = [(r, tuple(fit_m.loc[(r["currency"], r["date"]), ["alpha", "rho", "nu"]])) for r in mom[cols].to_dict("records")]
    sub = tasks[:: max(1, len(tasks) // GRID_CHECK)][:GRID_CHECK]
    with ProcessPoolExecutor(args.workers) as ex:
        fine = pd.DataFrame(list(ex.map(fine_identification, sub, chunksize=1)))
    m = fine.merge(ident, on=["date", "currency"], suffixes=("_fine", ""))
    widen, narrow, flips, n_cmp = 0.0, 0.0, 0, 0
    for name in ["i"] + [f"t{int(round(100 * t)):02d}" for t in THETAS]:
        lo_c, hi_c, lo_f, hi_f = f"k3_lo_{name}", f"k3_hi_{name}", f"k3_lo_{name}_fine", f"k3_hi_{name}_fine"
        if lo_f not in m:
            continue
        ok = m[lo_f].notna() & m[lo_c].notna()
        x = m[ok]
        width = (x[hi_c] - x[lo_c]).abs().replace(0, np.nan)
        widen = max(widen, float(np.nanmax(np.maximum(x[lo_c] - x[lo_f], x[hi_f] - x[hi_c]) / width)))
        narrow = min(narrow, float(np.nanmin(np.minimum(x[lo_c] - x[lo_f], x[hi_f] - x[hi_c]) / width)))
        coarse_id = (x[hi_c] < 0) | (x[lo_c] > 0)
        fine_id = (x[hi_f] < 0) | (x[lo_f] > 0)
        flips += int((coarse_id != fine_id).sum())
        n_cmp += int(ok.sum())
    dtheta = float((m.theta_star_fine - m.theta_star).abs().max())
    lines += ["## 2. Grid of means in the identified sets", "",
              f"- {len(m)} currency-months (the step-point check's subset), {n_cmp} ranges compared; means: {N_T_FINE} "
              f"refined around the {N_REFINE_FINE} best, against 25 refined around 2; breakdown search {N_T_BREAK_FINE} means against 15.",
              f"- Largest widening of a third-central-moment range, relative to its width: {widen:.2e} (largest narrowing {-narrow:.2e}); "
              f"ranges whose identification of the sign changes: {flips}; largest change in theta*: {dtheta:.3f}.", ""]

    # 3. Breakdown with and without the 5-delta prices on the same months
    w = pd.read_csv(ROOT / "data" / "private" / "results" / "2026-09-27" / "wing_test_1m.csv", parse_dates=["date"])
    same = w[w.theta_star_5.notna()]
    comp = same.merge(ident[["date", "currency", "theta_star"]], on=["date", "currency"], how="left")
    q = lambda s: ", ".join(f"{v:.2f}" for v in s.quantile([0.25, 0.5, 0.75]))
    lines += ["## 3. Breakdown point with and without the 5-delta prices", "",
              f"- Currency-months with theta* under the 5-delta constraints: {len(same)}; quartiles with the constraints "
              f"{q(same.theta_star_5)}; without them on the same months {q(same.theta_star_no5d)}; "
              f"without them on all {int(w.theta_star_no5d.notna().sum())} months {q(w.theta_star_no5d.dropna())}; "
              f"composite smiles on the same months ({int(comp.theta_star.notna().sum())}) {q(comp.theta_star.dropna())}; "
              f"median reduction by the 5-delta prices {float((same.theta_star_no5d - same.theta_star_5).median()):.2f}, "
              f"lowered in {float((same.theta_star_no5d > same.theta_star_5).mean()):.1%} of them.", "",
              f"Run time {time.time() - t0:.0f} s.", ""]
    (d / "reported_conditions_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
