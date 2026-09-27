"""Out-of-sample test of the tail hypotheses with one-month 5Δ quotes, post hoc.

The plan was recorded in the research log on 27 September 2026 before the 5Δ
history was retrieved; this script implements it. Data: the Fenics one-month 5Δ
risk reversals and butterflies (<CCY>1MR5=FN, <CCY>1MB5=FN) of the retrieval
of 27 September 2026, sampled at New York month-ends under the month-end rule
of the research design (section 4), from July 2022 to August 2026. The
boundary is the Fenics market-reading SABR smile (smile_panel_fn.csv, status
ok) with 10Δ strikes in each pair's convention, as in the moment analysis.

Conversion: σ_5P = σ_ATM + BF5 − RR5/2 and σ_5C = σ_ATM + BF5 + RR5/2 (the
smile-strangle reading, risk reversal as call minus put), with σ_ATM the
Fenics ATM mid; strikes from the flat-volatility delta inversion of R5 in each
pair's convention; undiscounted Garman–Kohlhagen out-of-the-money prices.

For each currency-month and tail (qef.fx.identification): whether the 5Δ price
lies in its identified interval under the weakest hypothesis of the class;
the largest exponent at which it does, as a fraction θ_5 of the maximal one;
the sharp identified set of the variance and third central moment with both
5Δ prices as equality constraints at the largest exponents they allow; and the
breakdown point of the sign of the skewness with the 5Δ prices as constraints,
along γ = θγ̄, η = 1 + θ(η̄ − 1) for θ up to min(θ_5 lower, θ_5 upper). Also
recorded: the sign agreement of RR5 with the Fenics RR10, and, as a
supplementary comparison not in the plan, the breakdown point without the 5Δ
prices on the same smiles and the boundary check of setting (i).

Writes wing_test_1m.csv and wing_test_summary.md to
data/private/results/2026-09-27/; only pooled statistics may be published.
"""

from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from qef.data.panel import ny_month_ends, read_raw, sample_month_ends
from qef.fx.conventions import G10
from qef.fx.gk import CALL, PUT, forward_premium, strike_from_delta, strike_from_delta_smile
from qef.fx.identification import SharpMomentSet, boundary, max_consistent_exponent, skewness_sign_breakdown
from qef.fx.moments import middle_contracts
from qef.fx.sabr import sabr_vol

ROOT = Path(__file__).resolve().parents[1]
START, END = "2022-07-01", "2026-08-31"
WINDOW = 5  # business days for month-end substitution (research design, section 4)
THETA_MIN, TOL = 0.01, 0.01


def sign_identified(S):
    lo, hi = S.third_central_range()
    return bool(hi < 0 or lo > 0), lo, hi


def one(task):
    r, (alpha, rho, nu), rr5, bf5 = task
    ccy = r["currency"]
    F, tau, df_base, atm, usd = r["F"], r["tau"], r["df_base"], r["atm"], G10[ccy].usd_base
    conv = G10[ccy].delta
    vol = lambda K: sabr_vol(K, F, tau, alpha, rho, nu, 1.0)
    out = {"date": r["date"], "currency": ccy, "rr5_rr10_same_sign": bool(np.sign(rr5) == np.sign(r["rr10"]))}
    K_min = strike_from_delta_smile(-0.10, F, tau, PUT, conv, vol, df_base)
    K_max = strike_from_delta_smile(0.10, F, tau, CALL, conv, vol, df_base)
    bd = boundary(F, tau, vol, K_min, K_max)
    s5p, s5c = atm + bf5 - 0.5 * rr5, atm + bf5 + 0.5 * rr5
    K5p = strike_from_delta(-0.05, F, s5p, tau, PUT, conv, df_base)
    K5c = strike_from_delta(0.05, F, s5c, tau, CALL, conv, df_base)
    p5, c5 = float(forward_premium(F, K5p, s5p, tau, PUT)), float(forward_premium(F, K5c, s5c, tau, CALL))
    out.update({"beyond_lower": K5p < K_min, "beyond_upper": K5c > K_max, "gamma_bar": bd.gamma_bar, "eta_bar": bd.eta_bar,
                "smile_p5": float(forward_premium(F, K5p, float(vol(K5p)), tau, PUT)), "quote_p5": p5,
                "smile_c5": float(forward_premium(F, K5c, float(vol(K5c)), tau, CALL)), "quote_c5": c5})
    # Supplementary, not in the plan: the breakdown point without the 5Δ prices, on the same smile,
    # and whether the boundary check refutes setting (i).
    out["theta_star_no5d"] = skewness_sign_breakdown(F, tau, vol, usd, bd)[0]
    out["setting_i_refuted"] = bool(bd.gamma_i > bd.gamma_bar or bd.eta_i > bd.eta_bar)
    if not (out["beyond_lower"] and out["beyond_upper"]):
        return {**out, "status": "5d_strike_not_beyond_10d"}
    g5, th_lo, ok_lo = max_consistent_exponent(bd, K5p, p5, "lower")
    e5, th_hi, ok_hi = max_consistent_exponent(bd, K5c, c5, "upper")
    out.update({"consistent_lower": ok_lo, "consistent_upper": ok_hi, "gamma_5": g5, "eta_5": e5,
                "theta5_lower": th_lo, "theta5_upper": th_hi})
    if not (ok_lo and ok_hi) or not (g5 > 1e-6 and e5 > 1 + 1e-6):
        return {**out, "status": "inconsistent_or_degenerate"}
    mid = middle_contracts(F, K_min, K_max, tau, vol, usd)["values"]
    tie = dict(put_beyond=(K5p, p5), call_beyond=(K5c, c5), mid=mid)
    g, e = min(0.999 * g5, bd.gamma_bar), 1.0 + 0.999 * (e5 - 1.0)
    ident, lo, hi = sign_identified(SharpMomentSet(F, tau, vol, usd, bd, g, e, **tie))
    out.update({"sign_identified_at_max": ident, "k3_lo_at_max": lo, "k3_hi_at_max": hi})
    # Breakdown along the common path, with the 5Δ prices as constraints.
    th_max = 0.999 * min(th_lo, th_hi)
    at = lambda th: sign_identified(SharpMomentSet(F, tau, vol, usd, bd, *bd.exponents(th), **tie))[0]
    if th_max <= THETA_MIN or not at(th_max):
        star = np.nan
    elif at(THETA_MIN):
        star = THETA_MIN
    else:
        a, b = THETA_MIN, th_max
        while b - a > TOL:
            m = 0.5 * (a + b)
            a, b = (a, m) if at(m) else (m, b)
        star = b
    out.update({"theta_max_common": th_max, "theta_star_5": star, "status": "ok"})
    return out


def q(s):
    s = pd.Series(s).dropna()
    return "{:.3g}, {:.3g}, {:.3g}".format(*np.percentile(s, [25, 50, 75])) if len(s) else "none"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--workers", type=int, default=os.cpu_count())
    args = p.parse_args()
    res_dir = ROOT / "data" / "private" / "results" / "2026-09-23"
    raw = ROOT / "data" / "private" / "lseg" / "2026-09-27" / "raw"
    out_dir = ROOT / "data" / "private" / "results" / "2026-09-27"
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs = pd.read_csv(res_dir / "smile_inputs_fn.csv", parse_dates=["date"]).set_index(["currency", "date"])
    panel = pd.read_csv(res_dir / "smile_panel_fn.csv", parse_dates=["date"])
    fit = panel[(panel.reading == "market") & (panel.status == "ok")].set_index(["currency", "date"])
    month_ends = ny_month_ends(START, END)
    tasks, counts = [], {"month_ends": 0, "missing_quote": 0, "no_smile": 0}
    for ccy in G10:
        q5 = {code: sample_month_ends(read_raw(raw / f"vol_{name}" / f"{ccy}1M{code}=FN.csv"), month_ends, WINDOW)
              for code, name in (("R5", "rr5"), ("B5", "bf5"))}
        for t in month_ends:
            counts["month_ends"] += 1
            r5, b5 = q5["R5"].loc[t], q5["B5"].loc[t]
            if r5["status"] == "missing" or b5["status"] == "missing":
                counts["missing_quote"] += 1
                continue
            if (ccy, t) not in fit.index or (ccy, t) not in inputs.index:
                counts["no_smile"] += 1
                continue
            r = inputs.loc[(ccy, t)]
            row = {"date": t, "currency": ccy, "F": r.F, "tau": r.tau, "df_base": r.df_base, "atm": r.atm, "rr10": r.rr10}
            mid = lambda x: 0.5 * (x["bid"] + x["ask"]) / 100.0
            tasks.append((row, tuple(fit.loc[(ccy, t), ["alpha", "rho", "nu"]]), mid(r5), mid(b5)))
    t0 = time.time()
    with ProcessPoolExecutor(args.workers) as ex:
        res = pd.DataFrame(list(ex.map(one, tasks, chunksize=2)))
    res.to_csv(out_dir / "wing_test_1m.csv", index=False)
    ok = res[res.status == "ok"]
    cons = res[res.status.isin(["ok", "inconsistent_or_degenerate"])]
    L = ["# Out-of-sample test of the tail hypotheses with one-month 5Δ quotes (restricted; post hoc)", "",
         f"Month-ends {START[:7]} to {END[:7]}: {counts['month_ends']} currency-month-ends; missing 5Δ quote {counts['missing_quote']}; "
         f"no converged Fenics smile {counts['no_smile']}; tested {len(res)}; run time {time.time() - t0:.0f} s.", "",
         f"- RR5 has the sign of the Fenics RR10 in {res.rr5_rr10_same_sign.mean():.1%} of tested currency-months",
         f"- 5Δ strikes beyond the smile's 10Δ strikes: lower {res.beyond_lower.mean():.1%}, upper {res.beyond_upper.mean():.1%}",
         f"- 5Δ price relative to the smile's own extrapolation (quote/smile): put quartiles {q(res.quote_p5 / res.smile_p5)}; "
         f"call quartiles {q(res.quote_c5 / res.smile_c5)}",
         f"- consistent with the class under the weakest hypothesis: lower {cons.consistent_lower.mean():.1%}, upper "
         f"{cons.consistent_upper.mean():.1%}, both {(cons.consistent_lower & cons.consistent_upper).mean():.1%} (of {len(cons)})",
         f"- θ_5, the largest consistent fraction of the maximal exponent: lower quartiles {q(cons.theta5_lower)}, "
         f"upper quartiles {q(cons.theta5_upper)}",
         f"- share of consistent tails at the maximal exponent (θ_5 = 1): lower {(cons.theta5_lower >= 1 - 1e-9).mean():.1%}, "
         f"upper {(cons.theta5_upper >= 1 - 1e-9).mean():.1%}",
         f"- sign of skewness identified at the largest exponents the 5Δ prices allow, with the 5Δ prices as constraints: "
         f"{ok.sign_identified_at_max.mean():.1%} of {len(ok)}; negative {(ok.k3_hi_at_max < 0).mean():.1%}",
         f"- breakdown θ*_5 with the 5Δ constraints, where identified: quartiles {q(ok.theta_star_5)}; "
         f"identified at θ ≤ 0.5: {(ok.theta_star_5 <= 0.5).mean():.1%} of {len(ok)}",
         f"- sign of skewness identified along the common path at θ = 0.999 min(θ_5 lower, θ_5 upper): "
         f"{ok.theta_star_5.notna().mean():.1%} of {len(ok)}",
         "", "Supplementary comparison, not in the plan (same Fenics smiles and currency-months):", "",
         f"- breakdown θ* without the 5Δ prices: quartiles {q(res.theta_star_no5d)}",
         f"- where the sign is identified with the 5Δ prices, they lower θ* in "
         f"{(ok.theta_star_5 < ok.theta_star_no5d - 1e-9)[ok.theta_star_5.notna()].mean():.1%} of currency-months, by a median of "
         f"{(ok.theta_star_no5d - ok.theta_star_5).median():.3f}",
         f"- setting (i) refuted by the boundary check: {res.setting_i_refuted.mean():.1%}",
         "", "## By currency", ""]
    for c, g in res.groupby("currency"):
        gok = g[g.status == "ok"]
        L.append(f"- {c}: tested {len(g)}; consistent both tails {(g.consistent_lower.fillna(False) & g.consistent_upper.fillna(False)).mean():.1%}; "
                 f"median θ_5 lower {pd.Series(g.theta5_lower).median():.2f}, upper {pd.Series(g.theta5_upper).median():.2f}; "
                 f"sign identified at max {gok.sign_identified_at_max.mean() if len(gok) else float('nan'):.1%}")
    (out_dir / "wing_test_summary.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
