"""Sharp identification of one-month option-implied moments (R10), post hoc.

Research log, 26 September 2026. This analysis was added after the moment
intervals of estimate_moments.py were recorded, and it is reported as post hoc.
For every currency-month of moments_1m.csv with status ok (primary sample,
market reading, SABR smile with status ok or not_converged), with the same 10Δ
strikes, it computes from the smile at the boundary strikes:

- γ̄ = K_min G0/P0 − 1 and η̄ = 1 + K_max Ḡ0/C0, the largest exponents
  compatible with the quotes (R10(a)), and the setting (i) exponents γ_i, η_i;
  setting (i) is refuted by the quotes in a tail where γ_i > γ̄ or η_i > η̄;
- the sharp identified ranges of the variance and of the third central moment
  of y = ln(X_T/F_X) under Q^USD (qef.fx.identification) at γ = η = 2
  (setting (ii), where admissible), at the setting (i) exponents (where
  admissible), and at θ = 0.25, 0.5, 0.75 and 0.9 of the maximal exponents,
  γ = θ γ̄ and η = 1 + θ(η̄ − 1);
- the moments at θ = 1, where the tails are pure power laws and the moments
  are point-identified;
- the breakdown value θ*, the smallest θ at which the sign of the skewness is
  identified, to within 0.01 and not below 0.01.

A grid check reruns a subset of currency-months with four times as many step
points. Writes identification_1m.csv (one row per currency-month) and
identification_summary.md (pooled statistics) to data/private/results/; both
are LSEG-derived, and only pooled statistics may be published.
"""

from __future__ import annotations

import os

# One thread per worker: the workers run in parallel, and multithreaded BLAS in each
# of them oversubscribes the cores.
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from qef.fx.identification import SharpMomentSet, boundary, pareto_moments, skewness_sign_breakdown
from qef.fx.moments import middle_contracts
from qef.fx.sabr import sabr_vol

ROOT = Path(__file__).resolve().parents[1]
THETAS = (0.25, 0.5, 0.75, 0.9)
GRID_CHECK = 40  # currency-months rerun with four times as many step points


def third_central(m):
    return m[2] - 3 * m[0] * m[1] + 2 * m[0] ** 3


def one(task):
    r, (alpha, rho, nu), n_c = task
    vol = lambda K: sabr_vol(K, r["F"], r["tau"], alpha, rho, nu, 1.0)
    F, tau, usd = r["F"], r["tau"], bool(r["usd_base"])
    bd = boundary(F, tau, vol, r["K_min"], r["K_max"])
    out = {"date": r["date"], "currency": r["currency"], "gamma_bar": bd.gamma_bar, "eta_bar": bd.eta_bar,
           "gamma_i": bd.gamma_i, "eta_i": bd.eta_i,
           "refuted_lower": bd.gamma_i > bd.gamma_bar, "refuted_upper": bd.eta_i > bd.eta_bar}
    mid = middle_contracts(F, bd.K_min, bd.K_max, tau, vol, usd)["values"]
    settings = {"ii": (2.0, 2.0), "i": (bd.gamma_i, bd.eta_i),
                **{f"t{int(round(100 * t)):02d}": bd.exponents(t) for t in THETAS}}
    for name, (g, e) in settings.items():
        if not (0 < g <= bd.gamma_bar and 1 < e <= bd.eta_bar):
            out[f"status_{name}"] = "not_admissible"
            continue
        S = SharpMomentSet(F, tau, vol, usd, bd, g, e, n_c=n_c, mid=mid)
        v, k3 = S.variance_range(), S.third_central_range()
        out.update({f"var_lo_{name}": v[0], f"var_hi_{name}": v[1], f"k3_lo_{name}": k3[0], f"k3_hi_{name}": k3[1],
                    f"sign_identified_{name}": bool(k3[1] < 0 or k3[0] > 0), f"status_{name}": "ok"})
    m = pareto_moments(F, tau, vol, usd, bd, mid)
    var_p = m[1] - m[0] ** 2
    out.update({"var_pareto": var_p, "k3_pareto": third_central(m), "skew_pareto": third_central(m) / var_p ** 1.5})
    out["theta_star"], out["sign_star"] = skewness_sign_breakdown(F, tau, vol, usd, bd, n_c=n_c)
    return out


def summary(res: pd.DataFrame, mom: pd.DataFrame, check: pd.DataFrame | None, seconds: float) -> list[str]:
    q = lambda s: "{:.3g}, {:.3g}, {:.3g}".format(*np.percentile(s.dropna(), [25, 50, 75]))
    n = len(res)
    r4 = mom[["date", "currency", "smile_tail_ok_lower_i", "smile_tail_ok_upper_i", "status_i", "skew_lo_i",
              "skew_hi_i", "skew_lo_ii", "skew_hi_ii"]]
    j = res.merge(r4.rename(columns={c: f"r4_{c}" for c in r4.columns[2:]}), on=["date", "currency"], how="left")
    L = ["# Sharp identification of one-month option-implied moments (restricted; post hoc)", "",
         f"Currency-months: {n} ({res.date.min():%Y-%m} to {res.date.max():%Y-%m}); run time {seconds:.0f} s.", "",
         "## Admissible exponents and setting (i)", "",
         f"- setting (i) refuted by the quotes: lower tail {res.refuted_lower.mean():.1%}, upper tail "
         f"{res.refuted_upper.mean():.1%}, either {(res.refuted_lower | res.refuted_upper).mean():.1%}",
         f"- γ̄ quartiles: {q(res.gamma_bar)}; η̄ quartiles: {q(res.eta_bar)}",
         f"- γ_i/γ̄ quartiles: {q(res.gamma_i / res.gamma_bar)}; (η_i − 1)/(η̄ − 1) quartiles: "
         f"{q((res.eta_i - 1) / (res.eta_bar - 1))}"]
    for side in ("lower", "upper"):
        ok = j[f"r4_smile_tail_ok_{side}_i"].astype("boolean")
        ref = j[f"refuted_{side}"]
        L.append(f"- {side} tail, refuted by the quotes vs smile diagnostic of estimate_moments (smile elasticity "
                 f"stays above the boundary value within 2 ATM s.d.): refuted and diagnostic fails {(ref & ~ok).sum()}, "
                 f"refuted and diagnostic passes {(ref & ok).sum()}, not refuted and diagnostic fails {(~ref & ~ok).sum()}, "
                 f"neither {(~ref & ok).sum()}")
    L += ["", "## Sign of the skewness", "",
          "| Setting | Admissible | Sign identified (sharp) | Negative | Sign identified (R4 box, estimate_moments) |",
          "| --- | --- | --- | --- | --- |"]
    for name, label in (("ii", "(ii) γ = η = 2"), ("i", "(i) boundary elasticities"),
                        *((f"t{int(round(100 * t)):02d}", f"θ = {t}") for t in THETAS)):
        adm = res[f"status_{name}"] == "ok"
        ident = res.loc[adm, f"sign_identified_{name}"].astype(bool)
        neg = (res.loc[adm, f"k3_hi_{name}"] < 0)
        box = ""
        if name in ("i", "ii"):
            ok = j["r4_status_i"] == "ok" if name == "i" else j["r4_skew_lo_ii"].notna()
            ex = ~((j.loc[ok, f"r4_skew_lo_{name}"] <= 0) & (0 <= j.loc[ok, f"r4_skew_hi_{name}"]))
            box = f"{ex.mean():.1%} of {ok.sum()}"
        L.append(f"| {label} | {adm.mean():.1%} | {ident.mean():.1%} | {neg.mean():.1%} | {box} |")
    th = res.theta_star
    L += ["", f"- breakdown θ* (sign identified for θ ≥ θ*; floor 0.01, tolerance 0.01): quartiles {q(th)}; "
          f"share ≤ 0.25: {(th <= 0.25).mean():.1%}, ≤ 0.5: {(th <= 0.5).mean():.1%}, ≤ 0.75: {(th <= 0.75).mean():.1%}, "
          f"≤ 0.9: {(th <= 0.9).mean():.1%}",
          f"- sign at θ = 1 (power tails): negative {(res.sign_star < 0).mean():.1%}; skewness quartiles {q(res.skew_pareto)}",
          "", "## By currency: share of months with the sign identified at θ = 0.75, and median θ*", ""]
    for c, g in res.groupby("currency"):
        adm = g["status_t75"] == "ok"
        L.append(f"- {c}: {g.loc[adm, 'sign_identified_t75'].astype(bool).mean():.1%} of {adm.sum()}, "
                 f"median θ* {g.theta_star.median():.2f}")
    if check is not None and len(check):
        m = res.merge(check, on=["date", "currency"], suffixes=("", "_fine"))
        dev = []
        for name in ("t50", "t75", "t90"):
            ok = (m[f"status_{name}"] == "ok") & (m[f"status_{name}_fine"] == "ok")
            w = (m.loc[ok, f"k3_hi_{name}_fine"] - m.loc[ok, f"k3_lo_{name}_fine"]).abs()
            for e in ("lo", "hi"):
                dev.append(((m.loc[ok, f"k3_{e}_{name}"] - m.loc[ok, f"k3_{e}_{name}_fine"]).abs() / w).max())
        L += ["", f"## Grid check ({len(m)} currency-months with four times as many step points)", "",
              f"- largest change of a third-central-moment bound, relative to the interval width: {np.nanmax(dev):.2e}",
              f"- largest change of θ*: {(m.theta_star - m.theta_star_fine).abs().max():.3f}"]
    return L


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--workers", type=int, default=os.cpu_count())
    p.add_argument("--limit", type=int, default=0, help="first n currency-months only (testing)")
    args = p.parse_args()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    mom = load("moments_1m.csv")
    mom = mom[mom.status == "ok"].sort_values(["date", "currency"])
    if args.limit:
        mom = mom.head(args.limit)
    panel = load("smile_panel.csv")
    fit = panel[panel.reading == "market"].set_index(["currency", "date"])
    cols = ["date", "currency", "F", "tau", "usd_base", "K_min", "K_max"]
    tasks = [(r, tuple(fit.loc[(r["currency"], r["date"]), ["alpha", "rho", "nu"]]), None)
             for r in mom[cols].to_dict("records")]
    t0 = time.time()
    with ProcessPoolExecutor(args.workers) as ex:
        res = pd.DataFrame(list(ex.map(one, [(r, s, 1000) for r, s, _ in tasks], chunksize=4)))
        sub = tasks[:: max(1, len(tasks) // GRID_CHECK)][:GRID_CHECK]
        check = pd.DataFrame(list(ex.map(one, [(r, s, 4000) for r, s, _ in sub], chunksize=1)))
    seconds = time.time() - t0
    res.to_csv(d / "identification_1m.csv", index=False)
    keep = ["date", "currency", "theta_star"] + [c for c in check.columns if c.startswith(("k3_", "status_"))]
    lines = summary(res, mom, check[keep].rename(columns={c: f"{c}_fine" for c in keep[2:]}), seconds)
    (d / "identification_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
