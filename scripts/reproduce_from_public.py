"""Recompute the paper's portfolio-level results from the public month series.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/reproduce_from_public.py \\
        [--public-dir data/public/fx_carry_portfolio_series] \\
        [--cboe-dir data/private/cboe/<retrieval date>/raw] \\
        [--verdelhan-file data/private/verdelhan/<retrieval date>/CurrencyPortfolios.xls] \\
        [--private-dir data/private/results/2026-09-23] \\
        [--out data/private/results/public_reproduction.md] [--workers 8]

Inputs. The three CSV files of data/public/fx_carry_portfolio_series/, its
manifest.json, whose SHA-256 hashes are checked before anything else, and the
mean theta0 of the 10-delta hedge stated in its README.md; and, each used only
if present, the Cboe VX settlement files downloaded by
scripts/acquire_cboe_vx.py (E5) and Verdelhan's CurrencyPortfolios.xls
downloaded by scripts/acquire_verdelhan.py (the comparison with published
portfolios, which needs the xlrd package). Without --cboe-dir or
--verdelhan-file the script looks where those scripts write
(data/private/cboe/<date>/raw, data/private/verdelhan/<date>/) and takes the
latest retrieval. It reads no LSEG data.

Computation. Every statistic is recomputed with the project's own functions,
imported unchanged: qef.stats.hac and qef.stats.bootstrap; regime_summary and
supplementary of estimate_e1.py; stambaugh_bootstrap, clark_west,
theta_confidence_set and e3_summary of estimate_stage4.py; e2_secondary of
estimate_moments.py; ar1_se, dummy_diff, joint_bootstrap, accepted_set and
GRID of estimate_regime_attribution.py; delta_method,
rejected_range, inversion_accepts, studentised_bootstrap, skewness and
t_fixed_lag of estimate_design_checks.py; window_stats and WINDOWS of
estimate_long_atm.py (whose theta_set applies the limit rule for
boundedness); vix_rolldown and ratio_set of estimate_e5.py with
qef.data.vix and qef.data.smile_inputs.option_dates; summarise of
robustness.py; and arch_difference of check_bootstrap_arch.py when the arch
package is installed. Where a script computes a statistic inside its main(),
which reads private files, I repeat those few lines here on the public frames;
each such place says so. No estimation script was changed. Bootstrap draws and
seeds are those behind the paper: 9,999 draws for E1, the E2 bias correction,
the realised-return split, the regime attribution and the bootstrap-t, and
1,999 for the secondary moment predictors and the robustness grid.

Comparisons. (1) Each recomputed number is compared with the number printed in
paper/main.tex, hard-coded in PAPER below with the paper's precision and section;
it passes if the two differ by at most half a unit in the last printed digit.
A few statements are bounds ("less than 4 per cent") or exact (counts, the kind
of a confidence set). (2) If my private results exist (--private-dir, default
data/private/results/2026-09-23), each recomputed number is also compared with
the value in my private summaries: the CSV summaries at full precision and the
markdown summaries at their printed precision. The tolerance is half a unit in
the summary's last printed digit plus PRIVATE_RTOL times the value, which
allows for the rounding of the public series to six significant figures
(PRIVATE_RTOL_EXTENDED for values that use extended.csv, rounded to four); for a
regime difference, and the statistics built from one, the relative part is
taken of the larger regime mean instead (private_scale), and for a bootstrap p
value a few draws (PRIVATE_P_DRAWS/B) are added.

E5 (paper, Section 5.4) is computed as estimate_e5.main computes it, for both
samples: the months are those with status ok of extended.csv and primary.csv,
in date order; the FX-volatility innovation is the residual of one AR(1)
fitted to their sigma_fx; and the VIX roll-down comes from the Cboe files.
Every E5 number is gated like any other.

Results that the public files cannot reproduce are listed in NOT_REPRODUCIBLE,
with the reason, and printed in the report. Among them are every result that
uses the oriented 10-delta risk reversal rr_or and the mean payoff and premium
of the long sample's at-the-money hedge, which are not published.

Exit status: 0 if every comparison passes, 1 if any fails, 2 if an input is
missing or a hash in the manifest does not match.
"""

from __future__ import annotations

import os

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import platform  # noqa: E402
import re  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from dataclasses import dataclass  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import estimate_long_atm as long_atm  # noqa: E402
from estimate_design_checks import (FIXED_LAGS, delta_method, inversion_accepts, rejected_range,  # noqa: E402
                                    skewness, studentised_bootstrap, t_fixed_lag)
from estimate_e1 import REGIME_BREAK, regime_summary, supplementary  # noqa: E402
from estimate_e5 import ratio_set, vix_rolldown  # noqa: E402
from estimate_moments import e2_secondary  # noqa: E402
from estimate_regime_attribution import GRID, accepted_set, ar1_se, dummy_diff, joint_bootstrap  # noqa: E402
from estimate_stage4 import clark_west, e3_summary, stambaugh_bootstrap, theta_confidence_set  # noqa: E402
from robustness import summarise  # noqa: E402

from qef.data.smile_inputs import option_dates  # noqa: E402
from qef.data.vix import load_vx  # noqa: E402
from qef.stats.bootstrap import bootstrap_distribution, optimal_block_length  # noqa: E402
from qef.stats.hac import mean_and_se, ols_hac  # noqa: E402

PUBLIC = ROOT / "data" / "public" / "fx_carry_portfolio_series"
PRIVATE = ROOT / "data" / "private" / "results" / "2026-09-23"
B_MAIN = 9999  # E1, E2, realised-return split, regime attribution, bootstrap-t (paper, Sections 4.4 and 5)
B_GRID = 1999  # secondary moment predictors and robustness grid (moments_summary.md; paper, Table 6)
COMBOS = (("market", 0.10),)  # the only E1 combination published; the 25-delta and smile-reading series are withheld
HEDGES = ("H10", "H25", "Hatm")
LA_WINDOWS = (*long_atm.WINDOWS, "n = 2 months", "n = 3 months")
ROB_VARIANTS = ("base", "hedge 25d", "hedge ATM")
PRIVATE_RTOL = 1e-4
PRIVATE_RTOL_EXTENDED = 2e-2  # values that use extended.csv (four significant figures); ratios of small means amplify the rounding
PRIVATE_P_DRAWS = 3


def dlab(delta: float) -> str:
    return str(round(100 * delta))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------------------------------------------------------------------
# Public inputs


def load_public(pub: Path) -> dict:
    man_path = pub / "manifest.json"
    if not man_path.exists():
        raise FileNotFoundError(f"{man_path} not found")
    man = json.loads(man_path.read_text())
    for name, meta in man["files"].items():
        if sha256(pub / name) != meta["sha256"]:
            raise ValueError(f"{name}: SHA-256 differs from manifest.json")
    m = re.search(r"theta0_mean_10d = ([0-9.eE+-]+)", (pub / "README.md").read_text())
    if m is None:
        raise ValueError("README.md states no theta0_mean_10d")
    read = lambda n: pd.read_csv(pub / n, parse_dates=["date"])
    return {"primary": read("primary.csv"), "extended": read("extended.csv"), "long_atm": read("long_atm.csv"),
            "theta0": float(m.group(1)), "manifest": man}


def e1_frame(primary: pd.DataFrame, reading: str, delta: float) -> pd.DataFrame:
    """One E1 series as estimate_e1.e1_series returns it: date, status, C_skew, FD, phi."""
    name = f"{reading}_{dlab(delta)}d"
    s = primary[["date", f"C_skew_{name}", "FD", f"phi_{name}"]].rename(
        columns={f"C_skew_{name}": "C_skew", f"phi_{name}": "phi"})
    return s.assign(status=np.where(s.C_skew.notna(), "ok", "missing")).sort_values("date").reset_index(drop=True)


def stage4_frame(data: dict) -> pd.DataFrame:
    """The columns of stage4_months.csv that the estimation code reads, for both samples.

    theta0 is the published mean in every primary month (only its mean enters any
    statistic); spot_sub is not published and only feeds a count that is not compared.
    rr_or is not published, so nothing below uses it. sigma_fx is carried for E5.
    """
    p, e = data["primary"], data["extended"]
    ret = ["U", "H10", "H25", "Hatm", "c1", "c2", "c3", "sigma_fx"]
    prim = p[["date", "status", *ret]].assign(
        phi=p.phi_market_10d, C_skew=p.C_skew_market_10d, FD=p.FD, sample="primary", theta0=data["theta0"])
    ext = e.loc[e.status == "ok", ["date", "status", "phi", "C_skew", "FD", *ret]].assign(
        sample="extended", theta0=np.nan)
    out = pd.concat([ext, prim], ignore_index=True).sort_values("date").reset_index(drop=True)
    return out.assign(spot_sub=np.nan)


def primary_ok(s4: pd.DataFrame) -> pd.DataFrame:
    return s4[(s4["sample"] == "primary") & (s4.status == "ok")].sort_values("date").reset_index(drop=True)


def prefixed(prefix: str, d: dict, skip=("variable",)) -> dict:
    return {f"{prefix}.{k}": v for k, v in d.items() if k not in skip}


# ---------------------------------------------------------------------------
# Tasks: each returns a dict of results keyed by a dotted name


def task_e1_combo(primary, reading, delta, B):
    s = e1_frame(primary, reading, delta)
    out = {}
    for col in ("phi", "C_skew", "FD"):
        out.update(prefixed(f"e1.{reading}.{dlab(delta)}.{col}", regime_summary(s, col, B)))
    return out


def task_e1_supplementary(primary, extended, B):
    s = e1_frame(primary, "market", 0.10)
    supp = supplementary(s, B)
    out = prefixed("e1.supp", supp)
    out["e1.supp.n_zero_rate_2020_2021"] = int((s.date < REGIME_BREAK).sum()) - supp["n_zero_rate_excluding_2020_2021"]
    ok = extended[extended.status == "ok"]
    out.update({"e1.ext.n_ok": len(ok), "e1.ext.mean_C_skew": ok.C_skew.mean(), "e1.ext.mean_FD": ok.FD.mean(),
                "e1.ext.mean_phi": ok.phi.mean(), "e1.ext.first": f"{ok.date.min():%Y-%m}",
                "e1.ext.last": f"{ok.date.max():%Y-%m}"})
    return out


def task_realised_split(s4, B):
    pm = primary_ok(s4)
    out = {}
    for col in ("U", "H10"):
        out.update(prefixed(f"split.{col}", regime_summary(pm, col, B)))
    return out


def task_regime_attribution(primary, B):
    """The statistics of estimate_regime_attribution.main, repeated on the public frames.

    Its third variable, the oriented risk reversal rr_or_vol, is not published (NOT_REPRODUCIBLE).
    """
    e1 = e1_frame(primary, "market", 0.10)
    e1 = e1[e1.status == "ok"].assign(C_skew_bp=lambda x: 1e4 * x.C_skew, FD_bp=lambda x: 1e4 * x.FD)
    out = {}
    for frame, col in ((e1, "C_skew_bp"), (e1, "FD_bp")):
        s = frame.assign(status="ok")
        r = regime_summary(s, col, B)
        a, b = s.loc[s.date < REGIME_BREAK, col].to_numpy(), s.loc[s.date >= REGIME_BREAK, col].to_numpy()
        (se_a, rho_a), (se_b, rho_b) = ar1_se(a), ar1_se(b)
        r.update(se_difference_ar1=float(np.hypot(se_a, se_b)), ar1_zero_rate=rho_a, ar1_hiking=rho_b)
        r["t_ar1"] = r["difference"] / r["se_difference_ar1"]
        out.update(prefixed(f"ra.{col}", r))
    hiking = (e1.date >= REGIME_BREAK).to_numpy()
    C, FD = e1.C_skew_bp.to_numpy(), e1.FD_bp.to_numpy()
    dC, dF = C[hiking].mean() - C[~hiking].mean(), FD[hiking].mean() - FD[~hiking].mean()
    nw = np.array([abs(dummy_diff(C - g * FD, hiking)[0]) <= 1.96 * dummy_diff(C - g * FD, hiking)[1] for g in GRID])
    draws, blocks = joint_bootstrap(C[~hiking], FD[~hiking], C[hiking], FD[hiking], B)
    z = draws[:, [0]] - GRID[None, :] * draws[:, [1]]
    lo, hi = np.percentile(z, [2.5, 97.5], axis=0)
    boot = (lo <= 0) & (0 <= hi)
    one = dummy_diff(C - FD, hiking)
    z1 = draws[:, 0] - draws[:, 1]
    kind_nw, lo_nw, hi_nw = accepted_set(nw)
    kind_b, lo_b, hi_b = accepted_set(boot)
    out.update({"ra.ratio": dC / dF, "ra.set_nw": kind_nw, "ra.set_nw_lo": lo_nw, "ra.set_nw_hi": hi_nw,
                "ra.set_boot": kind_b, "ra.set_boot_lo": lo_b, "ra.set_boot_hi": hi_b,
                "ra.boot_block_zero_rate": blocks[0], "ra.boot_block_hiking": blocks[1],
                "ra.one_for_one_diff_bp": one[0], "ra.one_for_one_se_nw": one[1], "ra.one_for_one_t_nw": one[0] / one[1],
                "ra.one_for_one_p_boot": float(2 * min((z1 >= 0).mean(), (z1 <= 0).mean())),
                "ra.corr_C_FD_months": float(np.corrcoef(C, FD)[0, 1])})
    return out


def task_e2(s4, name, B):
    """E2 of estimate_stage4.main for one predictor: phi (rr_or is not published, see NOT_REPRODUCIBLE)."""
    pm = primary_ok(s4)
    allm = s4[s4.status == "ok"].sort_values("date")
    reg = ols_hac(pm.U.to_numpy(), np.column_stack([np.ones(len(pm)), pm[name]]))
    r = {"b": reg["beta"][1], "se_hac": reg["se"][1], "t_hac": reg["beta"][1] / reg["se"][1], "hac_lag": reg["lag"],
         **stambaugh_bootstrap(pm.U.to_numpy(), pm[name].to_numpy(), B)}
    if name == "phi":
        r.update(clark_west(allm.date.to_numpy(), allm.U.to_numpy(), allm[name].to_numpy()))
    out = prefixed(f"e2.{name}", r)
    out["e2.n"] = len(pm)
    return out


def task_moments(primary, s4, B):
    """estimate_moments.e2_secondary on the published setting (i) predictors, merged as in its main()."""
    p = primary[primary.var_lo.notna()]
    pm = pd.DataFrame({"date": p.date, "var_lo_i": p.var_lo, "var_hi_i": p.var_hi,
                       "oskew_lo_i": p.oskew_lo, "oskew_hi_i": p.oskew_hi})
    for x in ("var", "oskew"):
        pm[f"{x}_mid_i"] = 0.5 * (pm[f"{x}_lo_i"] + pm[f"{x}_hi_i"])
    pm = pm.merge(primary_ok(s4)[["date", "U"]], on="date", how="inner").sort_values("date")
    out = {"mom.n": len(pm)}
    for _, r in e2_secondary(pm, B).iterrows():
        out.update(prefixed(f"mom.{r.predictor}.{r.endpoint}", r.to_dict(), skip=("predictor", "endpoint")))
    return out


def task_e3(s4):
    pm = primary_ok(s4)
    ext = s4[(s4["sample"] == "extended") & (s4.status == "ok")].sort_values("date")
    out = {}
    for sample, df in (("primary", pm), ("extended", ext)):
        r = e3_summary(df)
        r.pop("spot_substitutions")  # not published
        if sample == "extended":
            r.pop("theta0_reference")  # theta0 of the extended sample is not published
        out.update(prefixed(f"e3.{sample}", r))
    out["e3.primary.t_U"] = out["e3.primary.mean_U_bp"] / out["e3.primary.se_U_bp"]
    return out


def task_design(s4, hedge, B):
    """Section 3 of estimate_design_checks.delta_check for one hedge, from its building blocks."""
    pm = primary_ok(s4)
    theta0 = pm.theta0.mean()
    U, H = pm.U.to_numpy(), pm[hedge].to_numpy()
    r = delta_method(H, U)
    out = {"theta": r["theta"], "se_delta": r["se_delta"], "lag": r["lag"], "lo": r["lo"], "hi": r["hi"],
           "se_explicit_2x2": r["se_explicit_2x2"], "delta_contains_0": bool(r["lo"] <= 0 <= r["hi"])}
    kind, _, _ = theta_confidence_set(H, U)
    out["set_kind"] = kind
    if kind != "interval":
        a, b, _ = rejected_range(H, U)
        out.update(rej_lo=a, rej_hi=b)
    acc0, t0 = inversion_accepts(H, U, 0.0)
    cost = H - U
    c = mean_and_se(cost)
    bt = studentised_bootstrap(cost, B)
    drop = mean_and_se(np.delete(cost, np.argmax(cost)))
    out.update(acc0=acc0, t0=t0, mean_cost=c["mean"], se_cost=c["se"], lag_cost=c["lag"],
               t_delta=-np.sign(r["mean_U"]) * r["theta"] / r["se_delta"],
               se_ratio=c["se"] / (r["se_delta"] * abs(r["mean_U"])), boot_block=bt["block"], boot_q_lo=bt["q_lo"],
               boot_q_hi=bt["q_hi"], boot_lo=bt["lo"], boot_hi=bt["hi"], boot_share_below=bt["share_below"],
               skew=skewness(cost), t_drop=drop["mean"] / drop["se"])
    out.update({f"tfix0.L{k}": t_fixed_lag(cost, k) for k in FIXED_LAGS})
    if hedge == "H10":
        acc, t_th0 = inversion_accepts(H, U, theta0)
        x0 = H - (1 - theta0) * U
        ts = [t_fixed_lag(x0, k) for k in FIXED_LAGS]
        out.update(theta0=theta0, delta_contains_theta0=bool(r["lo"] <= theta0 <= r["hi"]), acc_th0=acc, t_th0=t_th0,
                   lag_th0=mean_and_se(x0)["lag"], lags_rej_th0=[k for k, t in zip(FIXED_LAGS, ts) if abs(t) > 1.96])
        out.update({f"tfixth0.L{k}": t for k, t in zip(FIXED_LAGS, ts)})
        out["e3_gap"] = float(np.max(np.abs(cost - (pm.c1 + pm.c2 + pm.c3).to_numpy())))
    res = prefixed(f"dc.{hedge}", out)
    if hedge == "H10":
        ext = s4[(s4["sample"] == "extended") & (s4.status == "ok")]
        re_ = mean_and_se((ext.H10 - ext.U).to_numpy())
        res["dc.ext_t"] = re_["mean"] / re_["se"]
    return res


def task_long_atm(la, window):
    """estimate_long_atm.window_stats for one window.

    window_stats also summarises the hedge's payoff and premium, which are not
    published; it is given Hatm - U in their place, and those entries are dropped.
    """
    ok = la[la.status == "ok"]
    ok = ok.assign(theta0=np.nan, theta0_pa=np.nan, payoff=ok.Hatm - ok.U, premium=ok.Hatm - ok.U)
    if window in long_atm.WINDOWS:
        df = ok[long_atm.WINDOWS[window](ok.date)]
    else:
        df = ok[ok.n_legs == int(window.split()[2])]
    st = long_atm.window_stats(df)
    for k in ("theta0", "theta0_pa", "CS pieces"):  # theta0 is not published; the pieces are compared through lo, hi
        st.pop(k)
    st = {k: v for k, v in st.items() if not k.startswith(("payoff ", "premium "))}
    return prefixed(f"la.{window}", st)


def task_robustness(primary, variant, B):
    """robustness.summarise on a variant frame built from the public series, as run_variant builds it."""
    col = {"base": "market_10d"}.get(variant)  # the 25-delta skew price is withheld: that row is checked through U and H
    hedge = {"base": "H10", "hedge 25d": "H25", "hedge ATM": "Hatm"}[variant]
    p = primary.sort_values("date")
    C = p[f"C_skew_{col}"] if col else 0.0 * p.FD
    df = pd.DataFrame({"date": p.date, "status": "ok", "C": C, "FD": p.FD, "U": p.U, "H": p[hedge], "c3": -C,
                       "arb_fail_legs": 0, "phi": p[f"phi_{col}"] if col else np.nan})
    r = summarise(df, B)
    if variant == "hedge 25d":  # its skew price is withheld: only the returns and theta_UB are recomputed
        r = {k: v for k, v in r.items() if k.startswith(("U_", "H_", "theta_UB", "theta_")) or k in ("n", "n_months")}
    return prefixed(f"rob.{variant}", r)


def task_arch(primary, reading, delta):
    """The comparison of check_bootstrap_arch.main for one E1 combination."""
    from check_bootstrap_arch import B as ARCH_B
    from check_bootstrap_arch import arch_difference
    s = e1_frame(primary, reading, delta)
    out = {}
    for col in ("phi", "C_skew", "FD"):
        g = s[s.status == "ok"].dropna(subset=[col])
        a, b = g.loc[g.date < REGIME_BREAK, col].to_numpy(), g.loc[g.date >= REGIME_BREAK, col].to_numpy()
        mine = np.percentile(bootstrap_distribution(lambda x, y: y.mean() - x.mean(), [a, b], B=ARCH_B), [2.5, 97.5])
        theirs = arch_difference(a, b)
        scale = 1e4 if col != "phi" else 1.0
        key = f"arch.{reading}.{dlab(delta)}.{col}"
        out.update({f"{key}.lo_mine": scale * mine[0], f"{key}.hi_mine": scale * mine[1],
                    f"{key}.lo_arch": scale * theirs[0], f"{key}.hi_arch": scale * theirs[1],
                    f"{key}.b_mine_zero_rate": optimal_block_length(a), f"{key}.b_mine_hiking": optimal_block_length(b)})
    return out


def fx_vol_innovation(s4: pd.DataFrame) -> tuple[pd.DataFrame, float, np.ndarray]:
    """The E5 months and FX-volatility innovation of estimate_e5.main.

    The months are those with a realised return (status ok) of the extended and
    primary samples, in date order; the innovation is the residual of one AR(1)
    fitted by least squares to sigma_fx over all of them (NaN in the first month).
    Returns the months, the AR(1) coefficient and the innovation.
    """
    months = s4[s4.status == "ok"].sort_values("date").reset_index(drop=True)
    x = months.sigma_fx.to_numpy()
    c, rho = np.linalg.lstsq(np.column_stack([np.ones(len(x) - 1), x[:-1]]), x[1:], rcond=None)[0]
    return months, float(rho), np.r_[np.nan, x[1:] - c - rho * x[:-1]]


def task_e5(s4, cboe_dir):
    """E5 of estimate_e5.main, both samples, repeated on the public frames.

    The factors are the VIX roll-down from the Cboe files and the FX-volatility
    innovation of fx_vol_innovation; as there, the regressions run on each
    sample's months with both factors.
    """
    contracts, rescaled = load_vx(Path(cboe_dir))
    months, rho, innovation = fx_vol_innovation(s4)
    rows = []
    for _, m in months.iterrows():
        end = option_dates(m.date, "EUR")[2]
        r_vix, _ = vix_rolldown(contracts, m.date, end)
        rows.append({"date": m.date, "r_vix": r_vix})
    f = pd.DataFrame(rows).assign(d_sigma_fx=innovation)
    df = months.merge(f, on="date")
    out = {}
    for sample in ("primary", "extended"):  # the loop of estimate_e5.main
        s = df[df["sample"] == sample].dropna(subset=["r_vix", "d_sigma_fx"])
        Fm = s[["r_vix", "d_sigma_fx"]].to_numpy()
        X = np.column_stack([np.ones(len(s)), Fm])
        res = {"n": len(s)}
        for name in ("U", "H10"):
            r = ols_hac(s[name].to_numpy(), X)
            for k, lab in enumerate(("alpha_bp", "beta_vix", "beta_fxvol")):
                scale = 1e4 if k == 0 else 1.0
                res[f"{name}_{lab}"] = scale * r["beta"][k]
                res[f"{name}_{lab}_t"] = r["beta"][k] / r["se"][k]
            fit = X @ r["beta"]
            res[f"{name}_R2"] = 1 - np.var(s[name] - fit) / np.var(s[name])
        res["theta_UB_raw"] = 1 - s.H10.mean() / s.U.mean()
        res["theta_UB_adjusted"] = 1 - res["H10_alpha_bp"] / res["U_alpha_bp"]
        kind, lo, hi = ratio_set(s.H10.to_numpy(), s.U.to_numpy(), Fm)
        res.update(theta_adj_set=kind, theta_adj_lo=lo, theta_adj_hi=hi)
        res["corr_U_rvix"] = float(np.corrcoef(s.U, s.r_vix)[0, 1])
        res["corr_H10_rvix"] = float(np.corrcoef(s.H10, s.r_vix)[0, 1])
        out.update(prefixed(f"e5.{sample}", res))
        # for the paper's statement that neither portfolio earns a significant alpha
        out[f"e5.{sample}.max_abs_alpha_t"] = max(abs(res["U_alpha_bp_t"]), abs(res["H10_alpha_bp_t"]))
    prim = df[(df["sample"] == "primary")].dropna(subset=["r_vix"])  # the factor-free statistics, primary months
    out["e5.n_months"] = len(months)
    out["e5.n_r_vix"] = len(prim)
    out["e5.corr_U_rvix"] = float(np.corrcoef(prim.U, prim.r_vix)[0, 1])
    out["e5.corr_H10_rvix"] = float(np.corrcoef(prim.H10, prim.r_vix)[0, 1])
    out["e5.ar1_rho"] = rho
    out["e5.vix_loading_reduction"] = 1 - out["e5.primary.H10_beta_vix"] / out["e5.primary.U_beta_vix"]
    out["e5.rescaling_first_date"] = min(d for _, d in rescaled) if rescaled else ""
    return out


def task_verdelhan(primary, xls):
    """The statistics of compare_verdelhan.main, repeated with the public HML^U."""
    lrv = pd.read_excel(xls, sheet_name="Developed currencies")
    lrv = lrv[pd.to_datetime(lrv["Dates"], errors="coerce").notna()][["Dates", "HML = P5 - P1"]]
    lrv = lrv.rename(columns={"Dates": "lrv_date", "HML = P5 - P1": "lrv_hml"}).astype({"lrv_hml": float})
    lrv["key"] = pd.to_datetime(lrv.lrv_date).dt.to_period("M")
    s4 = primary.loc[primary.status == "ok", ["date", "U"]].copy()
    s4["key"] = s4.date.dt.to_period("M") + 1  # the return from month-end t is paired with their next month-end
    m = s4.merge(lrv, on="key", how="inner").sort_values("date")
    u, h = m.U.to_numpy(), m.lrv_hml.to_numpy()
    mu, mh, md = mean_and_se(u), mean_and_se(h), mean_and_se(u - h)
    slope = ols_hac(u, np.column_stack([np.ones(len(h)), h]))
    return {"vd.n": len(m), "vd.first": f"{pd.Timestamp(m.lrv_date.min()):%Y-%m}",
            "vd.last": f"{pd.Timestamp(m.lrv_date.max()):%Y-%m}",
            "vd.mean_u_bp": 1e4 * mu["mean"], "vd.se_u_bp": 1e4 * mu["se"], "vd.mean_h_bp": 1e4 * mh["mean"],
            "vd.se_h_bp": 1e4 * mh["se"], "vd.mean_d_bp": 1e4 * md["mean"], "vd.se_d_bp": 1e4 * md["se"],
            "vd.corr": float(np.corrcoef(u, h)[0, 1]), "vd.same_sign": float((np.sign(u) == np.sign(h)).mean()),
            "vd.slope": slope["beta"][1], "vd.slope_se": slope["se"][1],
            "vd.sd_u_bp": 1e4 * u.std(ddof=1), "vd.sd_h_bp": 1e4 * h.std(ddof=1)}


def counts(data: dict) -> dict:
    p, e, la = data["primary"], data["extended"], data["long_atm"]
    ok = la[la.status == "ok"]
    return {"count.primary": len(p), "count.primary_zero_rate": int((p.date < REGIME_BREAK).sum()),
            "count.primary_hiking": int((p.date >= REGIME_BREAK).sum()), "count.primary_returns": int((p.status == "ok").sum()),
            "count.extended_ok": int((e.status == "ok").sum()), "count.long_atm": len(ok),
            "count.long_atm_n2": int((ok.n_legs == 2).sum()), "count.long_atm_first": f"{ok.date.min():%Y-%m}",
            "count.long_atm_last": f"{ok.date.max():%Y-%m}"}


def derived(R: dict) -> dict:
    out = {}
    base = [(R[f"rob.base.{a}"], b) for a, b in (("phi_diff", R["e1.market.10.phi.difference"]),
                                                 ("U_bp", R["e3.primary.mean_U_bp"]), ("H_bp", R["e3.primary.mean_H10_bp"]),
                                                 ("skew_term_bp", R["e3.primary.mean_c3_bp"]))]
    out["rob.base.reproduction_gap"] = max(abs(x - y) for x, y in base)
    if "arch.market.10.phi.lo_mine" in R:
        rows = [(r, d, c) for r, d in COMBOS for c in ("phi", "C_skew", "FD")]
        gap = lambda cols: max(abs(R[f"arch.{r}.{dlab(d)}.{c}.{e}_mine"] - R[f"arch.{r}.{dlab(d)}.{c}.{e}_arch"])
                               for r, d, c in rows if c in cols for e in ("lo", "hi"))
        out["arch.max_phi"] = gap(("phi",))
        out["arch.max_bp"] = gap(("C_skew", "FD"))
        excl = lambda r, d, c, w: R[f"arch.{r}.{dlab(d)}.{c}.lo_{w}"] > 0 or R[f"arch.{r}.{dlab(d)}.{c}.hi_{w}"] < 0
        out["arch.agree"] = all(excl(r, d, c, "mine") == excl(r, d, c, "arch") for r, d, c in rows)
    return out


# ---------------------------------------------------------------------------
# The paper's printed numbers (paper/main.tex, with the corrections of 1 October 2026)


@dataclass
class Expect:
    key: str
    printed: object  # the printed number as text, or an exact value
    where: str
    label: str
    scale: float = 1.0  # multiplies the recomputed value into the paper's units
    kind: str = "round"  # round, exact, le (rounded value at most the bound) or lt (strictly below)


def _e1_rows():
    rows = []
    fields = ("mean_zero_rate", "se_zero_rate", "mean_hiking", "se_hiking", "difference", "se_difference_hac",
              "boot_lo", "boot_hi")
    table = {("market", "10", "phi", 1.0): ("0.809", "0.185", "0.454", "0.066", "-0.355", "0.204", "-0.814", "0.014"),
             ("market", "10", "C_skew", 1e4): ("11.24", "0.82", "11.19", "0.83", "-0.05", "1.30", "-2.44", "2.20"),
             ("market", "10", "FD", 1e4): ("18.61", "1.77", "27.05", "1.41", "8.44", "2.30", "3.48", "13.65")}
    for (reading, d, col, scale), vals in table.items():
        for f, v in zip(fields, vals):
            rows.append(Expect(f"e1.{reading}.{d}.{col}.{f}", v, "Table 1 (Section 5.1)",
                               f"{col}, {reading} reading, {d}-delta: {f}", scale))
    rows += [Expect("e1.market.10.phi.n_zero_rate", 104, "Table 1 (Section 5.1)", "zero-rate month-ends", kind="exact"),
             Expect("e1.market.10.phi.n_hiking", 56, "Table 1 (Section 5.1)", "hiking month-ends", kind="exact")]
    return rows


def _supp_rows():
    w = "Section 5.1 (supplementary, post hoc)"
    vals = [("ar1_zero_rate", "0.85"), ("ar1_hiking", "0.76"), ("effective_n_zero_rate", "8"), ("effective_n_hiking", "8"),
            ("se_difference_ar1_plugin", "0.27"), ("t_difference_ar1_plugin", "-1.31"),
            ("difference_excluding_2020_2021", "-0.03"), ("median_zero_rate", "0.53"), ("median_hiking", "0.38"),
            ("trimmed10_zero_rate", "0.63"), ("trimmed10_hiking", "0.41"), ("mean_log_phi_zero_rate", "-0.46"),
            ("mean_log_phi_hiking", "-0.89"), ("ratio_of_means_zero_rate", "0.60"), ("ratio_of_means_hiking", "0.41"),
            ("se_difference_lag_0", "0.079"), ("se_difference_lag_9", "0.204"), ("se_difference_lag_20", "0.236"),
            ("bootstrap_p_two_sided_percentile", "0.064"), ("difference_drop_top1_zero_rate", "-0.32"),
            ("difference_drop_top3_zero_rate", "-0.28"), ("difference_drop_top5_zero_rate", "-0.23"),
            ("difference_drop_top10_zero_rate", "-0.14")]
    rows = [Expect(f"e1.supp.{k}", v, w, k.replace("_", " ")) for k, v in vals]
    rows.append(Expect("e1.supp.n_zero_rate_2020_2021", 24, w, "zero-rate month-ends in 2020-2021", kind="exact"))
    return rows


def _paper():
    S32, S51, S52, S53, S54, S46, S6 = ("Section 3.2", "Section 5.1", "Section 5.2", "Section 5.3", "Table 5 (Section 5.4)",
                                        "Section 4.6", "Table 6 (Section 6)")
    E = Expect
    rows = [E("count.primary", 160, S32, "primary month-ends", kind="exact"),
            E("count.primary_zero_rate", 104, S32, "zero-rate month-ends", kind="exact"),
            E("count.primary_hiking", 56, S32, "hiking month-ends", kind="exact"),
            E("count.primary_returns", 159, S32, "primary months with realised returns", kind="exact"),
            E("count.extended_ok", 72, S32, "extended month-ends with four currencies", kind="exact"),
            E("count.long_atm", 379, S32, "long ATM months with realised returns", kind="exact"),
            E("count.long_atm_n2", 36, S32, "long ATM months with two legs per side", kind="exact"),
            E("count.long_atm_first", "1995-01", S32, "first long ATM month-end", kind="exact"),
            *_e1_rows(), *_supp_rows(),
            E("ra.C_skew_bp.difference", "-0.05", S51, "regime change of C_skew (bp)"),
            E("ra.C_skew_bp.se_difference_hac", "1.30", S51, "its Newey-West s.e."),
            E("ra.C_skew_bp.se_difference_ar1", "1.29", S51, "its AR(1) plug-in s.e."),
            E("ra.FD_bp.t_ar1", "2.80", S51, "t of the FD change with the AR(1) plug-in s.e."),
            E("ra.ratio", "-0.006", S51, "Wald ratio dC_skew/dFD"),
            E("ra.set_nw", "interval", S51, "kind of the Newey-West set", kind="exact"),
            E("ra.set_nw_lo", "-0.26", S51, "Newey-West set, lower end"), E("ra.set_nw_hi", "0.49", S51, "Newey-West set, upper end"),
            E("ra.set_boot_lo", "-0.23", S51, "bootstrap set, lower end"), E("ra.set_boot_hi", "0.48", S51, "bootstrap set, upper end"),
            E("ra.one_for_one_t_nw", "-2.66", S51, "t of the one-for-one test"),
            E("ra.one_for_one_p_boot", "0.009", S51, "bootstrap p of the one-for-one test"),
            E("e1.ext.n_ok", 72, S51, "extended month-ends", kind="exact"),
            E("e1.ext.mean_C_skew", "22.3", S51, "extended sample mean C_skew (bp)", 1e4),
            E("e1.ext.mean_phi", "0.64", S51, "extended sample mean phi"),
            E("split.U.mean_zero_rate", "7.4", "Sections 1 and 5.3", "HML^U, zero-rate mean (bp)", 1e4),
            E("split.U.se_zero_rate", "12.4", S53, "its s.e. (bp)", 1e4),
            E("split.U.mean_hiking", "37.0", "Sections 1 and 5.3", "HML^U, hiking mean (bp)", 1e4),
            E("split.U.se_hiking", "17.8", S53, "its s.e. (bp)", 1e4),
            E("split.U.difference", "29.6", "Section 1", "HML^U regime difference (bp)", 1e4),
            E("split.U.se_difference_hac", "23.8", "Section 1", "its Newey-West s.e. (bp)", 1e4),
            E("split.U.boot_lo", "-19.5", S53, "bootstrap 95% lower bound (bp)", 1e4),
            E("split.U.boot_hi", "80.0", S53, "bootstrap 95% upper bound (bp)", 1e4),
            E("split.H10.mean_zero_rate", "-2.7", S53, "HML^H (10-delta), zero-rate mean (bp)", 1e4),
            E("split.H10.se_zero_rate", "15.7", S53, "its s.e. (bp)", 1e4),
            E("split.H10.mean_hiking", "29.3", S53, "HML^H (10-delta), hiking mean (bp)", 1e4),
            E("split.H10.se_hiking", "16.0", S53, "its s.e. (bp)", 1e4),
            E("split.H10.boot_lo", "-16.8", S53, "bootstrap 95% lower bound (bp)", 1e4),
            E("split.H10.boot_hi", "80.3", S53, "bootstrap 95% upper bound (bp)", 1e4)]
    T2 = "Table 2 (Section 5.2)"
    for name, vals in (("phi", ("0.0062", "5.95", "0.00055", "0.00049", "0.0057", "0.021")),):
        for f, v in zip(("b", "t_hac", "bias_bootstrap", "bias_analytic", "b_bias_corrected", "p_one_sided_b_gt_0"), vals):
            rows.append(E(f"e2.{name}.{f}", v, T2, f"{name}: {f}"))
    rows += [E("e2.n", 159, T2, "months in sample", kind="exact"),
             E("e2.phi.n_forecasts", 116, T2, "out-of-sample forecasts", kind="exact"),
             E("e2.phi.oos_r2", "3.7", T2, "phi: out-of-sample R2 (%)", 100),
             E("e2.phi.cw_t", "1.639", T2, "phi: Clark-West t"), E("e2.phi.cw_p_one_sided", "0.051", T2, "phi: its p"),
             E("e2.phi.cw_t_nw", "1.50", T2, "phi: Clark-West t, Newey-West"),
             E("e2.phi.cw_p_one_sided_nw", "0.067", T2, "phi: its p"),
             E("e2.phi.rho_predictor", "0.86", S52, "AR(1) of phi"),
             E("e2.phi.corr_u_v", "-0.39", S52, "phi: correlation of innovations with returns"),
             E("rob.base.CW_t", "0.96", S52, "Clark-West t with training from the primary start")]
    T3 = "Table 3 (Section 5.2)"
    for pred, b, t, p, pk in (("var", ("4.98", "4.81", "4.64"), ("0.87", "0.90", "0.93"), ("0.185", "0.171", "0.163"), "p_boot_b_gt_0"),
                              ("oskew", ("-0.024", "-0.028", "-0.028"), ("-2.74", "-2.59", "-1.87"), ("0.013", "0.015", "0.045"),
                               "p_boot_b_lt_0")):
        for i, end in enumerate(("lo", "mid", "hi")):
            rows += [E(f"mom.{pred}.{end}.b", b[i], T3, f"{pred} at {end}: b"),
                     E(f"mom.{pred}.{end}.t_hac", t[i], T3, f"{pred} at {end}: Newey-West t"),
                     E(f"mom.{pred}.{end}.{pk}", p[i], T3, f"{pred} at {end}: one-sided bootstrap p")]
    rows.append(E("mom.n", 159, T3, "months", kind="exact"))
    T4 = "Table 4 (Section 5.3)"
    for name, m, s in (("U", "17.7", "12.1"), ("H10", "8.3", "11.9"), ("H25", "9.1", "11.3"), ("Hatm", "9.7", "8.9"),
                       ("c1", "0.5", "3.8"), ("c2", "1.4", "1.4"), ("c3", "-11.2", "0.65")):
        rows += [E(f"e3.primary.mean_{name}_bp", m, T4, f"mean {name} (bp)"), E(f"e3.primary.se_{name}_bp", s, T4, f"s.e. {name} (bp)")]
    rows += [E("e3.primary.n_months", 159, T4, "months", kind="exact"),
             E("e3.primary.t_U", "1.47", S53, "t of mean HML^U"),
             E("dc.H10.mean_cost", "-9.3", S53, "mean H10 - U (bp)", 1e4), E("dc.H10.se_cost", "4.0", S53, "its s.e. (bp)", 1e4),
             E("dc.H10.t0", "-2.35", S53, "its t"), E("dc.H10.skew", "3.6", S53, "skewness of H10 - U"),
             E("dc.H10.boot_q_lo", "-2.78", S53, "bootstrap-t 2.5% quantile"),
             E("dc.H10.boot_lo", "-16.4", S53, "bootstrap-t 95% interval, lower (bp)", 1e4),
             E("dc.H10.boot_hi", "1.7", S53, "bootstrap-t 95% interval, upper (bp)", 1e4),
             E("dc.H10.boot_share_below", "0.04", S53, "one-sided bootstrap p"),
             E("dc.H10.t_drop", "-3.55", S53, "t without the month with the largest payoff"),
             E("dc.ext_t", "-0.25", S53, "t of mean H10 - U, extended sample"),
             E("e3.primary.theta_UB_H10", "0.53", S53, "theta_UB, 10-delta"),
             E("e3.primary.theta_UB_H25", "0.49", S53, "theta_UB, 25-delta"),
             E("e3.primary.theta_UB_Hatm", "0.45", S53, "theta_UB, ATM"),
             E("e3.primary.theta0_reference", "0.10", S53, "theta0, 10-delta"),
             E("e3.primary.theta_UB_H10_set", "unbounded", S53, "kind of the 10-delta set", kind="exact"),
             E("dc.H25.set_kind", "unbounded", S53, "kind of the 25-delta set", kind="exact"),
             E("dc.Hatm.set_kind", "unbounded", S53, "kind of the ATM set", kind="exact"),
             E("dc.H10.acc0", False, S53, "the 10-delta set excludes 0 (normal reference)", kind="exact"),
             E("dc.H10.rej_lo", "-1.23", S53, "10-delta set: line less [a, b], a"),
             E("dc.H10.rej_hi", "0.09", S53, "10-delta set: line less [a, b], b"),
             E("dc.H10.acc_th0", True, S53, "the 10-delta set retains theta0", kind="exact"),
             E("dc.H10.t_th0", "-1.93", S53, "t at theta0"), E("dc.H10.lag_th0", 6, S53, "its automatic lag", kind="exact"),
             E("dc.H10.lags_rej_th0", [0, 1, 2, 3, 4, 5], S53, "fixed lags at which theta0 is rejected", kind="exact"),
             E("dc.H10.lo", "-0.21", S53, "delta-method interval, 10-delta, lower"),
             E("dc.H10.hi", "1.27", S53, "delta-method interval, 10-delta, upper"),
             E("dc.H25.lo", "-0.30", S53, "delta-method interval, 25-delta, lower"),
             E("dc.H25.hi", "1.27", S53, "delta-method interval, 25-delta, upper"),
             E("dc.Hatm.lo", "-0.30", S53, "delta-method interval, ATM, lower"),
             E("dc.Hatm.hi", "1.20", S53, "delta-method interval, ATM, upper"),
             E("e3.extended.n_months", 72, S53, "extended months", kind="exact"),
             E("e3.extended.mean_c1_bp", "16.7", S53, "extended: term (i) (bp)"),
             E("e3.extended.se_c1_bp", "10.2", S53, "extended: its s.e. (bp)"),
             E("e3.extended.mean_c3_bp", "-22.3", S53, "extended: skew term (bp)"),
             E("e3.extended.se_c3_bp", "3.0", S53, "extended: its s.e. (bp)"),
             E("e3.extended.mean_U_bp", "8.4", S53, "extended: mean HML^U (bp)"),
             E("e3.extended.se_U_bp", "62.7", S53, "extended: its s.e. (bp)"),
             E("la.full.months", 379, S53, "long ATM months", kind="exact"),
             E("la.full.U mean", "37.2", S53, "long ATM: mean HML^U (bp)"), E("la.full.U se", "13.2", S53, "its s.e."),
             E("la.full.U t", "2.83", S53, "its t"), E("la.full.Hatm mean", "16.3", S53, "long ATM: mean Hatm (bp)"),
             E("la.full.Hatm se", "7.7", S53, "its s.e."), E("la.full.Hatm-U mean", "-20.9", S53, "long ATM: Hatm - U (bp)"),
             E("la.full.Hatm-U se", "7.9", S53, "its s.e."), E("la.full.theta_UB", "0.56", S53, "long ATM: theta_UB"),
             E("la.full.CS", "interval", S53, "long ATM: the set is bounded", kind="exact"),
             E("la.full.CS lo", "0.27", S53, "long ATM: set, lower end"), E("la.full.CS hi", "0.94", S53, "long ATM: set, upper end"),
             E("la.full without 2008.CS lo", "0.33", S53, "without 2008: set, lower end"),
             E("la.full without 2008.CS hi", "0.85", S53, "without 2008: set, upper end"),
             E("la.2008.Hatm mean", "-99", S53, "2008: mean Hatm (bp)"), E("la.2008.U mean", "-222", S53, "2008: mean HML^U (bp)")]
    for name, v in (("U_alpha_bp", "-7.3"), ("U_alpha_bp_t", "-0.61"), ("U_beta_vix", "0.049"), ("U_beta_vix_t", "7.41"),
                    ("U_beta_fxvol", "-3.11"), ("U_beta_fxvol_t", "-2.16"), ("U_R2", "0.33"), ("H10_alpha_bp", "-10.6"),
                    ("H10_alpha_bp_t", "-0.79"), ("H10_beta_vix", "0.039"), ("H10_beta_vix_t", "4.69"),
                    ("H10_beta_fxvol", "-1.76"), ("H10_beta_fxvol_t", "-1.21"), ("H10_R2", "0.22"),
                    ("theta_UB_adjusted", "-0.45")):
        rows.append(E(f"e5.primary.{name}", v, S54, name))
    rows += [E("e5.primary.n", 159, S54, "months", kind="exact"), E("e5.ar1_rho", "0.76", S54, "AR(1) of the FX-volatility level"),
             E("e5.n_r_vix", 159, S54, "months with the VIX roll-down", kind="exact"),
             E("e5.corr_U_rvix", "0.55", "Section 5.4", "correlation of HML^U with the VIX roll-down"),
             E("e5.primary.theta_adj_set", "unbounded", "Section 5.4", "kind of the exposure-adjusted set", kind="exact"),
             E("e5.vix_loading_reduction", "0.2", "Sections 1 and 5.4", "the hedge lowers the VIX loading by about a fifth"),
             E("e5.extended.n", 65, "Section 5.4", "extended sample: months", kind="exact"),
             E("e5.extended.U_beta_vix", "0.124", "Section 5.4", "extended sample: VIX loading of HML^U"),
             E("e5.extended.U_beta_vix_t", "5.39", "Section 5.4", "extended sample: its t"),
             E("e5.extended.H10_beta_vix", "0.112", "Section 5.4", "extended sample: VIX loading of HML^H"),
             E("e5.extended.H10_beta_vix_t", "4.53", "Section 5.4", "extended sample: its t"),
             E("e5.extended.max_abs_alpha_t", "1.96", "Section 5.4",
               "extended sample: neither alpha is significant (largest absolute t below 1.96)", kind="lt"),
             E("e5.rescaling_first_date", "2007-03-26", "Section 3.5", "Cboe archive rescaled before 26 March 2007", kind="exact"),
             E("vd.n", 96, S46, "overlapping months", kind="exact"), E("vd.first", "2013-06", S46, "first overlapping month", kind="exact"),
             E("vd.mean_u_bp", "7.8", S46, "mean HML^U (bp)"), E("vd.se_u_bp", "16.1", S46, "its s.e."),
             E("vd.mean_h_bp", "7.2", S46, "mean developed HML (bp)"), E("vd.se_h_bp", "26.5", S46, "its s.e."),
             E("vd.corr", "0.51", S46, "correlation"), E("vd.same_sign", "77", S46, "months with the same sign (%)", 100),
             E("vd.slope", "0.41", S46, "slope of HML^U on their HML"), E("vd.slope_se", "0.12", S46, "its s.e."),
             E("arch.max_phi", "0.0093", S46, "arch moves no phi bound by more than 0.0093", kind="le"),
             E("arch.max_bp", "0.031", S46, "arch moves no C_skew or FD bound by more than 0.031 bp", kind="le"),
             E("arch.agree", True, S46, "arch changes no conclusion", kind="exact"),
             E("rob.base.reproduction_gap", "1e-9", S46, "the robustness base reproduces E1 and Stage 4", kind="lt")]
    for v, vals in (("base", ("-0.355", "0.204", "-0.82", "0.01", "0.0058", "0.017", "0.96", "17.7", "1.47", "8.3", "0.70", "-11.2", "0.53")),
                    ("hedge 25d", (None,) * 7 + ("17.7", "1.47", "9.1", "0.81", None, "0.49")),
                    ("hedge ATM", (None,) * 7 + ("17.7", "1.47", "9.7", "1.09", "0", "0.45"))):
        for f, x in zip(("phi_diff", "phi_diff_se", "phi_diff_boot_lo", "phi_diff_boot_hi", "E2_b_bc", "E2_p", "CW_t",
                         "U_bp", "U_t", "H_bp", "H_t", "skew_term_bp", "theta_UB"), vals):
            if x is not None:
                rows.append(E(f"rob.{v}.{f}", x, S6, f"{v}: {f}"))
    return rows


PAPER = _paper()

RR_REASON = ("rr_or is not published: it is a plain signed average of six raw 10-delta risk-reversal quotes in quote "
             "units, the least transformed of the portfolio-level series")
NOT_REPRODUCIBLE = [
    ("Abstract; Section 5.2; Table ident; Section 7", "Identification of the sign of one-month risk-neutral skewness: "
     "admissible exponents (medians 59.7 and 67.8), setting (i) refuted in 75.3% of currency-months, sign identified in "
     "93.3%, 26.7% and 70.3%, breakdown quartiles 0.49, 0.63, 0.80, median skewness -0.34",
     "computed per currency-month from each calibrated smile (estimate_identification.py); no per-currency value is published"),
    ("Section 5.2", "Per-currency moment intervals: tail exponents of setting (i) (median 61 and 67), median skewness "
     "interval [-0.70, 0.14], 32% excluding zero, setting (ii) about fifty times wider, smile extrapolation contradicting "
     "setting (i) in about three quarters", "per currency-month (moments_1m.csv); only the portfolio averages are published"),
    ("Section 5.2; Table wing", "Out-of-sample test of the tail hypotheses with one-month 5-delta quotes (99.1%, 99.6%, "
     "theta_5 quartiles, 90.6% of 446, breakdown quartiles with and without the 5-delta prices)",
     "needs per-currency Fenics 5-delta quotes and smiles (estimate_wing_test.py)"),
    ("Table 2, Panels A and B; Section 5.2", "Every E2 result of the oriented 10-delta risk reversal: in sample b 0.397 "
     "(t 4.13), biases 0.0183 and 0.0148, corrected b 0.379 (p 0.004); out of sample R2 -5.5%, Clark-West -1.11 (p 0.87), "
     "Newey-West -1.18 (p 0.88); its AR(1) coefficient 0.52, its correlation 0.56 with phi and the correlation -0.53 of "
     "its innovations with returns", RR_REASON),
    ("Section 5.1", "Regime change of the oriented 10-delta risk reversal (0.21 volatility points, s.e. 0.29, bootstrap "
     "interval [-0.34, 0.67])", RR_REASON),
    ("Section 5.3", "Mean payoff (502 bp per month) and mean premium (378 bp) of the at-the-money hedge in 2008",
     "the long sample's payoff and premium are not published: in 42 of the 379 months exactly one leg's ATM option "
     "ends in the money, so the payoff would depend on that leg alone; Hatm - U, which always mixes all legs, is "
     "published through U and Hatm, and the 2008 means of both are reproduced"),
    ("Section 3.5; Section 5.4", "Cboe second-month settlements equal LSEG's second-month continuation at 223 of 223 "
     "month-ends", "needs LSEG's VX continuation series"),
    ("Section 5.3", "theta0 of the at-the-money hedge in the long sample (0.50 with the pips or premium-adjusted forward "
     "delta; every sub-period's set contains it)", "theta0 is published only as the mean of the 10-delta hedge"),
    ("Section 3.3", "Leg-by-leg splice of Fenics and composite skew costs (correlations 0.995 and 0.997), the second "
     "retrieval (243 instruments byte-identical) and the audit outcomes", "per-currency quotes"),
    ("Sections 3.4 and 4.1", "Calibration diagnostics: convergence, residuals below 1e-13, static-arbitrage checks, the "
     "R5(c) slope condition, Jacobian condition numbers and the held-out 10-delta fit of the two butterfly readings",
     "per-currency smiles and quotes"),
    ("Section 4.6", "SOFR OIS check (0.013%, 1.1e-4, 0.002 bp, E1 difference -0.3548 under both rates) and named-broker "
     "check (96 EUR month-ends, 0.02 points, 0.07 bp on 2.7 bp, 0.03 bp on 11.5 bp)",
     "recalibration of per-currency smiles from quotes (estimate_design_checks.py, checks 1 and 2)"),
    ("Section 4.6", "Closed-form checks of the moment code (1.6e-15 in all 1,440 currency-months) and the C++ kernel's "
     "parity and timing", "per currency-month and synthetic inputs respectively"),
    ("Section 4.6; Table models", "Model validation of the moment code (Merton and Heston models)",
     "synthetic and already public: scripts/validate_moments_models.py and tests/test_moments_models.py"),
    ("Table 1 and Sections 5.1 and 6", "E4 (the smile-strangle reading of the butterfly) and the 25-delta rows of E1: their "
     "regime means, differences, standard errors and bootstrap intervals of phi (0.613, 0.341, -0.272; 0.781, 0.437, -0.345; "
     "0.595, 0.330, -0.265), E4 moving phi by less than 4 per cent, and the phi, E2 and skew-term cells of the 25-delta and "
     "smile-strangle rows of Table 6",
     "the 25-delta and smile-reading skew-price series are withheld: with them a worst-case attacker could narrow single "
     "risk-reversal quotes, without them it cannot (scripts/check_public_reversibility.py and the public README)"),
    ("Table 6 (Section 6)", "Robustness variants other than the base, 25-delta and ATM rows (Fenics quotes, stale "
     "butterflies missing, the two exclusions, implementable returns, dollar carry, ten currencies, log returns, "
     "previous-day spot, vanna-volga), the hedged return (8.7, t 0.73) and theta_UB (0.51) of the smile-strangle row, and "
     "the counts behind them (109 stale currency-months, 43 months, 40 months)",
     "each needs quotes, per-leg premia or a different set of legs; the hedged mean and theta_UB of the 25-delta row "
     "are reproduced"),
    ("Section 6.2", "Vanna-volga smiles (99.7% arbitrage-free, held-out errors 0.10 and 0.13 against 0.13 and 0.16)",
     "per-currency smiles"),
    ("Table 7 (Section 6.3)", "Three-month tenor", "the three-month series are not published"),
    ("Section 3.2", "Substitution and exclusion counts and the coverage of the long sample's series before 1995",
     "quote-level audit records"),
]


# ---------------------------------------------------------------------------
# My private summaries (optional)


@dataclass
class Ref:
    value: object
    half_ulp: float
    source: str


def _num(tok: str):
    try:
        return float(tok)
    except ValueError:
        return tok


def _half_ulp_sig(v, sig: int) -> float:
    if not isinstance(v, float) or not math.isfinite(v) or v == 0:
        return 0.0 if not isinstance(v, float) or v != 0 else 0.5 * 10.0 ** (-sig)
    return 0.5 * 10.0 ** (math.floor(math.log10(abs(v))) - sig + 1)


def _half_ulp_text(tok: str) -> float:
    m = re.fullmatch(r"[-+]?\d*\.?(\d*)(?:[eE]([-+]?\d+))?", tok.strip())
    if not m:
        return 0.0
    return 0.5 * 10.0 ** (int(m.group(2) or 0) - len(m.group(1)))


def _block(text: str, heading: str) -> list[str]:
    """Lines of the first fenced block after a heading that starts with `heading`."""
    i = text.index(heading)
    a = text.index("```", i) + 3
    b = text.index("```", a)
    return [ln for ln in text[a:b].splitlines() if ln.strip()]


def load_private(d: Path) -> dict:
    refs: dict[str, Ref] = {}

    def put(key, v, half, src):
        refs[key] = Ref(v, half, src)

    e1 = pd.read_csv(d / "e1_summary.csv")
    for r in e1.itertuples(index=False):
        for f in ("n_zero_rate", "mean_zero_rate", "se_zero_rate", "n_hiking", "mean_hiking", "se_hiking", "difference",
                  "se_difference_hac", "hac_lag", "boot_lo", "boot_hi"):
            put(f"e1.{r.reading}.{dlab(r.delta)}.{r.variable}.{f}", getattr(r, f), 0.0, "e1_summary.csv")
    md = (d / "e1_summary.md").read_text()
    for ln in _block(md, "## Supplementary"):
        k, v = ln.split()
        put(f"e1.supp.{k}", float(v), _half_ulp_sig(float(v), 6), "e1_summary.md")
    ext = _block(md, "## Extended sample")
    cols = ext[0].split()
    mean = dict(zip(cols, map(float, next(ln for ln in ext if ln.split()[0] == "mean").split()[1:])))
    for c in ("C_skew", "FD", "phi"):
        put(f"e1.ext.mean_{c}", mean[c], _half_ulp_sig(mean[c], 6), "e1_summary.md")

    s4 = (d / "stage4_summary.md").read_text()
    lines = _block(s4, "## E2")
    head = lines[0].split()
    for ln in lines[1:]:
        tok = ln.split()
        if tok[0] == "phi":  # the rr_or rows are not compared: rr_or is not published
            for c, v in zip(head, tok[1:]):
                if v != "NaN" and c != "value":
                    put(f"e2.{tok[0]}.{c}", float(v), _half_ulp_sig(float(v), 6), "stage4_summary.md")
    for ln in _block(s4, "## E3")[1:]:
        k, a, b = ln.split()
        for sample, v in (("primary", a), ("extended", b)):
            v = _num(v)
            put(f"e3.{sample}.{k}", v, _half_ulp_sig(v, 6), "stage4_summary.md")

    mom = (d / "moments_summary.md").read_text()
    lines = _block(mom, "## E2 secondary")
    head = lines[0].split()[:11]
    for ln in lines[1:]:
        tok = ln.split()[:11]
        for c, v in zip(head[2:], tok[2:]):
            put(f"mom.{tok[0]}.{tok[1]}.{c}", float(v), _half_ulp_sig(float(v), 4), "moments_summary.md")

    ra = pd.read_csv(d / "regime_attribution.csv")
    for r in ra.to_dict("records"):
        for k, v in r.items():
            if k != "variable":
                put(f"ra.{r['variable']}.{k}", v, 0.0, "regime_attribution.csv")
        put(f"ra.{r['variable']}.t_ar1", r["difference"] / r["se_difference_ar1"], 0.0, "regime_attribution.csv")
    rr = pd.read_csv(d / "regime_attribution_ratio.csv").iloc[0].to_dict()
    for k, key in (("ratio_dC_over_dFD", "ratio"), ("set_nw", "set_nw"), ("set_nw_lo", "set_nw_lo"), ("set_nw_hi", "set_nw_hi"),
                   ("set_boot", "set_boot"), ("set_boot_lo", "set_boot_lo"), ("set_boot_hi", "set_boot_hi"),
                   ("one_for_one_diff_bp", "one_for_one_diff_bp"), ("one_for_one_se_nw", "one_for_one_se_nw"),
                   ("one_for_one_t_nw", "one_for_one_t_nw"), ("one_for_one_p_boot_two_sided", "one_for_one_p_boot"),
                   ("corr_C_FD_months", "corr_C_FD_months")):
        put(f"ra.{key}", rr[k], 0.0, "regime_attribution_ratio.csv")

    e5 = (d / "e5_summary.md").read_text()
    lines = _block(e5, "# E5")
    head = lines[0].split()
    for ln in lines[1:]:
        k, *vals = ln.split()
        for sample, a in zip(head, vals):
            v = _num(a)
            put(f"e5.{sample}.{k}", v, _half_ulp_sig(v, 6), "e5_summary.md")
    for k in ("corr_U_rvix", "corr_H10_rvix"):
        refs[f"e5.{k}"] = refs[f"e5.primary.{k}"]
    put("e5.n_r_vix", int(refs["e5.primary.n"].value), 0.0, "e5_summary.md")
    m = re.search(r"AR\(1\) coefficient of the FX-volatility series: (\d+(?:\.\d+)?)", e5)
    put("e5.ar1_rho", float(m.group(1)), _half_ulp_text(m.group(1)), "e5_summary.md")

    la = (d / "long_atm_summary.md").read_text()
    labels = sorted(LA_WINDOWS, key=len, reverse=True)
    ret_cols = ["months", "first", "last", "U mean", "U se", "U t", "Hatm mean", "Hatm se", "Hatm t", "Hatm-U mean", "Hatm-U se"]
    theta_cols = ["theta_UB", "CS", "CS lo", "CS hi", "theta0", "premium mean", "premium se", "payoff mean", "payoff se",
                  "U lag", "Hatm lag"]
    for heading, cols in (("## HML^U and ATM-hedged", ret_cols), ("## theta_UB, hedge premium", theta_cols)):
        for ln in _block(la, heading)[1:]:
            label = next(x for x in labels if ln.startswith(x))
            tok = ln[len(label):].split()
            if cols is theta_cols:
                tok = [tok[0], " ".join(tok[1:-9]), *tok[-9:]]
            for c, v in zip(cols, tok):
                if c == "theta0":
                    continue
                v = {"+inf": np.inf, "-inf": -np.inf}.get(v, _num(v))
                if c in ("months", "U lag", "Hatm lag"):
                    v = int(v)
                put(f"la.{label}.{c}", v, _half_ulp_sig(v, 4) if isinstance(v, float) else 0.0, "long_atm_summary.md")

    dc = (d / "design_checks_summary.md").read_text()
    dc = dc[dc.index("## 3. Delta"):dc.index("## 4.")]
    hk = {"10-delta": "H10", "25-delta": "H25", "ATM": "Hatm"}
    src = "design_checks_summary.md"
    for m in re.finditer(r"^\| (10-delta|25-delta|ATM) \| (-?\d+(?:\.\d+)?) \| (\d+(?:\.\d+)?) \| (\d+) \| \[(-?\d+(?:\.\d+)?), (-?\d+(?:\.\d+)?)\] \| "
                         r"(True|False) \| (True|False|) \| (.+?) \| (True|False) \| (True|False|) \|$", dc, re.M):
        h = hk[m.group(1)]
        for i, f in ((2, "theta"), (3, "se_delta"), (5, "lo"), (6, "hi")):
            put(f"dc.{h}.{f}", float(m.group(i)), _half_ulp_text(m.group(i)), src)
        put(f"dc.{h}.lag", int(m.group(4)), 0.0, src)
        put(f"dc.{h}.delta_contains_0", m.group(7) == "True", 0.0, src)
        put(f"dc.{h}.acc0", m.group(10) == "True", 0.0, src)
        if m.group(8):
            put(f"dc.{h}.delta_contains_theta0", m.group(8) == "True", 0.0, src)
        s = re.search(r"grid minus \[(-?\d+(?:\.\d+)?), (-?\d+(?:\.\d+)?)\]", m.group(9))
        put(f"dc.{h}.set_kind", m.group(9).split(":")[0], 0.0, src)
        if s:
            put(f"dc.{h}.rej_lo", float(s.group(1)), _half_ulp_text(s.group(1)), src)
            put(f"dc.{h}.rej_hi", float(s.group(2)), _half_ulp_text(s.group(2)), src)
    for m in re.finditer(r"(10-delta|25-delta|ATM): mean (-?\d+(?:\.\d+)?) bp a month, t (-?\d+(?:\.\d+)?) by test inversion and "
                         r"(-?\d+(?:\.\d+)?) by the delta method; s.e. of the mean of H - U over that of H - \(1 - theta_UB\)U (\d+(?:\.\d+)?)", dc):
        h = hk[m.group(1)]
        put(f"dc.{h}.mean_cost", float(m.group(2)) / 1e4, _half_ulp_text(m.group(2)) / 1e4, src)
        for i, f in ((3, "t0"), (4, "t_delta"), (5, "se_ratio")):
            put(f"dc.{h}.{f}", float(m.group(i)), _half_ulp_text(m.group(i)), src)
    for m in re.finditer(r"(10-delta|25-delta|ATM): block length (\d+(?:\.\d+)?), 2.5% and 97.5% quantiles of t\* (-?\d+(?:\.\d+)?) and "
                         r"(-?\d+(?:\.\d+)?), 95% interval for the mean \[(-?\d+(?:\.\d+)?), (-?\d+(?:\.\d+)?)\] bp, share of t\* at or below the "
                         r"observed t of (-?\d+(?:\.\d+)?): (\d+(?:\.\d+)?)", dc):
        h = hk[m.group(1)]
        for i, f, sc in ((2, "boot_block", 1), (3, "boot_q_lo", 1), (4, "boot_q_hi", 1), (5, "boot_lo", 1e4),
                         (6, "boot_hi", 1e4), (8, "boot_share_below", 1)):
            put(f"dc.{h}.{f}", float(m.group(i)) / sc, _half_ulp_text(m.group(i)) / sc, src)
    for m in re.finditer(r"(10-delta|25-delta|ATM): skewness of H - U (-?\d+(?:\.\d+)?), t after dropping the month with the "
                         r"largest H - U (-?\d+(?:\.\d+)?)", dc):
        h = hk[m.group(1)]
        put(f"dc.{h}.skew", float(m.group(2)), _half_ulp_text(m.group(2)), src)
        put(f"dc.{h}.t_drop", float(m.group(3)), _half_ulp_text(m.group(3)), src)
    m = re.search(r"retains theta_0 = (\d+(?:\.\d+)?) at the automatic lag \(t (-?\d+(?:\.\d+)?), lag (\d+)\); the fixed lags from 0 to 12 "
                  r"at which the test rejects it: ([\d, ]+|none)", dc)
    put("dc.H10.theta0", float(m.group(1)), _half_ulp_text(m.group(1)), src)
    put("dc.H10.acc_th0", True, 0.0, src)
    put("dc.H10.t_th0", float(m.group(2)), _half_ulp_text(m.group(2)), src)
    put("dc.H10.lag_th0", int(m.group(3)), 0.0, src)
    put("dc.H10.lags_rej_th0", [] if m.group(4) == "none" else [int(x) for x in m.group(4).split(",")], 0.0, src)
    m = re.search(r"the t of the mean of H10 - U is (-?\d+(?:\.\d+)?)", dc)
    put("dc.ext_t", float(m.group(1)), _half_ulp_text(m.group(1)), src)
    for m in re.finditer(r"^\| (10-delta|25-delta|ATM), v = (0|theta_0) \| (.+) \|$", dc, re.M):
        cells = m.group(3).split(" | ")
        kind = "tfix0" if m.group(2) == "0" else "tfixth0"
        for k, c in zip(FIXED_LAGS, cells):
            put(f"dc.{hk[m.group(1)]}.{kind}.L{k}", float(c), _half_ulp_text(c), src)

    vd = d / "verdelhan_comparison.md"
    if vd.exists():
        t = vd.read_text()
        src = "verdelhan_comparison.md"
        put("vd.n", int(re.search(r"Overlap: (\d+) months", t).group(1)), 0.0, src)
        for pat, keys in ((r"Mean HML\^U \(bp per month, Newey.West s.e.\) \| (-?\d+(?:\.\d+)?) \((-?\d+(?:\.\d+)?)\)", ("vd.mean_u_bp", "vd.se_u_bp")),
                          (r"Mean developed HML \(bp per month, Newey.West s.e.\) \| (-?\d+(?:\.\d+)?) \((-?\d+(?:\.\d+)?)\)", ("vd.mean_h_bp", "vd.se_h_bp")),
                          (r"Mean difference \(bp per month, Newey.West s.e.\) \| (-?\d+(?:\.\d+)?) \((-?\d+(?:\.\d+)?)\)", ("vd.mean_d_bp", "vd.se_d_bp")),
                          (r"Correlation of monthly returns \| (-?\d+(?:\.\d+)?)", ("vd.corr",)),
                          (r"Slope of HML\^U on their HML \(Newey.West s.e.\) \| (-?\d+(?:\.\d+)?) \((-?\d+(?:\.\d+)?)\)", ("vd.slope", "vd.slope_se")),
                          (r"Standard deviation, HML\^U and theirs \(bp per month\) \| (\d+(?:\.\d+)?), (\d+(?:\.\d+)?)", ("vd.sd_u_bp", "vd.sd_h_bp"))):
            m = re.search(pat, t)
            for i, k in enumerate(keys, start=1):
                put(k, float(m.group(i)), _half_ulp_text(m.group(i)), src)
        m = re.search(r"Share of months with the same sign \| (\d+(?:\.\d+)?)%", t)
        put("vd.same_sign", float(m.group(1)) / 100, _half_ulp_text(m.group(1)) / 100, src)

    rob = pd.read_csv(d / "robustness_summary.csv", index_col=0)
    for v in ROB_VARIANTS:
        for k, x in rob.loc[v].items():
            if k in ("n_phi", "n_returns", "arb_fail_legs", "CW_n"):
                x = int(x) if np.isfinite(x) else x
            put(f"rob.{v}.{k}", x, 0.0, "robustness_summary.csv")

    arch = d / "bootstrap_arch_check.md"
    if arch.exists():
        t = arch.read_text()
        lines = _block(t, "Interval bounds")
        head = lines[0].split()
        for ln in lines[1:]:
            tok = ln.split()
            key = f"arch.{tok[0]}.{dlab(float(tok[1]))}.{tok[2]}"
            for c, v in zip(head[3:], tok[3:]):
                if c in ("lo_mine", "hi_mine", "lo_arch", "hi_arch", "b_mine_zero_rate", "b_mine_hiking"):
                    put(f"{key}.{c}", float(v), _half_ulp_sig(float(v), 4), "bootstrap_arch_check.md")
        m = re.search(r"Largest difference between interval bounds: φ (\d+(?:\.\d+)?); C_skew and FD (\d+(?:\.\d+)?) bp", t)
        put("arch.max_phi", float(m.group(1)), _half_ulp_text(m.group(1)), "bootstrap_arch_check.md")
        put("arch.max_bp", float(m.group(2)), _half_ulp_text(m.group(2)), "bootstrap_arch_check.md")
    return refs


# ---------------------------------------------------------------------------
# Comparisons and report


def _decimals(text: str) -> int:
    return len(text.split(".")[1]) if "." in text else 0


def _is_p(key: str) -> bool:
    return any(s in key for s in ("p_one_sided", "p_boot", "_p_", "E2_p", "CW_p", "share_below", "percentile"))


def _draws(key: str) -> int:
    return B_GRID if key.startswith(("mom.", "rob.")) else B_MAIN


def check_paper(e: Expect, R: dict) -> tuple[str, str]:
    if e.key not in R:
        return "missing", ""
    v = R[e.key]
    if e.kind == "exact":
        return ("pass" if v == e.printed else "FAIL"), str(v)
    x = float(v) * e.scale
    p = float(e.printed)
    if e.kind == "lt":
        return ("pass" if x < p else "FAIL"), f"{x:.3g}"
    dec = _decimals(e.printed)
    shown = f"{x:.{dec + 2}f}"
    if e.kind == "le":
        return ("pass" if round(x, dec) <= p + 1e-12 else "FAIL"), shown
    ok = abs(x - p) <= 0.5 * 10.0 ** (-dec) * (1 + 1e-9) + 1e-12
    return ("pass" if ok else "FAIL"), shown


def private_scale(key: str, R: dict) -> float:
    """Magnitude against which the rounding of the public series is judged.

    A regime difference, or a statistic built from one, can be much smaller than
    the two means it subtracts, and the rounding error of the inputs is
    proportional to those means, not to the difference. For such keys the scale
    is the larger regime mean (divided by the standard error for a t statistic,
    and by the change in FD for the Wald ratio); otherwise it is the value itself.
    """
    stem, _, last = key.rpartition(".")
    if last in ("difference", "t_ar1") and f"{stem}.mean_zero_rate" in R:
        m = max(abs(R[f"{stem}.mean_zero_rate"]), abs(R[f"{stem}.mean_hiking"]))
        return m if last == "difference" else m / R[f"{stem}.se_difference_ar1"]
    if key == "ra.ratio":
        return abs(R["ra.C_skew_bp.mean_zero_rate"] / R["ra.FD_bp.difference"])
    return 0.0


def check_private(key: str, v, ref: Ref, scale: float = 0.0) -> tuple[bool, float]:
    """Whether the recomputed value agrees with the private summary; returns the absolute difference."""
    r = ref.value
    if isinstance(r, (bool, list, str)) or isinstance(v, (bool, list, str)):
        return (v == r), 0.0
    v, r = float(v), float(r)
    if not (math.isfinite(v) and math.isfinite(r)):
        return (v == r) or (math.isnan(v) and math.isnan(r)), 0.0
    rtol = PRIVATE_RTOL_EXTENDED if ("extended" in key or ".ext." in key) else PRIVATE_RTOL
    tol = ref.half_ulp + rtol * max(abs(r), scale) + 1e-12
    if _is_p(key):
        tol += PRIVATE_P_DRAWS / _draws(key)
    return abs(v - r) <= tol * (1 + 1e-9), abs(v - r)


def _fmt(v) -> str:
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def run_tasks(data: dict, s4: pd.DataFrame, cboe: Path | None, xls: Path | None, have_arch: bool, workers: int):
    p, e, la = data["primary"], data["extended"], data["long_atm"]
    tasks = [(f"E1 {r} {dlab(d)}", task_e1_combo, (p, r, d, B_MAIN)) for r, d in COMBOS]
    tasks += [("E1 supplementary", task_e1_supplementary, (p, e, B_MAIN)),
              ("realised-return split", task_realised_split, (s4, B_MAIN)),
              ("regime attribution", task_regime_attribution, (p, B_MAIN)),
              ("E2 phi", task_e2, (s4, "phi", B_MAIN)),
              ("moment predictors", task_moments, (p, s4, B_GRID)), ("E3", task_e3, (s4,))]
    tasks += [(f"design checks {h}", task_design, (s4, h, B_MAIN)) for h in HEDGES]
    tasks += [(f"long ATM {w}", task_long_atm, (la, w)) for w in LA_WINDOWS]
    tasks += [(f"robustness {v}", task_robustness, (p, v, B_GRID)) for v in ROB_VARIANTS]
    if have_arch:
        tasks += [(f"arch {r} {dlab(d)}", task_arch, (p, r, d)) for r, d in COMBOS]
    if cboe is not None:
        tasks.append(("E5", task_e5, (s4, str(cboe))))
    if xls is not None:
        tasks.append(("published portfolios", task_verdelhan, (p, str(xls))))
    R, timing = {}, {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [(name, ex.submit(_timed, f, args)) for name, f, args in tasks]
        for name, fut in futs:
            res, secs = fut.result()
            R.update(res)
            timing[name] = secs
    return R, timing


def _timed(f, args):
    t = time.perf_counter()
    res = f(*args)
    return res, time.perf_counter() - t


def _latest(pattern: str) -> Path | None:
    hits = sorted(ROOT.glob(pattern))
    return hits[-1] if hits else None


def main():
    ap = argparse.ArgumentParser(description="Recompute the paper's portfolio-level results from the public series.")
    ap.add_argument("--public-dir", type=Path, default=PUBLIC)
    ap.add_argument("--cboe-dir", type=Path, default=None, help="folder of Cboe VX settlement files (acquire_cboe_vx.py)")
    ap.add_argument("--verdelhan-file", type=Path, default=None, help="CurrencyPortfolios.xls (acquire_verdelhan.py)")
    ap.add_argument("--private-dir", type=Path, default=PRIVATE, help="my private results, compared if present")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "private" / "results" / "public_reproduction.md")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    for name in ("public_dir", "cboe_dir", "verdelhan_file", "private_dir", "out"):
        if getattr(args, name) is not None:
            setattr(args, name, getattr(args, name).resolve())
    t0 = time.perf_counter()
    workers = max(1, min(args.workers, 8))

    try:
        data = load_public(args.public_dir)
    except (FileNotFoundError, ValueError) as exc:
        print(f"reproduce_from_public: {exc}", file=sys.stderr)
        sys.exit(2)
    cboe = args.cboe_dir or _latest("data/private/cboe/*/raw")
    xls = args.verdelhan_file or _latest("data/private/verdelhan/*/CurrencyPortfolios.xls")
    for name, path in (("--cboe-dir", args.cboe_dir), ("--verdelhan-file", args.verdelhan_file)):
        if path is not None and not path.exists():
            print(f"reproduce_from_public: {name} {path} does not exist", file=sys.stderr)
            sys.exit(2)
    cboe = cboe if cboe is not None and cboe.exists() else None
    have_xlrd = importlib.util.find_spec("xlrd") is not None
    xls = xls if xls is not None and xls.exists() and have_xlrd else None
    have_arch = importlib.util.find_spec("arch") is not None

    s4 = stage4_frame(data)
    R, timing = run_tasks(data, s4, cboe, xls, have_arch, workers)
    R.update(counts(data))
    if "e5.ar1_rho" not in R:  # the AR(1) of the FX-volatility level needs no Cboe file
        R["e5.ar1_rho"] = fx_vol_innovation(s4)[1]
    R.update(derived(R))

    private = args.private_dir if (args.private_dir / "e1_summary.csv").exists() else None
    refs = load_private(private) if private else {}

    paper_rows = []
    for e in PAPER:
        status, shown = check_paper(e, R)
        if status == "missing":
            continue  # an optional input was absent; listed below
        priv = ""
        if e.key in refs:
            ok, _ = check_private(e.key, R[e.key], refs[e.key], private_scale(e.key, R))
            rv = refs[e.key].value
            rv = rv * e.scale if isinstance(rv, (float, int, np.floating, np.integer)) and not isinstance(rv, bool) else rv
            priv = f"{_fmt(float(rv)) if isinstance(rv, (float, np.floating)) else rv} ({'agrees' if ok else 'DIFFERS'})"
        paper_rows.append((e, status, shown, priv))
    priv_rows = []
    for key, ref in refs.items():
        if key not in R:
            continue
        ok, diff = check_private(key, R[key], ref, private_scale(key, R))
        priv_rows.append((key, R[key], ref, ok, diff))
    missing = [e for e in PAPER if e.key not in R]

    paper_fail = [r for r in paper_rows if r[1] != "pass"]
    priv_fail = [r for r in priv_rows if not r[3]]
    elapsed = time.perf_counter() - t0

    L = ["# Reproduction of the paper's portfolio-level results from the public series", "",
         "Written by `scripts/reproduce_from_public.py`. Pooled statistics only.", "",
         "## Inputs", "",
         f"- Public series: `{args.public_dir.relative_to(ROOT) if args.public_dir.is_relative_to(ROOT) else args.public_dir}`; "
         f"SHA-256 of {', '.join(data['manifest']['files'])} checked against manifest.json; theta0_mean_10d = {data['theta0']:.6g} "
         "(README.md).",
         f"- Cboe VX files: {cboe.relative_to(ROOT) if cboe and cboe.is_relative_to(ROOT) else cboe or 'absent, E5 not computed'}.",
         f"- Verdelhan's portfolios: {xls.relative_to(ROOT) if xls and xls.is_relative_to(ROOT) else xls or ('absent' if have_xlrd else 'xlrd not installed') + ', comparison not computed'}.",
         f"- arch package: {'present' if have_arch else 'absent, the arch check is not computed'}.",
         f"- Private summaries: {private.relative_to(ROOT) if private and private.is_relative_to(ROOT) else private or 'absent, no private comparison'}.",
         f"- Python {platform.python_version()}, numpy {np.__version__}, pandas {pd.__version__}; {workers} worker processes; "
         f"bootstrap draws {B_MAIN:,} and {B_GRID:,} as in the paper.", "",
         "## Outcome", "",
         f"- Paper numbers compared: {len(paper_rows)}; matching at the printed precision: "
         f"{len(paper_rows) - len(paper_fail)}; not matching: {len(paper_fail)}; not computed because an optional "
         f"input is absent: {len(missing)}.",
         f"- Private summary values compared: {len(priv_rows)}; agreeing: {len(priv_rows) - len(priv_fail)}; "
         f"differing: {len(priv_fail)}.",
         f"- Results listed as not reproducible from the public series: {len(NOT_REPRODUCIBLE)} (below).",
         f"- Run time: {elapsed:.1f} s.", "",
         f"Result: {'PASS' if not paper_fail and not priv_fail else 'FAIL'}.", ""]
    if priv_rows:
        rel = [(r[4] / abs(float(r[2].value)), r[0]) for r in priv_rows if r[3] and isinstance(r[2].value, float)
               and r[2].half_ulp == 0.0 and r[2].value != 0 and math.isfinite(r[2].value)]
        if rel:
            worst = max(rel)
            L += [f"Against the full-precision private CSV summaries, the largest relative difference is {worst[0]:.2e} "
                  f"({worst[1]}).", ""]
    L += ["## Paper numbers", "",
          "Recomputed values are shown to two digits beyond the paper's precision; private summary values are in the "
          "paper's units. Status: pass or FAIL.", "",
          "| Paper location | Result | Paper | Recomputed | Private summary | Status |", "| --- | --- | --- | --- | --- | --- |"]
    for e, status, shown, priv in paper_rows:
        prefix = {"lt": "< ", "le": "<= "}.get(e.kind, "")
        L.append(f"| {e.where} | {e.label} | {prefix}{e.printed} | {shown} | {priv} | {status} |")
    if priv_fail:
        L += ["", "## Private summary values that differ", "", "| Key | Recomputed | Private | Source |",
              "| --- | --- | --- | --- |"]
        for key, v, ref, _, _ in priv_fail:
            L.append(f"| {key} | {_fmt(v)} | {_fmt(ref.value)} | {ref.source} |")
    if missing:
        L += ["", "## Paper numbers not computed (optional input absent)", ""]
        L += [f"- {e.where}: {e.label}" for e in missing]
    L += ["", "## Not reproducible from the public series", "", "| Paper location | Result | Reason |", "| --- | --- | --- |"]
    L += [f"| {w} | {what} | {why} |" for w, what, why in NOT_REPRODUCIBLE]
    L += ["", "## Run time by task (seconds, in parallel)", "", "| Task | Seconds |", "| --- | --- |"]
    L += [f"| {k} | {v:.1f} |" for k, v in sorted(timing.items(), key=lambda kv: -kv[1])]
    L.append("")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(L))

    print(f"Paper numbers: {len(paper_rows)} compared, {len(paper_rows) - len(paper_fail)} matching, "
          f"{len(paper_fail)} not matching, {len(missing)} not computed (optional input absent)")
    print(f"Private summary values: {len(priv_rows)} compared, {len(priv_rows) - len(priv_fail)} agreeing, "
          f"{len(priv_fail)} differing")
    print(f"Not reproducible from the public series: {len(NOT_REPRODUCIBLE)} results (listed in the report)")
    for e, status, shown, _ in paper_fail:
        print(f"  FAIL {e.where}: {e.label}: paper {e.printed}, recomputed {shown}")
    for key, v, ref, _, _ in priv_fail:
        print(f"  DIFFERS {key}: recomputed {_fmt(v)}, private {_fmt(ref.value)} ({ref.source})")
    print(f"Report: {args.out}; {elapsed:.1f} s")
    sys.exit(1 if paper_fail or priv_fail else 0)


if __name__ == "__main__":
    main()
