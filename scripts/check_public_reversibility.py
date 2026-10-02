"""Month-by-month check that the public portfolio-level series cannot be traced back to the LSEG quotes.

Purpose. I publish month-level series at portfolio level
(data/public/fx_carry_portfolio_series/primary.csv, extended.csv and
long_atm.csv) so that a reader can rerun most of the paper's portfolio-level tests
without an LSEG licence. The licence holder confirmed on 1 October 2026 that
manipulated or transformed data may be published and raw data may not. I
publish only output that meets the test of LSEG's Redistribution Guidelines
for compliant derived data: it must be "unrecognizable, non-reversible, and
cannot be traced back to the original content without exceptional effort".
This script checks that test month by month for the final publication set. It
reads private data, writes a pooled summary (public_reversibility.md) and a
month-level diagnostic table without any quote
(public_reversibility_months.csv) to data/private/results/2026-09-23/, and
never writes into data/public/.

Final publication set, rounded to six significant figures:

- primary.csv: C_skew and phi of the market reading at 10 delta (E1), FD, U,
  H10, H25, Hatm, c1, c2, c3, var_lo, var_hi, oskew_lo, oskew_hi and sigma_fx;
- extended.csv: C_skew and phi (market reading, 10 delta), FD, U, H10, H25,
  Hatm, c1, c2, c3 and sigma_fx;
- long_atm.csv: n_legs, U and Hatm.

Leg identities and per-currency values are never published. The skew price
and phi of the other three E1 combinations (market reading at 25 delta, smile
reading at 10 and 25 delta) are withheld: with them, the linear attacker of
part (f) could narrow single risk-reversal quotes to within a quarter of a
volatility point in about a fifth of primary legs, and without them it cannot.
Two columns of an earlier build are also withdrawn: rr_or (a plain signed average of the six
quoted 10-delta risk reversals in quote units, the least transformed column)
and the long at-the-money payoff and premium (in some months exactly one
leg's ATM option ends in the money, so the payoff would depend on that leg
alone; Hatm - U always mixes all legs).

Adversary. I assume the worst case: the adversary knows everything about a
month except the smile quotes, that is, the leg identities and orientations,
spot, forward points, rates, dates, the spot at expiry and the realised
volatility forecast, and knows the code. The unknowns of a month are the
calibration quotes q = (ATM, RR25, BF25) of each leg (three per leg). In the
long at-the-money sample they are the ATM volatilities only. A real reader
does not know the leg identities, so this adversary is stronger than any
reader of the files.

Published smile-dependent values. Every published smile-dependent number of a
month is an average over the legs, with orientation, of per-leg contributions,
and each leg's contribution depends only on its own quotes. Exact identities
add no information: c3 = -C_skew (market reading, 10-delta), H10 = U + c1 +
c2 + c3 and phi = C_skew/FD, with U and FD known to the adversary; and every
long at-the-money month from the primary start has the same legs and quotes
as the primary month, so its U and Hatm repeat the primary month's. The
independent values are:

- primary sample (composite quotes, 3 and 3 legs): C_skew (market, 10-delta),
  H25, Hatm, c1, c2 and the setting (i) portfolio moment predictors var_lo,
  var_hi, oskew_lo and oskew_hi. That is 9 values for 18 unknowns; the last
  formation month has no realised return and gives 5;
- extended sample (Fenics quotes, 2 and 2 legs): C_skew (market, 10-delta),
  H25, Hatm, c1 and c2, so 5 values for 12 unknowns;
- long at-the-money sample before the primary start: Hatm, 1 value for 4 or
  6 ATM unknowns.

FD and U depend on forward points and spot, and sigma_fx (primary and
extended) on daily spot only; part (d) treats them. theta0 is not published
per month, only its mean over months. In the constructive part (c) I
nevertheless hold each month's theta0 (10-delta, stage 4) and the ATM-hedge
theta0 (long at-the-money file) fixed, so that any mean of them that I
publish is reproduced as well. Fenics and composite quotes
are distinct LSEG content and are treated as distinct unknowns.

Method.

(a) Per-leg forward map. For each leg, leg_terms maps its quotes to its
    contributions to every month-level quantity with the project's own
    functions: SABR calibration under both butterfly readings
    (qef.fx.smile.calibrate_sabr, warm-started as scripts/calibrate_smiles.py
    does), protective-option strikes and premia (qef.fx.crash.leg_skew_cost),
    payoffs, forward returns and the ATM hedge as in scripts/estimate_stage4.py
    and scripts/estimate_long_atm.py, and the moment intervals of
    scripts/estimate_moments.py (qef.fx.moments.implied_moment_intervals, with
    the leg averages of portfolio_predictors). At the true quotes the averaged
    contributions must reproduce the private month-level values of
    e1_series.csv, stage4_months.csv, long_atm_months.csv and the portfolio
    predictors; the summary reports the largest relative gap. I run this
    validation on every month before anything else (--validate-only runs it
    alone).
(b) Counting and local rank. For each month the central-difference Jacobian
    of the independent published values with respect to the unknowns (step
    1e-5 in volatility; steps of 1e-4 and 1e-6 are compared on every eighth
    month) is scaled by rows to units of the publication's rounding (the half
    unit of the sixth significant figure) and by columns to volatility
    points. I report its singular values, its numerical rank and the
    dimension of its null space.
(c) Constructive non-reversibility. For each month and each target distance
    d of 0.5, 1 and 2 volatility points (and 4 if 2 succeeds), a linear
    programme picks a direction v in the null space of the Jacobian with
    v_j = +-1 for one quote j, |v_i| <= 1 for all i, the smallest weighted
    L1 norm, and a linear prediction q + d v inside the plausibility box. From
    q + d v, chord Gauss-Newton steps with the pseudo-inverse (quote j held
    fixed, so the max-norm distance stays at least d; the Jacobian is
    recomputed when progress stalls) return to the level set until every
    independent value and every pinned theta0 is within 0.01 half units of
    its true value. The alternative is accepted only if, recomputed with the
    pipeline's own warm starts, every published column of the month (the
    long_atm.csv U and Hatm included for months from the primary start)
    rounds to the same six significant figures as the true one, both SABR
    calibrations
    of every leg converge (largest residual below 1e-8) and pass both
    arbitrage checks of scripts/calibrate_smiles.py (+-4 and +-10 ATM
    standard deviations), the moment intervals have status ok, ATM is
    positive and no larger than the largest ATM seen in the sample, |RR25|
    does not exceed the largest |RR25| seen in the sample and BF25 lies
    within the sample's range. Up to four directions are tried
    per target. When the full run would take too long, --construct-sample N
    runs this part on N months per sample, evenly spaced in formation date
    (a sample stratified by date), and the summary says so.
(d) FD and U are two equations in the six forward points when spot is
    known, but forward points also enter every smile-dependent value
    (through the strikes, the payoffs and D_b = F D_q / S). I therefore check, on six primary months, the joint
    problem with both the smile quotes and the forward points unknown, and
    report the rank of the forward block for an adversary who already holds
    the smile quotes. sigma_fx is the window mean of the cross-sectional mean
    absolute daily log spot change: reversing the order of each currency's
    daily changes between dates whose spot enters any other published value
    leaves every published value unchanged and moves the intermediate spots;
    I check this on months of both the primary and the extended sample.
    rr_or is no longer published (see above) and is not analysed.
(e) Formal argument. A continuous map f from an open set U of R^n into R^m
    with m < n is not injective. If it were, composing it with the inclusion
    of R^m in R^n as R^m x {0} would give a continuous injective map from U
    into R^n whose image lies in R^m x {0}, which has empty interior in R^n;
    invariance of domain says that such an image is open, a contradiction.
    Each month's published values are such a map of the quotes, with m < n
    in every sample, so the per-month inverse problem has no unique solution
    whatever the adversary's effort, and rounding to six significant figures
    only enlarges the set of quote vectors consistent with them. Where the
    Jacobian has full row rank m, the implicit function theorem makes that
    set a smooth manifold of dimension n - m near the truth; part (c) shows
    that it reaches quotes that differ materially from the truth, not just
    infinitesimally.
(f) The best linear attacker. Given the exact Jacobian of a month's published
    values at the true quotes, which no real attacker has, the legs, spot and
    forwards, and a normal prior per quote type with the spread of a
    public-information predictor (ATM: the realised-volatility forecast;
    RR25: zero; BF25: the pooled median), the posterior standard deviation of
    each quote measures how closely the published values pin it down; it is
    computed for every month.

Discreteness. Real quotes lie on a grid: composite quotes on multiples of
0.0005 volatility points and Fenics quotes on multiples of about 0.0125 (the
grid is measured from the data and reported). An adversary who knows the grid
looks only for grid points in the thin set of quote vectors consistent with
the rounded values. The argument of (e) is for continuous unknowns and does
not settle that question, so for each month I report the linearised expected
number of grid points in that set inside a null-space box of side 2
volatility points, 2^m L^k / (prod s_i prod g_j), where s_i are the scaled
singular values, k = n - m and g_j the grid step of quote j. A value well
below one means that, on the grid, the truth may be the only consistent point
nearby; the protection then rests on the size of the search, not on
non-uniqueness. I also report how that number changes with fewer significant
figures and when one equation is dropped.

Runtime. On eight processes the whole run takes about a quarter of an hour.
Reproduction, after the scripts that write the inputs (scripts/reproduce.sh,
steps 2 to 7 and 18):

    .venv/bin/python scripts/check_public_reversibility.py --workers 8
"""

from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402
from collections import Counter, defaultdict  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
from dataclasses import dataclass, field  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import linprog  # noqa: E402
from scipy.stats import norm  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from estimate_e1 import PRIMARY_START  # noqa: E402
from estimate_e5 import fx_abs_returns  # noqa: E402
from estimate_moments import portfolio_predictors  # noqa: E402
from estimate_stage4 import daily_spot_mid, realised_vol, spot_on  # noqa: E402

from qef.data.panel import ny_month_ends  # noqa: E402
from qef.data.smile_inputs import build_month_end_inputs, option_dates  # noqa: E402
from qef.fx.arbitrage import check_smile  # noqa: E402
from qef.fx.conventions import G10  # noqa: E402
from qef.fx.crash import forward_discount, forward_return, leg_skew_cost, option_payoff, protective_option, rank_legs  # noqa: E402
from qef.fx.gk import atm_dns_strike, d_plus_minus, forward_premium  # noqa: E402
from qef.fx.moments import implied_moment_intervals  # noqa: E402
from qef.fx.sabr import sabr_vol  # noqa: E402
from qef.fx.smile import SmileQuotes, calibrate_sabr  # noqa: E402

SIG_FIGS = 6  # significant figures of the public file being analysed (set per month in run_month)
SIG_BY_SAMPLE = {"primary": 6, "extended": 4, "long": 6}  # primary.csv, extended.csv, long_atm.csv
WITHHELD_TAGS = ("_m25", "_s10", "_s25")  # E1 combinations whose C_skew and phi are not published
VP = 0.01  # one volatility point in decimal volatility
H_FD = 1e-5  # central-difference step, decimal volatility
H_CHECK = (1e-4, 1e-6)  # steps compared with H_FD
H_LOGF = 1e-6  # central-difference step in log forward (part d)
BP = 1e-4  # one basis point of log forward
TARGETS = (0.5, 1.0, 2.0)  # target max-norm distances, volatility points
STRETCH = 4.0  # tried only when 2 volatility points succeed
FWD_TARGETS = (1.0, 5.0, 10.0)  # forward-point targets of part (d), basis points of log F
GN_TOL = 0.01  # convergence: largest residual in half units of the published rounding
GN_MAX_IT = 40
GN_MAX_REFRESH = 3
MAX_ATTEMPTS = 4  # null-space directions tried per target
RANK_RTOL = 1e-8  # singular values below this share of the largest count as zero
GRID_BOX = 2.0  # side of the null-space box of the grid count, volatility points
CALIB_TOL = 1e-8  # calibrate_sabr success threshold (largest residual, volatility units)
LP_MARGIN = {"atm_lo": 0.5, "rr": 0.05, "bf": 0.02}  # volatility points kept inside the box by the programme
BF_WEIGHT = 3.0  # weight of BF25 moves in the programme's L1 objective
GRID_CANDIDATES = (0.05, 0.025, 0.0125, 0.00625, 0.005, 0.0025, 0.001, 0.0005, 0.00025, 0.0001, 0.00005)
GRID_SHARE = 0.995
COMBOS = (("market", 0.10, "m10"), ("market", 0.25, "m25"), ("smile", 0.10, "s10"), ("smile", 0.25, "s25"))
KEYS = {
    "primary": ("C_m10", "H25", "Hatm", "c1", "c2", "var_lo", "var_hi", "oskew_lo", "oskew_hi"),
    "extended": ("C_m10", "H25", "Hatm", "c1", "c2"),
    "long": ("Hatm",),
}
PINNED = {"primary": ("theta0", "theta0_atm"), "extended": ("theta0",), "long": ("theta0_atm",)}
# published columns whose removal removes each independent equation (identities included)
DROP_COLUMNS = {
    "C_m10": "C_skew and phi (market, 10-delta), c3 and H10", "C_m25": "C_skew and phi (market, 25-delta)",
    "C_s10": "C_skew and phi (smile, 10-delta)", "C_s25": "C_skew and phi (smile, 25-delta)",
    "H25": "H25", "c1": "c1 and H10", "c2": "c2 and H10", "Hatm": "Hatm",
    "var_lo": "var_lo", "var_hi": "var_hi", "oskew_lo": "oskew_lo", "oskew_hi": "oskew_hi",
}
DROP_COLUMNS_PRIM = {"Hatm": "Hatm in primary.csv and in long_atm.csv"}
DROP_COLUMNS_EXT = {"C_m10": "C_skew, phi, c3 and H10", "Hatm": "Hatm"}


# ---------------------------------------------------------------------------
# Month and leg records (private; held in memory only)


@dataclass(frozen=True)
class Leg:
    currency: str
    is_long: bool
    n: int  # legs per side
    S: float
    F: float
    tau: float
    df_quote: float
    df_base: float
    q: tuple  # true (ATM, RR25, BF25), decimal; RR25 and BF25 are NaN in the long sample
    S_end: float  # NaN when the return is not realised
    sig_p: float  # realised volatility forecast, NaN if unavailable (stage 4 then uses ATM)
    warm: tuple  # pipeline warm starts (market, smile), as in calibrate_smiles.py


@dataclass
class Month:
    sample: str  # primary, extended or long
    date: pd.Timestamp
    legs: list
    FD: float
    truth: dict  # published column -> private value
    pinned: dict  # theta0 keys -> private value
    box: dict  # plausibility bounds, volatility points
    flags: dict = field(default_factory=dict)
    others_fd: dict = field(default_factory=dict)  # forward discounts of the other available currencies


def per_leg(sample):
    return 1 if sample == "long" else 3


def fmt_sig(v):
    return "nan" if not np.isfinite(v) else f"{v:.{SIG_FIGS - 1}e}"


def half_unit(v):
    """Half a unit in the last published significant figure of v (zero for v = 0)."""
    if not np.isfinite(v) or v == 0:
        return 0.0
    return 0.5 * 10.0 ** (np.floor(np.log10(abs(v))) - (SIG_FIGS - 1))


# ---------------------------------------------------------------------------
# (a) Per-leg forward map


def leg_terms(leg: Leg, q, sample: str, x0=None, F=None, verify=False):
    """Contributions of one leg to every month-level quantity of its sample.

    q = (ATM, RR25, BF25) in decimal volatility (ATM only is used in the long
    sample). ``x0`` maps a reading to a warm start (log alpha, atanh rho,
    log nu); ``F`` overrides the forward (part d), with D_b = F D_q / S.
    Returns (terms, info); info holds the calibrated parameters and, with
    ``verify``, the convergence, arbitrage and cold-start diagnostics.
    """
    atm, rr, bf = (float(v) for v in q)
    c = leg.currency
    conv = G10[c].delta
    if F is None or F == leg.F:
        F, df_base = leg.F, leg.df_base
    else:
        df_base = F * leg.df_quote / leg.S
    tau, n, w = leg.tau, leg.n, 1.0 / leg.n
    phi = protective_option(c, leg.is_long)
    t, info = {}, {}
    if not atm > 0:
        raise ValueError("ATM volatility not positive")
    # ATM hedge (estimate_stage4 'Hatm' leg and estimate_long_atm)
    k_atm = atm_dns_strike(F, atm, tau, conv)
    v_atm = float(forward_premium(F, k_atm, atm, tau, phi)) / F
    d1, _ = d_plus_minus(F, k_atm, atm, tau)
    t["premium"] = w * v_atm
    t["theta0_atm"] = float(norm.cdf(phi * float(d1))) / (2 * n)
    t["fd"] = (1.0 if leg.is_long else -1.0) * forward_discount(leg.S, F, G10[c].usd_base) / n
    realised = np.isfinite(leg.S_end)
    if realised:
        rx = forward_return(c, leg.is_long, leg.S_end, F)
        pay_atm = option_payoff(c, leg.is_long, k_atm, F, leg.S_end)
        t["U"] = w * rx
        t["payoff"] = w * pay_atm
        t["Hatm"] = w * (rx + pay_atm - v_atm)
    if sample == "long":
        return t, info
    sq = SmileQuotes(F, tau, atm, rr, bf, 0.25, df_base, conv)
    readings = ("market", "smile") if sample == "primary" else ("market",)
    vols = {}
    for reading in readings:
        start = None if x0 is None else x0.get(reading)
        res = calibrate_sabr(sq, reading, x0=start)
        sm = res.smile
        info[f"x_{reading}"] = np.array([np.log(sm.alpha), np.arctanh(sm.rho), np.log(sm.nu)])
        info[f"ok_{reading}"] = bool(res.success)
        info[f"resid_{reading}"] = float(np.max(np.abs(res.residuals)))
        vols[reading] = lambda K, a=sm.alpha, r_=sm.rho, nu=sm.nu: sabr_vol(K, F, tau, a, r_, nu, 1.0)
        if verify:
            sd = atm * np.sqrt(tau)
            arb = []
            for width in (4.0, 10.0):
                rep = check_smile(sm.vol, F, tau, -width * sd, width * sd, n=2001)
                arb.append(rep.convex_ok and rep.slope_bounds_ok and rep.density_ok)
            info[f"arb_{reading}"] = bool(all(arb))
            cold = calibrate_sabr(sq, reading).smile
            info[f"cold_{reading}"] = float(max(abs(cold.alpha / sm.alpha - 1), abs(cold.rho - sm.rho), abs(cold.nu / sm.nu - 1)))
    vol = vols["market"]
    K10, vs10, vf10 = leg_skew_cost(vol, F, tau, df_base, c, atm, leg.is_long, 0.10)
    K25, vs25, vf25 = leg_skew_cost(vol, F, tau, df_base, c, atm, leg.is_long, 0.25)
    t["C_m10"] = (vs10 - vf10) / n
    t["C_m25"] = (vs25 - vf25) / n
    if realised:  # estimate_stage4.month_rows
        sig_p = leg.sig_p if np.isfinite(leg.sig_p) else atm
        vp10 = float(forward_premium(F, K10, sig_p, tau, phi)) / F
        pay10 = option_payoff(c, leg.is_long, K10, F, leg.S_end)
        pay25 = option_payoff(c, leg.is_long, K25, F, leg.S_end)
        t["H10"] = w * (rx + pay10 - vs10)
        t["H25"] = w * (rx + pay25 - vs25)
        t["c1"] = w * (pay10 - vp10)
        t["c2"] = -w * (vf10 - vp10)
        t["c3"] = -w * (vs10 - vf10)
        dd1, _ = d_plus_minus(F, K10, float(vol(K10)), tau)
        t["theta0"] = float(norm.cdf(phi * float(dd1))) / (2 * n)
    if sample == "primary":
        for delta, key in ((0.10, "C_s10"), (0.25, "C_s25")):
            _, a, b = leg_skew_cost(vols["smile"], F, tau, df_base, c, atm, leg.is_long, delta)
            t[key] = (a - b) / n
        mom = implied_moment_intervals(F, tau, vol, G10[c].usd_base, conv, df_base, 0.10)
        info["moment_status"] = mom["status_i"]
        if mom["status_i"] == "ok":  # estimate_moments.portfolio_predictors, setting (i)
            m = 2 * n
            t["var_lo"] = mom["var_lo_i"] / m
            t["var_hi"] = mom["var_hi_i"] / m
            t["oskew_lo"] = (mom["skew_lo_i"] if leg.is_long else -mom["skew_hi_i"]) / m
            t["oskew_hi"] = (mom["skew_hi_i"] if leg.is_long else -mom["skew_lo_i"]) / m
    return t, info


def quotes_of(month: Month, X, i):
    """Quotes of leg i from the unknown vector X (volatility points)."""
    if month.sample == "long":
        return (X[i] * VP, np.nan, np.nan)
    return tuple(X[3 * i: 3 * i + 3] * VP)


def forwards_of(month: Month, X):
    """Forwards when X carries log-forward coordinates (basis points) after the quotes (part d)."""
    nq = per_leg(month.sample) * len(month.legs)
    if len(X) == nq:
        return [None] * len(month.legs)
    return [leg.F * np.exp(X[nq + i] * BP) for i, leg in enumerate(month.legs)]


def month_sums(month: Month, X, x0s=None, verify=False):
    """Month-level quantities: sums of the leg contributions, legs in E1 order."""
    S, infos = defaultdict(float), []
    Fs = forwards_of(month, X)
    for i, leg in enumerate(month.legs):
        t, info = leg_terms(leg, quotes_of(month, X, i), month.sample, None if x0s is None else x0s[i], Fs[i], verify)
        for k, v in t.items():
            S[k] += v
        infos.append(info)
    return S, infos


def published_columns(month: Month, S) -> dict:
    """Every public column of the month that depends on the quotes or forwards, from month sums."""
    FD = S["fd"] if "fd" in S else month.FD
    out = {}
    for col in month.truth:
        if col.startswith("C_skew_"):
            out[col] = S["C_" + col[len("C_skew_"):]]
        elif col.startswith("phi_"):
            out[col] = S["C_" + col[len("phi_"):]] / FD
        elif col == "FD":
            out[col] = FD
        elif col.startswith("long_"):
            out[col] = S[col[len("long_"):]]
        else:
            out[col] = S[col]
    return out


def warm_of(infos, prev=None):
    """Warm starts for the next evaluation: a calibration that did not converge keeps the previous warm start.

    As in calibrate_smiles.py, where the warm start advances only after a converged calibration; a
    diverged calibration (for example nu overflowing) would otherwise seed the next one.
    """
    out = []
    for i, info in enumerate(infos):
        w = {}
        for r in ("market", "smile"):
            if f"x_{r}" not in info:
                continue
            if prev is None or (info[f"ok_{r}"] and info[f"resid_{r}"] < CALIB_TOL):
                w[r] = info[f"x_{r}"]
            elif prev[i] and prev[i].get(r) is not None:
                w[r] = prev[i][r]
        out.append(w)
    return out


# ---------------------------------------------------------------------------
# (b) Jacobian, rank and grid count


def jacobian(month: Month, X, keys, x0s, h=H_FD, fwd=False):
    """Central differences of the month-level quantities ``keys``: per volatility point, and per basis point of log F."""
    nq = per_leg(month.sample) * len(month.legs)
    J = np.zeros((len(keys), len(X)))
    Fs = forwards_of(month, X)
    for i, leg in enumerate(month.legs):
        cols = [i * per_leg(month.sample) + j for j in range(per_leg(month.sample))]
        if fwd:
            cols.append(nq + i)
        for col in cols:
            vals = []
            for sgn in (1.0, -1.0):
                Xp = X.copy()
                Xp[col] += sgn * (h / BP if col >= nq else h / VP)
                F = leg.F * np.exp(Xp[nq + i] * BP) if fwd else Fs[i]
                t, _ = leg_terms(leg, quotes_of(month, Xp, i), month.sample, x0s[i], F)
                vals.append(t)
            step = 2 * (h / BP if col >= nq else h / VP)
            J[:, col] = [(vals[0][k] - vals[1][k]) / step for k in keys]
    return J


def svd_stats(A):
    s = np.linalg.svd(A, compute_uv=False)
    rank = int(np.sum(s > RANK_RTOL * s[0])) if s.size and s[0] > 0 else 0
    return s, rank


def grid_stats(A, grid):
    """Linearised density of grid points consistent with the rounded values.

    Returns (log10 E1, k): E1 is the expected number of grid points of the
    quotes in the set {x : |A (x - x_true)| <= 1 row by row} per unit
    k-volume (volatility points) of its null-space extent, with k = n - rank,
    so a null-space box of side L holds E1 L^k of them, and consistent grid
    points are typically L* = E1^(-1/k) volatility points apart.
    """
    s, rank = svd_stats(A)
    return rank * np.log10(2.0) - np.sum(np.log10(s[:rank])) - np.sum(np.log10(grid)), A.shape[1] - rank


def spacing(log10_e1, k):
    return 10.0 ** (-log10_e1 / k) if k > 0 else np.nan


# ---------------------------------------------------------------------------
# (c) Constructive non-reversibility


def box_arrays(month: Month, nq):
    lo, hi = np.full(nq, -np.inf), np.full(nq, np.inf)
    b = month.box
    if month.sample == "long":
        lo[:], hi[:] = LP_MARGIN["atm_lo"], b["atm_hi"] - LP_MARGIN["atm_lo"]
        return lo, hi
    for i in range(nq // 3):
        lo[3 * i], hi[3 * i] = LP_MARGIN["atm_lo"], b["atm_hi"] - LP_MARGIN["atm_lo"]
        lo[3 * i + 1], hi[3 * i + 1] = -b["rr_abs"] + LP_MARGIN["rr"], b["rr_abs"] - LP_MARGIN["rr"]
        lo[3 * i + 2], hi[3 * i + 2] = b["bf_lo"] + LP_MARGIN["bf"], b["bf_hi"] - LP_MARGIN["bf"]
    return lo, hi


def lp_direction(N, j, s, X, d, lo, hi, weights):
    """Null-space direction v = N c with v_j = s, |v| <= 1, minimal weighted L1 norm and X + d v in the box."""
    n, k = N.shape
    I = np.eye(n)
    Z = np.zeros((n, n))
    A_ub = [np.hstack([N, -I]), np.hstack([-N, -I])]
    b_ub = [np.zeros(n), np.zeros(n)]
    fh, fl = np.isfinite(hi), np.isfinite(lo)
    if fh.any():
        A_ub.append(np.hstack([d * N[fh], Z[fh]]))
        b_ub.append(hi[fh] - X[fh])
    if fl.any():
        A_ub.append(np.hstack([-d * N[fl], Z[fl]]))
        b_ub.append(X[fl] - lo[fl])
    res = linprog(np.r_[np.zeros(k), weights], A_ub=np.vstack(A_ub), b_ub=np.concatenate(b_ub),
                  A_eq=np.r_[N[j], np.zeros(n)][None, :], b_eq=[s],
                  bounds=[(None, None)] * k + [(0.0, 1.0)] * n, method="highs")
    if res.status != 0:
        return None
    return N @ res.x[:k], float(res.fun)


def project(month, X0, keys, target, half, A0, fixed, x0s, fwd=False):
    """Chord Gauss-Newton back to the level set, coordinate ``fixed`` held; Jacobian refreshed when progress stalls."""
    X = X0.copy()
    free = np.array([i for i in range(len(X)) if i != fixed])
    A, refresh, prev, last_dX = A0, 0, np.inf, None
    for it in range(GN_MAX_IT):
        S, infos = month_sums(month, X, x0s)
        r = np.array([(S[k] - target[k]) / half[k] for k in keys])
        err = float(np.max(np.abs(r)))
        if not np.isfinite(err):
            return X, it, refresh, "evaluation_error"
        x0s = warm_of(infos, x0s)
        if err <= GN_TOL:
            return X, it, refresh, "converged"
        if err > prev and last_dX is not None and refresh >= GN_MAX_REFRESH:
            X[free] -= 0.5 * last_dX  # backtrack once the refreshes are used up
            last_dX = 0.5 * last_dX
            continue
        if err > 0.5 * prev and refresh < GN_MAX_REFRESH:
            A = jacobian(month, X, keys, x0s, fwd=fwd) / half_vec(half, keys)[:, None]
            refresh += 1
        dX = np.linalg.lstsq(A[:, free], -r, rcond=None)[0]
        X[free] += dX
        last_dX, prev = dX, err
    return X, GN_MAX_IT, refresh, "gn_not_converged"


def half_vec(half, keys):
    return np.array([half[k] for k in keys])


def verify(month: Month, X, S_true, pinned_keys, fwd=False):
    """Recompute with the pipeline's warm starts; return (accepted, reasons, diagnostics)."""
    reasons = []
    try:
        S, infos = month_sums(month, X, [{"market": leg.warm[0], "smile": leg.warm[1]} for leg in month.legs], verify=True)
    except (ValueError, ArithmeticError) as exc:  # ArithmeticError: overflow, zero division, floating point
        return False, [f"evaluation_error:{type(exc).__name__}"], {}
    readings = ("market", "smile") if month.sample == "primary" else ("market",) if month.sample == "extended" else ()
    diag = {"cold_max": 0.0, "currency_range_ok": True}
    for i, info in enumerate(infos):
        for r in readings:
            if not (info[f"ok_{r}"] and info[f"resid_{r}"] < CALIB_TOL):
                reasons.append("calibration_not_converged")
            if not info[f"arb_{r}"]:
                reasons.append("arbitrage_check_failed")
            diag["cold_max"] = max(diag["cold_max"], info[f"cold_{r}"])
        if month.sample == "primary" and info.get("moment_status") != "ok":
            reasons.append("moment_status_not_ok")
        atm, rr, bf = (v / VP for v in quotes_of(month, X, i))
        b = month.box
        if not 0 < atm <= b["atm_hi"]:
            reasons.append("implausible_atm")
        if month.sample != "long":
            if abs(rr) > b["rr_abs"]:
                reasons.append("implausible_rr25")
            if not b["bf_lo"] <= bf <= b["bf_hi"]:
                reasons.append("implausible_bf25")
            lo, hi = b["ccy"][month.legs[i].currency]
            diag["currency_range_ok"] &= bool(np.all((np.array([atm, rr, bf]) >= lo) & (np.array([atm, rr, bf]) <= hi)))
    cols = {c: v for c, v in published_columns(month, S).items()
            if not c.endswith(WITHHELD_TAGS)}  # the withheld skew-price series need not be reproduced
    diag["max_resid_half_units"] = max((abs(cols[c] - month.truth[c]) / half_unit(month.truth[c])
                                        for c in cols if half_unit(month.truth[c]) > 0), default=0.0)
    if any(fmt_sig(cols[c]) != fmt_sig(month.truth[c]) for c in cols):
        reasons.append("published_value_mismatch")
    for k in pinned_keys:
        if abs(S[k] - S_true[k]) > GN_TOL * half_unit(S_true[k]):
            reasons.append("theta0_moved")
    if fwd:
        fd = dict(month.others_fd)
        for leg, F in zip(month.legs, forwards_of(month, X)):
            fd[leg.currency] = forward_discount(leg.S, F, G10[leg.currency].usd_base)
        n = month.legs[0].n
        longs, shorts = rank_legs(fd, n)
        if set(longs) != {lg.currency for lg in month.legs if lg.is_long} or set(shorts) != {lg.currency for lg in month.legs if not lg.is_long}:
            reasons.append("ranking_changed")
    diag["pinned_shift"] = {k: (S[k] - S_true[k]) for k in pinned_keys}
    return not reasons, sorted(set(reasons)), diag


def alternatives(month, X_true, keys, S_true, half, A, x_true, coords, targets, fwd=False, stretch=True):
    """Part (c): for each target distance, the first accepted alternative among up to MAX_ATTEMPTS directions."""
    s, rank = svd_stats(A)
    _, _, Vt = np.linalg.svd(A)
    N = Vt[rank:].T
    nq = per_leg(month.sample) * len(month.legs)
    lo, hi = np.full(len(X_true), -np.inf), np.full(len(X_true), np.inf)
    lo[:nq], hi[:nq] = box_arrays(month, nq)
    weights = np.ones(len(X_true))
    if month.sample != "long":
        weights[2:nq:3] = BF_WEIGHT
    out = []
    pinned_keys = [k for k in keys if k.startswith("theta0")]
    todo = list(targets)
    while todo:
        d = todo.pop(0)
        cands = []
        for j in coords:
            for sgn in (1.0, -1.0):
                r = lp_direction(N, j, sgn, X_true, d, lo, hi, weights)
                if r is not None:
                    cands.append((r[1], j, sgn, r[0]))
        cands.sort(key=lambda c: c[0])
        rec = {"target": d, "ok": False, "attempts": 0, "reasons": [], "distance": np.nan}
        if not cands:
            rec["reasons"] = ["no_plausible_null_direction"]
        for _, j, sgn, v in cands[:MAX_ATTEMPTS]:
            rec["attempts"] += 1
            try:
                X, it, refresh, status = project(month, X_true + d * v, keys, S_true, half, A, j, x_true, fwd)
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:  # e.g. SABR nu overflowing far from the truth
                rec["reasons"].append(f"evaluation_error:{type(exc).__name__}")
                continue
            if status != "converged":
                rec["reasons"].append(status)
                continue
            ok, reasons, diag = verify(month, X, S_true, pinned_keys, fwd)
            if not ok:
                rec["reasons"] += reasons
                continue
            dist = np.abs(X - X_true)
            rec.update(ok=True, iterations=it, refreshes=refresh, coord=j % per_leg(month.sample) if j < nq else "fwd",
                       distance=float(dist[:nq].max()) if j < nq else float(dist[nq:].max()),
                       quote_distance=float(dist[:nq].max()), n_moved=int(np.sum(dist[:nq] > 1e-3)),
                       max_resid_half_units=diag["max_resid_half_units"], cold_max=diag["cold_max"],
                       currency_range_ok=diag["currency_range_ok"], pinned_shift=diag["pinned_shift"])
            break
        out.append(rec)
        if stretch and rec["ok"] and d == targets[-1] and STRETCH not in [o["target"] for o in out]:
            todo.append(STRETCH)
    return out


# ---------------------------------------------------------------------------
# Worker


def run_month(month: Month, opts: dict) -> dict:
    global SIG_FIGS
    SIG_FIGS = SIG_BY_SAMPLE[month.sample]  # one month at a time per process
    t0 = time.perf_counter()
    res = {"sample": month.sample, "date": month.date, "n_legs": len(month.legs)}
    try:
        X = np.concatenate([np.array(leg.q[: per_leg(month.sample)]) / VP for leg in month.legs])
        S_true, infos = month_sums(month, X, [{"market": leg.warm[0], "smile": leg.warm[1]} for leg in month.legs])
        x_true = warm_of(infos)
        # (a) validation against the private month-level values
        cols = published_columns(month, S_true)
        gaps = {c: abs(cols[c] - month.truth[c]) / max(abs(month.truth[c]), 1e-300) if month.truth[c] != 0 else abs(cols[c])
                for c in cols}
        gaps.update({f"pinned_{k}": abs(S_true[k] - v) / abs(v) for k, v in month.pinned.items()})
        res["gaps"] = gaps
        res["param_gap"] = max((float(np.max(np.abs(info[f"x_{r}"] - month.flags["panel_x"][i][r])))
                                for i, info in enumerate(infos) for r in ("market", "smile") if f"x_{r}" in info), default=0.0)
        if opts.get("validate_only"):
            res["seconds"] = time.perf_counter() - t0
            return res
        # (b) published independent values
        pub = [k for k in KEYS[month.sample] if k in S_true and month.flags["has"].get(k, False)]
        zero = [k for k in pub if S_true[k] == 0.0]
        pub = [k for k in pub if S_true[k] != 0.0]
        pins = [k for k in PINNED[month.sample] if k in month.pinned]
        keys_c = pub + pins
        half_c = {k: half_unit(S_true[k]) for k in keys_c}
        half = {k: half_c[k] for k in pub}
        Ac = jacobian(month, X, keys_c, x_true) / half_vec(half_c, keys_c)[:, None]
        A = Ac[: len(pub)]
        res["attacker"] = {"J": A * half_vec(half, pub)[:, None], "half": half_vec(half, pub), "X": X,
                           "k": per_leg(month.sample), "sig_p": [leg.sig_p for leg in month.legs]}
        s, rank = svd_stats(A)
        res.update(m=len(pub), m_zero=len(zero), n=len(X), rank=rank, null_dim=len(X) - rank, s_max=float(s[0]),
                   s_min=float(s[rank - 1]) if rank else np.nan, s_ratio=float(s[rank - 1] / s[0]) if rank else np.nan,
                   s_all=";".join(f"{v:.3e}" for v in s))
        grid = np.concatenate([month.flags["grid"]] * len(month.legs))
        res["grid_e1"], res["grid_k"] = grid_stats(A, grid)
        res["grid_gbar"] = float(10 ** np.mean(np.log10(grid)))
        res["drop"] = {k: grid_stats(np.delete(A, i, axis=0), grid) for i, k in enumerate(pub)} if len(pub) > 1 else {}
        if opts["steps"] and month.flags.get("step_check"):
            for h in H_CHECK:
                Ah = jacobian(month, X, pub, x_true, h=h) / half_vec(half, pub)[:, None]
                sh, rank_h = svd_stats(Ah)
                res[f"step_{h:g}_same_rank"] = rank_h == rank
                res[f"step_{h:g}_jac"] = float(np.linalg.norm(Ah - A) / np.linalg.norm(A))
                res[f"step_{h:g}_sv"] = float(np.max(np.abs(sh[:rank] / s[:rank] - 1)))
        # (c) alternatives (every month, or the stratified sample of --construct-sample)
        if opts["construct"] and month.flags.get("construct", True):
            sc, rank_c = svd_stats(Ac)
            res.update(m_c=len(keys_c), rank_c=rank_c, null_dim_c=len(X) - rank_c)
            res["alts"] = alternatives(month, X, keys_c, S_true, half_c, Ac, x_true, range(len(X)), TARGETS)
        # (d) forward points jointly with the quotes
        if opts["construct"] and month.flags.get("forward_check"):
            res["fwd"] = forward_check(month, X, S_true, pub, half, x_true)
    except Exception as exc:  # recorded, not dropped
        res["error"] = f"{type(exc).__name__}: {exc}"
        res["trace"] = traceback.format_exc(limit=3)
    res["seconds"] = time.perf_counter() - t0
    return res


def forward_check(month, Xq, S_true, pub, half, x_true):
    """Part (d): forward points unknown together with the quotes; FD and U added to the equations."""
    nl = len(month.legs)
    X = np.r_[Xq, np.zeros(nl)]
    keys = pub + ["fd", "U"]
    pins = [k for k in PINNED[month.sample] if k in month.pinned]
    keys_c = keys + pins
    hc = {**half, **{k: half_unit(S_true[k]) for k in ["fd", "U"] + pins}}
    Ac = jacobian(month, X, keys_c, x_true, fwd=True) / half_vec(hc, keys_c)[:, None]
    A = Ac[: len(keys)]
    s, rank = svd_stats(A)
    nq = len(Xq)
    sf, rank_f = svd_stats(A[:, nq:])
    alts = alternatives(month, X, keys_c, S_true, hc, Ac, x_true, range(nq, nq + nl), FWD_TARGETS, fwd=True, stretch=False)
    return {"m": len(keys), "n": len(X), "rank": rank, "null_dim": len(X) - rank, "fwd_block_rank": rank_f,
            "fwd_block_sv_ratio": float(sf[-1] / sf[0]), "alts": alts}


# ---------------------------------------------------------------------------
# Building the months (main process)


def grid_of(values):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)] / VP
    for g in GRID_CANDIDATES:
        share = np.mean(np.abs(x / g - np.round(x / g)) < 1e-6)
        if share >= GRID_SHARE:
            return g, float(share)
    return GRID_CANDIDATES[-1], float(np.mean(np.abs(x / GRID_CANDIDATES[-1] - np.round(x / GRID_CANDIDATES[-1])) < 1e-6))


def plausibility(inputs, lo, hi):
    s = inputs[(inputs.date >= lo) & (inputs.date <= hi) & inputs.complete_calib]
    ccy = {}
    for c, g in s.groupby("currency"):
        ccy[c] = (np.array([g.atm.min(), g.rr25.min(), g.bf25.min()]) / VP, np.array([g.atm.max(), g.rr25.max(), g.bf25.max()]) / VP)
    return {"rr_abs": float(s.rr25.abs().max() / VP), "bf_lo": float(s.bf25.min() / VP), "bf_hi": float(s.bf25.max() / VP),
            "atm_lo": float(s.atm.min() / VP), "atm_hi": float(s.atm.max() / VP), "ccy": ccy, "n": len(s)}


def warm_starts(panel):
    """Warm start of every (currency, reading, date) as in calibrate_smiles._calibrate_currency."""
    out, own = {}, {}
    for (c, reading), g in panel.sort_values("date").groupby(["currency", "reading"]):
        last = None
        for f in g.itertuples():
            out[(c, reading, f.date)] = last
            if f.status in ("ok", "not_converged"):
                own[(c, reading, f.date)] = np.array([np.log(f.alpha), np.arctanh(f.rho), np.log(f.nu)])
            if f.status == "ok":
                last = own[(c, reading, f.date)]
    return out, own


def build_months(d: Path, raw: Path):
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    e1, s4, la, mom = load("e1_series.csv"), load("stage4_months.csv"), load("long_atm_months.csv"), load("moments_1m.csv")
    comp_in, comp_pan = load("smile_inputs.csv"), load("smile_panel.csv")
    fn_in, fn_pan = load("smile_inputs_fn.csv"), load("smile_panel_fn.csv")
    spots = {c: daily_spot_mid(raw, c) for c in G10}
    long_in = build_month_end_inputs(raw, ny_month_ends("1995-01-01", "2026-08-31"))
    sel = lambda sample, reading, delta: e1[(e1["sample"] == sample) & (e1.reading == reading) & (e1.delta == delta)
                                            & (e1.status == "ok")].set_index("date")
    s4i = {s: s4[(s4["sample"] == s) & (s4.status == "ok")].set_index("date") for s in ("primary", "extended")}
    la_ok = la[la.status == "ok"].set_index("date")
    prim = {tag: sel("primary", r, dl) for r, dl, tag in COMBOS}
    pp = portfolio_predictors(prim["m10"].reset_index(), mom).set_index("date")
    boxes = {"primary": plausibility(comp_in, PRIMARY_START, pd.Timestamp("2026-08-31")),
             "extended": plausibility(fn_in, pd.Timestamp("2007-01-31"), PRIMARY_START - pd.Timedelta(days=1))}
    lw = long_in[long_in.date < PRIMARY_START]
    boxes["long"] = {"atm_lo": float(lw.atm.min() / VP), "atm_hi": float(lw.atm.max() / VP), "ccy": {}, "n": int(lw.atm.notna().sum())}
    cw = comp_in[(comp_in.date >= PRIMARY_START) & comp_in.complete_calib]
    fw = fn_in[(fn_in.date < PRIMARY_START) & fn_in.complete_calib]
    grids = {"primary": [grid_of(cw[c]) for c in ("atm", "rr25", "bf25")],
             "extended": [grid_of(fw[c]) for c in ("atm", "rr25", "bf25")],
             "long": [grid_of(lw.atm)]}
    warm_c, own_c = warm_starts(comp_pan)
    warm_f, own_f = warm_starts(fn_pan)

    def legs_for(row, inputs, n, warm, own, panel_x):
        idx = inputs.set_index(["currency", "date"])
        legs = []
        for c, is_long in [(c, True) for c in row.longs.split()] + [(c, False) for c in row.shorts.split()]:
            r = idx.loc[(c, row.name)]
            S_end, _ = spot_on(spots[c], option_dates(row.name, c)[2])
            q = (r.atm, r.rr25, r.bf25) if "rr25" in r and np.isfinite(r.rr25) else (r.atm, np.nan, np.nan)
            legs.append(Leg(c, is_long, n, float(r.S), float(r.F), float(r.tau), float(r.df_quote), float(r.df_base),
                            tuple(float(v) for v in q), float(S_end), realised_vol(spots[c], row.name),
                            (warm.get((c, "market", row.name)), warm.get((c, "smile", row.name)))))
            panel_x.append({k: own[(c, k, row.name)] for k in ("market", "smile") if (c, k, row.name) in own})
        return legs

    months = []
    comp_idx = comp_in
    for k, (t, row) in enumerate(prim["m10"].iterrows()):
        panel_x = []
        legs = legs_for(row, comp_idx, 3, warm_c, own_c, panel_x)
        truth, pinned, has = {}, {}, {}
        for _, _, tag in COMBOS:
            truth[f"C_skew_{tag}"], truth[f"phi_{tag}"] = prim[tag].loc[t, "C_skew"], prim[tag].loc[t, "phi"]
            has[f"C_{tag}"] = True
        truth["FD"] = row.FD
        for col in ("var_lo", "var_hi", "oskew_lo", "oskew_hi"):
            truth[col] = pp.loc[t, f"{col}_i"]
            has[col] = True
        if t in s4i["primary"].index:
            r4 = s4i["primary"].loc[t]
            truth.update({c: r4[c] for c in ("U", "H10", "H25", "Hatm", "c1", "c2", "c3")})
            pinned["theta0"] = r4.theta0
            has.update(H25=True, Hatm=True, c1=True, c2=True)
        if t in la_ok.index:  # same legs and quotes: its U and Hatm repeat the primary month's
            rl = la_ok.loc[t]
            if {*rl.longs.split()} != {*row.longs.split()} or {*rl.shorts.split()} != {*row.shorts.split()}:
                raise RuntimeError(f"long at-the-money legs differ from the primary legs in {t:%Y-%m}")
            truth.update({f"long_{c}": rl[c] for c in ("U", "Hatm")})
            pinned["theta0_atm"] = rl.theta0
        months.append(Month("primary", t, legs, row.FD, truth, pinned, boxes["primary"],
                            {"has": has, "panel_x": panel_x, "grid": np.array([g for g, _ in grids["primary"]]),
                             "step_check": k % 8 == 0}))
    ext = sel("extended", "market", 0.10)
    for k, (t, row) in enumerate(ext.iterrows()):
        panel_x = []
        legs = legs_for(row, fn_in, 2, warm_f, own_f, panel_x)
        truth = {"C_skew_m10": row.C_skew, "phi_m10": row.phi, "FD": row.FD}
        pinned, has = {}, {"C_m10": True}
        if t in s4i["extended"].index:
            r4 = s4i["extended"].loc[t]
            truth.update({c: r4[c] for c in ("U", "H10", "H25", "Hatm", "c1", "c2", "c3")})
            pinned["theta0"] = r4.theta0
            has.update(H25=True, Hatm=True, c1=True, c2=True)
        months.append(Month("extended", t, legs, row.FD, truth, pinned, boxes["extended"],
                            {"has": has, "panel_x": panel_x, "grid": np.array([g for g, _ in grids["extended"]]),
                             "step_check": k % 8 == 0}))
    lo_months = la_ok[la_ok.index < PRIMARY_START]
    for k, (t, row) in enumerate(lo_months.iterrows()):
        panel_x = []
        legs = legs_for(row, long_in, int(row.n_legs), {}, {}, panel_x)
        truth = {c: row[c] for c in ("U", "Hatm")}
        months.append(Month("long", t, legs, row.FD, truth, {"theta0_atm": row.theta0}, boxes["long"],
                            {"has": {"Hatm": True}, "panel_x": panel_x,
                             "grid": np.array([grids["long"][0][0]]), "step_check": k % 8 == 0}))
    # part (d) forward check on six primary months with realised returns
    full = [m for m in months if m.sample == "primary" and "U" in m.truth]
    for i in np.linspace(0, len(full) - 1, 6).round().astype(int):
        m = full[i]
        m.flags["forward_check"] = True
        idx = comp_in[(comp_in.date == m.date) & comp_in.complete_calib]
        legs_c = {lg.currency for lg in m.legs}
        m.others_fd = {r.currency: forward_discount(r.S, r.F, G10[r.currency].usd_base) for r in idx.itertuples() if r.currency not in legs_c}
    meta = {"boxes": boxes, "grids": grids, "n_long_overlap": int((la_ok.index >= PRIMARY_START).sum()),
            "n_long_total": len(la_ok), "spots": spots, "stage4": s4i, "e5": load("e5_series.csv")}
    return months, meta


# ---------------------------------------------------------------------------
# (d) spot-only check (main process)


def sigma_fx_check(months, meta, raw, n_months=6):
    """sigma_fx: reverse each currency's daily log changes between anchor dates; every published value is unchanged.

    Runs on ``n_months`` months of each of the primary and extended samples, evenly spaced in date (0 = every month).
    """
    e5 = meta["e5"].set_index("date")
    absr_true = fx_abs_returns(raw)
    logs = {}
    for c in G10:
        sp = meta["spots"][c]
        sp = sp[sp.index.dayofweek < 5]
        logs[c] = np.log(sp)
    todo = []
    for sample in ("primary", "extended"):
        sm = [m for m in months if m.sample == sample and "U" in m.truth and m.date in e5.index]
        n = len(sm) if n_months <= 0 else min(n_months, len(sm))
        todo += [sm[i] for i in sorted(set(np.linspace(0, len(sm) - 1, n).round().astype(int)))]
    out = []
    for m in todo:
        end = option_dates(m.date, "EUR")[2]
        anchors = {m.date, end}
        # every month-end, option expiry and realised-volatility window start near the window
        for t in ny_month_ends(m.date - pd.Timedelta(days=70), end + pd.Timedelta(days=70)):
            anchors.add(t)
            anchors |= {option_dates(t, c)[2] for c in G10}
            for c in G10:
                s = logs[c].loc[:t]
                if len(s) > 21:
                    anchors.add(s.index[-22])
        alt_abs, max_move, n_seg = {}, 0.0, 0
        for c in G10:
            s = logs[c].loc[(logs[c].index >= m.date - pd.Timedelta(days=10)) & (logs[c].index <= end + pd.Timedelta(days=10))]
            dates = s.index
            # a missing anchor date is replaced by the last quote before it, as estimate_stage4.spot_on does
            anchors_c = {dates[dates <= a][-1] for a in anchors if (dates <= a).any()}
            r = np.diff(s.to_numpy())
            cuts = [0] + [k for k in range(1, len(dates) - 1) if dates[k] in anchors_c] + [len(dates) - 1]
            r_alt = r.copy()
            for a, b in zip(cuts[:-1], cuts[1:]):
                if dates[a] >= m.date and dates[b] <= end:  # segments inside this month's window only
                    r_alt[a:b] = r[a:b][::-1]
                    n_seg += 1
            path = s.iloc[0] + np.r_[0.0, np.cumsum(r_alt)]
            max_move = max(max_move, float(np.max(np.abs(np.exp(path - s.to_numpy()) - 1))))
            alt_abs[c] = pd.Series(np.abs(r_alt), index=dates[1:])
        frame = pd.DataFrame(alt_abs)
        win = lambda x: x.loc[(x.index > m.date) & (x.index <= end)]
        alt = float(win(frame).mean(axis=1, skipna=True).mean())
        true = float(win(absr_true).mean())
        out.append({"sample": m.sample, "same_6sf": fmt_sig(alt) == fmt_sig(e5.loc[m.date, "sigma_fx"]),
                    "rel_gap": abs(alt / true - 1), "true_gap": abs(true / e5.loc[m.date, "sigma_fx"] - 1),
                    "max_spot_move": max_move, "segments": n_seg})
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Summary


def q(x, ps=(0.0, 0.05, 0.5, 0.95, 1.0), f="{:.3g}"):
    x = pd.Series(x, dtype=float).dropna()
    if x.empty:
        return "n/a"
    return ", ".join(f.format(x.quantile(p)) for p in ps)


def summarise(results, meta, sfx, runtime, opts):
    L = ["# Reversibility of the public portfolio-level series (restricted; pooled statistics only)", "",
         "Produced by scripts/check_public_reversibility.py. Adversary: knows everything about a month except the "
         "smile quotes (leg identities and orientations, spot, forward points, rates, dates, spot at expiry, realised "
         "volatility forecast). Unknowns: ATM, RR25 and BF25 of each leg (ATM only in the long at-the-money sample). "
         f"Published precision: six significant figures in primary.csv and long_atm.csv, four in extended.csv. Final publication set: primary.csv without rr_or, "
         "extended.csv with sigma_fx, long_atm.csv with U and Hatm only (payoff and premium withdrawn).", "",
         f"Run time {runtime / 60:.1f} minutes on {opts['workers']} processes. Months processed: "
         + ", ".join(f"{s} {sum(r['sample'] == s for r in results)}" for s in ("primary", "extended", "long"))
         + f". Long at-the-money months from the primary start ({meta['n_long_overlap']} of {meta['n_long_total']}) "
         "have the same legs and quotes as primary months (checked) and are analysed jointly with them: their U and "
         "Hatm repeat the primary month's, so they add no equation but must round to the same published values; the "
         "remaining long at-the-money months are analysed alone.", ""]
    if opts.get("validate_only"):
        L[-2] = L[-2].replace("Run time", "Validation of the forward map only (part a). Run time")
    errs = [r for r in results if "error" in r]
    L += [f"Months with an error in the check itself: {len(errs)}"
          + ("" if not errs else " (" + "; ".join(sorted({r['error'][:80] for r in errs})) + ")"), ""]
    ok = [r for r in results if "error" not in r]
    by = lambda s: [r for r in ok if r["sample"] == s]

    L += ["## (a) Forward map: reproduction of the private month-level values at the true quotes", "",
          "| sample | months | columns checked | largest relative gap | column with the largest gap | largest gap of recalibrated SABR parameters |",
          "| --- | --- | --- | --- | --- | --- |"]
    for s in ("primary", "extended", "long"):
        R = by(s)
        if not R:
            continue
        g = defaultdict(float)
        for r in R:
            for c, v in r["gaps"].items():
                g[c] = max(g[c], v)
        worst = max(g, key=g.get)
        L.append(f"| {s} | {len(R)} | {len(g)} | {g[worst]:.2e} | {worst} | {max(r['param_gap'] for r in R):.2e} |")
    pooled = defaultdict(float)
    for r in ok:
        for c, v in r["gaps"].items():
            pooled[c] = max(pooled[c], v)
    L += ["", "Largest relative gap by column, pooled over samples (absolute where the true value is zero): "
          + "; ".join(f"{c} {v:.1e}" for c, v in sorted(pooled.items())), ""]
    if opts.get("validate_only"):
        return "\n".join(L)

    L += ["## (b) Counting and local rank", "",
          "Jacobian rows: independent published smile-dependent values, scaled to half units of the sixth significant "
          "figure; columns: quotes, per volatility point. Exact identities (c3, H10, phi, and the long_atm.csv U and "
          "Hatm of months from the primary start) are excluded. A row whose true value is exactly zero would be "
          f"excluded from the rank; rank threshold {RANK_RTOL:g} of the largest singular value.", "",
          "| sample | months | (values m, unknowns n): months | rank = m in | null-space dimension: months | s_min/s_max (min, 5%, median, 95%, max) | s_max (median) |",
          "| --- | --- | --- | --- | --- | --- | --- |"]
    for s in ("primary", "extended", "long"):
        R = by(s)
        if not R:
            continue
        mn = Counter((r["m"], r["n"]) for r in R)
        nd = Counter(r["null_dim"] for r in R)
        L.append(f"| {s} | {len(R)} | " + ", ".join(f"({a}, {b}): {v}" for (a, b), v in sorted(mn.items()))
                 + f" | {sum(r['rank'] == r['m'] for r in R)}/{len(R)} | " + ", ".join(f"{k}: {v}" for k, v in sorted(nd.items()))
                 + f" | {q([r['s_ratio'] for r in R], f='{:.2e}')} | {np.median([r['s_max'] for r in R]):.3g} |")
    zero = sum(r["m_zero"] for r in ok)
    L += ["", f"Rows excluded because the true value is exactly zero: {zero}.", ""]
    st = [r for r in ok if "step_0.0001_jac" in r]
    if st:
        L += [f"Step-size stability on {len(st)} months (every eighth month of each sample): relative Frobenius "
              f"difference of the scaled Jacobian from the step 1e-5 one, largest {max(r['step_0.0001_jac'] for r in st):.2e} "
              f"(step 1e-4) and {max(r['step_1e-06_jac'] for r in st):.2e} (step 1e-6); largest relative change of a "
              f"nonzero singular value {max(r['step_0.0001_sv'] for r in st):.2e} and {max(r['step_1e-06_sv'] for r in st):.2e}; "
              f"rank unchanged in {sum(r['step_0.0001_same_rank'] and r['step_1e-06_same_rank'] for r in st)}/{len(st)} "
              "of them.", ""]
    L += ["### Grid-aware adversary (linearised)", "",
          "Quote grid measured from the data (coarsest step covering at least "
          f"{GRID_SHARE:.1%} of the quotes; volatility points): "
          + "; ".join(f"{s}: " + ", ".join(f"{n} {g:g} ({sh:.1%})" for n, (g, sh) in zip(('ATM', 'RR25', 'BF25'), meta['grids'][s]))
                      for s in ("primary", "extended", "long")) + ".", "",
          "Linearised count of grid points consistent with the rounded published values: E1 per unit null-space volume "
          "(volatility points^k), so a null-space box of side L around the truth holds E1 L^k of them, and consistent grid "
          "points lie typically L* = E1^(-1/k) volatility points apart along the level set. An adversary who knows the grid "
          "can single out the truth only with a prior that confines the quotes to within about L* in every null-space "
          "direction, and must then still search about (L*/g)^k grid combinations of k free quotes (g: geometric mean grid "
          "step). Statistics are (min, 5%, median, 95%, max) over months.", "",
          "| sample | k | L*, vol. points | log10 E1 L^k for L = 1, 2, 4 (medians) | months with fewer than 1 point in a box of side 2 | log10 (L*/g)^k | median L* at 5, 4, 3 s.f. |",
          "| --- | --- | --- | --- | --- | --- | --- |"]
    for s in ("primary", "extended", "long"):
        R = by(s)
        if not R:
            continue
        e1 = np.array([r["grid_e1"] for r in R])
        k = np.array([r["grid_k"] for r in R])
        mm = np.array([r["rank"] for r in R])
        ls = np.array([spacing(a, b) for a, b in zip(e1, k)])
        L.append(f"| {s} | " + ", ".join(f"{a}: {b}" for a, b in sorted(Counter(k).items())) + f" | {q(ls)} | "
                 + ", ".join(f"{np.median(e1 + k * np.log10(side)):.1f}" for side in (1, 2, 4))
                 + f" | {np.mean(e1 + k * np.log10(2.0) < 0):.1%} | {q(k * np.log10(ls / np.array([r['grid_gbar'] for r in R])), f='{:.1f}')} | "
                 + ", ".join(f"{np.median([spacing(a + b * (SIG_BY_SAMPLE[s] - p), c) for a, b, c in zip(e1, mm, k)]):.3g}" for p in (5, 4, 3)) + " |")
    L += ["", "Fewer significant figures widen every interval tenfold per figure removed, so log10 E1 rises by m per figure. "
          "Effect of dropping one independent value (months with every value of their sample; median L* after dropping it, "
          "and the published columns that would have to go with it because of the identities):", ""]
    for s in ("primary", "extended"):
        full = len(KEYS[s])
        P = [r for r in by(s) if r["m"] == full]
        if not P:
            continue
        base = np.median([spacing(r["grid_e1"], r["grid_k"]) for r in P])
        L += [f"{s} ({len(P)} months; median L* with every value {base:.3g}):", "",
              "| value dropped | median L* after dropping | months with at least 1 point in a box of side 2 after dropping | columns to remove |",
              "| --- | --- | --- | --- |"]
        for key in KEYS[s]:
            D = [r["drop"][key] for r in P if key in r["drop"]]
            if D:
                L.append(f"| {key} | {np.median([spacing(a, b) for a, b in D]):.3g} | "
                         f"{np.mean([a + b * np.log10(2.0) >= 0 for a, b in D]):.1%} | "
                         f"{(DROP_COLUMNS_EXT if s == 'extended' else DROP_COLUMNS_PRIM).get(key) or DROP_COLUMNS[key]} |")
        L += [""]

    if opts["construct"]:
        L += ["## (c) Constructive non-reversibility", "",
              "Each alternative reproduces every published column of its month to its file's significant figures (six; four in extended.csv) (recomputed "
              "with the pipeline's own warm starts), keeps the month's theta0 values within 0.01 half units, "
              "calibrates under every reading used (largest residual below 1e-8), passes both arbitrage checks, has "
              "moment status ok (primary), ATM > 0 and no larger than the sample's largest ATM, |RR25| and BF25 within "
              "the sample's range. Distance: largest "
              "absolute difference from the true quotes, volatility points. The 4-point target is tried only after "
              "2 points succeed. "
              + ("Run on every month." if not opts.get("construct_sample") else
                 f"Run on a stratified sample of {opts['construct_sample']} months per sample (evenly spaced in formation "
                 "date, first and last month included; every month of a sample with fewer months), because the full "
                 "run would take too long; parts (a) and (b) above cover every month."), "",
              "| sample | months | equations incl. pinned theta0: months | null dim. (with theta0): months | >= 0.5 | >= 1 | >= 2 | >= 4 | largest distance (min, 5%, median, 95%, max) |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        for s in ("primary", "extended", "long"):
            R = [r for r in by(s) if "alts" in r]
            if not R:
                continue
            share = lambda d: np.mean([any(a["ok"] and a["target"] >= d for a in r["alts"]) for r in R])
            best = [max([a["distance"] for a in r["alts"] if a["ok"]], default=0.0) for r in R]
            L.append(f"| {s} | {len(R)} | " + ", ".join(f"{k}: {v}" for k, v in sorted(Counter(r['m_c'] for r in R).items()))
                     + " | " + ", ".join(f"{k}: {v}" for k, v in sorted(Counter(r['null_dim_c'] for r in R).items()))
                     + f" | {share(0.5):.1%} | {share(1.0):.1%} | {share(2.0):.1%} | {share(STRETCH):.1%} | {q(best, f='{:.2f}')} |")
        L += ["", "Accepted alternatives, pooled:", ""]
        A = [a for r in ok if "alts" in r for a in r["alts"] if a["ok"]]
        if A:
            coord = Counter({0: "ATM", 1: "RR25", 2: "BF25"}.get(a["coord"], a["coord"]) for a in A)
            L += [f"- {len(A)} alternatives; quote moved by the full distance: " + ", ".join(f"{k} {v}" for k, v in coord.items())
                  + f"; quotes moved by more than 0.001 volatility points: median {np.median([a['n_moved'] for a in A]):.0f}, "
                  f"min {min(a['n_moved'] for a in A)}, max {max(a['n_moved'] for a in A)}",
                  f"- largest residual over all published columns: {max(a['max_resid_half_units'] for a in A):.2e} half units "
                  f"of the sixth significant figure; Gauss-Newton iterations median {np.median([a['iterations'] for a in A]):.0f} "
                  f"(max {max(a['iterations'] for a in A)}), Jacobian refreshes max {max(a['refreshes'] for a in A)}",
                  f"- largest difference between the warm-started and a cold-started calibration of the alternative quotes: "
                  f"{max(a['cold_max'] for a in A):.2e} (relative for alpha and nu, absolute for rho)",
                  f"- alternatives whose quotes also lie inside each currency's own observed range: "
                  f"{np.mean([a['currency_range_ok'] for a in A if a['coord'] != 'fwd']):.1%} (primary and extended only: "
                  f"{np.mean([a['currency_range_ok'] for r in ok if r['sample'] != 'long' and 'alts' in r for a in r['alts'] if a['ok']]):.1%})"]
        L += ["", "Failures (months without an accepted alternative at the target), by target and reason (a month can have several reasons):", ""]
        for s in ("primary", "extended", "long"):
            R = [r for r in by(s) if "alts" in r]
            for d in (*TARGETS, STRETCH):
                fails = [a for r in R for a in r["alts"] if a["target"] == d and not a["ok"]]
                if not fails:
                    continue
                reasons = Counter(x.split(":")[0] for a in fails for x in set(a["reasons"]))
                years = Counter(r["date"].year for r in R for a in r["alts"] if a["target"] == d and not a["ok"])
                L.append(f"- {s}, {d:g} points: {len(fails)} months; reasons " + ", ".join(f"{k} {v}" for k, v in reasons.most_common())
                         + "; by formation year " + ", ".join(f"{y}: {v}" for y, v in sorted(years.items())))
        L += [""]

    L += ["## (d) Quotes that enter linearly, or only through forwards and spot", ""]
    L += ["- rr_or (a plain signed average of the six quoted 10-delta risk reversals, in quote units) is withdrawn from "
          "primary.csv as the least transformed column and is not analysed."]
    F = [r["fwd"] for r in ok if "fwd" in r]
    if F:
        L += [f"- FD and U: two equations in the six forward points when spot is known. Forward points also enter every "
              f"smile-dependent value, so I solved the joint problem (quotes and forward points unknown) on {len(F)} primary "
              f"months: equations {Counter(f['m'] for f in F)}, unknowns {Counter(f['n'] for f in F)}, rank "
              f"{Counter(f['rank'] for f in F)}, null-space dimension {Counter(f['null_dim'] for f in F)}. Alternatives that "
              "keep every published value (FD, phi and U included), the pinned theta0, the leg ranking against the other "
              "currencies and the plausibility and calibration conditions, with the forwards moved by at least: "
              + ", ".join(f"{d:g} bp of log F in {sum(any(a['ok'] and a['target'] >= d for a in f['alts']) for f in F)}/{len(F)} months"
                          for d in FWD_TARGETS)
              + ". Failure reasons: " + (", ".join(f"{k} {v}" for k, v in Counter(x.split(':')[0] for f in F for a in f['alts'] if not a['ok'] for x in set(a['reasons'])).most_common()) or "none")
              + f". For an adversary who already held the month's smile quotes, the forward block alone has rank "
              f"{Counter(f['fwd_block_rank'] for f in F)} out of 6 (smallest to largest singular value "
              f"{q([f['fwd_block_sv_ratio'] for f in F], ps=(0.0, 0.5, 1.0), f='{:.1e}')}), so the published values "
              "would pin down the forward points locally; that adversary already holds LSEG quotes for the month."]
    if sfx is not None and len(sfx):
        L += ["- sigma_fx: one number per month from about 9 x 21 absolute daily log spot changes. Reversing each "
              "currency's daily changes between consecutive dates whose spot enters any other published value (month-ends, "
              "option expiries, starts of the realised-volatility windows) keeps those spots, every window's set of "
              "changes and hence every published value. "
              + "; ".join(f"{smp}, {len(g)} months: sigma_fx identical at six significant figures in "
                          f"{int(g.same_6sf.sum())}/{len(g)} (largest relative gap {g.rel_gap.max():.1e}; my sigma_fx at the true "
                          f"spots against the private value: largest relative gap {g.true_gap.max():.1e}), segments reversed per "
                          f"month (all currencies) {q(g.segments, ps=(0.0, 0.5, 1.0), f='{:.0f}')}, largest move of an intermediate "
                          f"daily spot inside the window {q(g.max_spot_move, ps=(0.0, 0.5, 1.0), f='{:.2%}')} (min, median, max)"
                          for smp, g in sfx.groupby("sample", sort=False)) + "."]
    L += ["", "## (e) Formal argument", "",
          "A continuous map from an open subset of R^n into R^m with m < n is not injective (by invariance of domain: "
          "composed with the inclusion of R^m in R^n it would be a continuous injection whose image, inside R^m x {0}, "
          "has empty interior, yet must be open). Every month above has m < n, so no adversary can recover the quotes "
          "uniquely from the exact values, let alone from values rounded to six significant figures; part (c) shows the "
          "non-uniqueness is material. The argument assumes continuous quotes; the grid-aware count in (b) is the "
          "relevant measure for an adversary who also knows the quote grid.", ""]
    return "\n".join(L)


# ---------------------------------------------------------------------------


def linear_attacker(results) -> str:
    """(f) The best linear attacker: posterior standard deviation of each quote given the month's published values.

    The attacker knows everything except the smile quotes and is given the exact Jacobian J of the published
    values at the true quotes, which no real attacker has. Prior: independent normal errors per quote type with
    the standard deviation of a public-information predictor over the sample (ATM: the realised-volatility
    forecast; RR25: zero; BF25: the pooled median). Posterior covariance S0 - S0 J'(J S0 J' + E)^-1 J S0, with E
    the variance of the rounding to six significant figures. Pooled output only.
    """
    lines = ["## (f) The best linear attacker", "",
             "Given the exact Jacobian at the true quotes (which presupposes knowing them), the legs, spot and forwards, "
             "and a normal prior per quote type with the spread of a public-information predictor, the attacker's "
             "posterior standard deviation per quote, in volatility points:", "",
             "| sample | months | published values per month | unknown quotes | quote | prior sd | posterior sd, median | "
             "10th percentile | share below 0.10 | share below 0.25 |", "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for smp in ("primary", "extended", "long"):
        rs = [r["attacker"] for r in results if r.get("sample") == smp and "attacker" in r]
        if not rs:
            continue
        k = rs[0]["k"]
        bf_med = float(np.median(np.concatenate([a["X"][2::3] for a in rs]))) if k == 3 else None
        err = {j: [] for j in range(k)}
        for a in rs:
            for i, sp in enumerate(a["sig_p"]):
                prior = [sp / VP if np.isfinite(sp) else np.nan, 0.0, bf_med][:k]
                for j in range(k):
                    if np.isfinite(prior[j]):
                        err[j].append(a["X"][i * k + j] - prior[j])
        sd0 = np.array([np.sqrt(np.mean(np.square(err[j]))) for j in range(k)])
        post = {j: [] for j in range(k)}
        for a in rs:
            n = len(a["X"])
            S0 = np.diag(np.tile(sd0 ** 2, n // k))
            J, E = a["J"], np.diag(a["half"] ** 2 / 3.0)
            S1 = S0 - S0 @ J.T @ np.linalg.solve(J @ S0 @ J.T + E, J @ S0)
            sd1 = np.sqrt(np.clip(np.diag(S1), 0.0, None))
            for i in range(n // k):
                for j in range(k):
                    post[j].append(sd1[i * k + j])
        m_eq = int(np.median([len(a["half"]) for a in rs]))
        for j, name in enumerate(("ATM", "RR25", "BF25")[:k]):
            q = np.array(post[j])
            lines.append(f"| {smp} | {len(rs)} | {m_eq} | {len(rs[0]['X'])} | {name} | {sd0[j]:.2f} | {np.median(q):.2f} | "
                         f"{np.quantile(q, 0.1):.2f} | {np.mean(q < 0.10):.1%} | {np.mean(q < 0.25):.1%} |")
    return "\n".join(lines) + "\n"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--limit", type=int, default=0, help="months per sample, for testing (0 = all)")
    p.add_argument("--no-construct", action="store_true", help="skip parts (c) and the forward check of (d)")
    p.add_argument("--construct-sample", type=int, default=0,
                   help="run part (c) on this many months per sample, evenly spaced in date (0 = every month)")
    p.add_argument("--validate-only", action="store_true", help="run part (a) alone on every month and print it")
    p.add_argument("--sigma-fx-months", type=int, default=0,
                   help="months per sample in the sigma_fx check, evenly spaced (0 = every month)")
    p.add_argument("--out-dir", default="", help="directory for the outputs (default: the private results folder)")
    args = p.parse_args()
    if args.workers > 8:
        raise SystemExit("at most 8 processes")
    t0 = time.perf_counter()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    raw = ROOT / "data" / "private" / "lseg" / args.retrieval_date / "raw"
    out_dir = Path(args.out_dir) if args.out_dir else d
    if (ROOT / "data" / "public") in [out_dir.resolve(), *out_dir.resolve().parents]:
        raise SystemExit("this script never writes into data/public")
    months, meta = build_months(d, raw)
    if args.limit:
        keep = []
        for s in ("primary", "extended", "long"):
            sm = [m for m in months if m.sample == s]
            keep += [sm[i] for i in np.linspace(0, len(sm) - 1, min(args.limit, len(sm))).round().astype(int)]
        months = keep
    if args.construct_sample:
        for s in ("primary", "extended", "long"):
            sm = [m for m in months if m.sample == s]
            chosen = set(np.linspace(0, len(sm) - 1, min(args.construct_sample, len(sm))).round().astype(int))
            for i, m in enumerate(sm):
                m.flags["construct"] = i in chosen
    opts = {"steps": not args.validate_only, "construct": not (args.no_construct or args.validate_only),
            "construct_sample": args.construct_sample, "validate_only": args.validate_only, "workers": args.workers}
    print(f"{len(months)} months; built in {time.perf_counter() - t0:.0f} s", file=sys.stderr, flush=True)
    weight = {"primary": 0, "extended": 1, "long": 2}
    order = sorted(months, key=lambda m: weight[m.sample])
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(run_month, m, opts) for m in order]
        for i, f in enumerate(as_completed(futs), 1):
            results.append(f.result())
            if i % 25 == 0 or i == len(futs):
                print(f"{i}/{len(futs)} months done, {time.perf_counter() - t0:.0f} s", file=sys.stderr, flush=True)
    results.sort(key=lambda r: (weight[r["sample"]], r["date"]))
    if args.validate_only:  # printed only; the outputs of a full run are left as they are
        print(summarise(results, meta, None, time.perf_counter() - t0, opts))
        return
    att = linear_attacker(results)
    sfx = sigma_fx_check(months, meta, raw, args.sigma_fx_months)
    print(f"sigma_fx check done, {time.perf_counter() - t0:.0f} s", file=sys.stderr, flush=True)
    text = summarise(results, meta, sfx, time.perf_counter() - t0, opts) + "\n" + att
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "public_reversibility.md").write_text(text)
    rows = []
    for r in results:
        row = {k: r.get(k) for k in ("sample", "date", "n_legs", "m", "m_zero", "n", "rank", "null_dim", "s_min", "s_max",
                                     "s_ratio", "s_all", "grid_e1", "grid_k", "m_c", "rank_c", "null_dim_c", "seconds", "error")}
        row["construct"] = "alts" in r
        row["largest_gap_a"] = max(r["gaps"].values()) if "gaps" in r else np.nan
        for a in r.get("alts", []):
            row[f"ok_{a['target']:g}"] = a["ok"]
            row[f"distance_{a['target']:g}"] = a["distance"]
            row[f"reasons_{a['target']:g}"] = ";".join(sorted(set(a["reasons"]))) if not a["ok"] else ""
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_dir / "public_reversibility_months.csv", index=False)
    print(text)


if __name__ == "__main__":
    main()
