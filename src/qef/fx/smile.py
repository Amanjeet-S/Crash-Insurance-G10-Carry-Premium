"""Smile construction from ATM, risk-reversal and butterfly quotes (result R8).

Quotes are in volatility units (0.10 = 10%). The risk reversal is σ(call) −
σ(put) at the stated delta on the quoted pair. Two butterfly readings are
supported:

``"market"``  market strangle: the single volatility σ_ATM + BF, at the strikes
              that volatility implies, prices a strangle whose value the smile
              must reproduce at those same strikes (Reiswich and
              Wystup, 2012).
``"smile"``   smile strangle: [σ(K_25C) + σ(K_25P)]/2 − σ_ATM = BF at the
              smile's own delta strikes.

The system of three equations in (α, ρ, ν) at fixed β is solved by
Levenberg–Marquardt in the unconstrained variables (ln α, atanh ρ, ln ν).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq, least_squares

from .gk import CALL, PUT, DeltaConvention, atm_dns_strike, forward_premium, strike_from_delta, strike_from_delta_smile, vega
from .sabr import sabr_vol


@dataclass(frozen=True)
class SmileQuotes:
    F: float
    tau: float
    sigma_atm: float
    rr: float
    bf: float
    delta: float = 0.25
    df_base: float = 1.0  # base-currency discount factor, needed for spot delta
    convention: DeltaConvention = field(default_factory=DeltaConvention)


@dataclass(frozen=True)
class SabrSmile:
    F: float
    tau: float
    alpha: float
    rho: float
    nu: float
    beta: float = 1.0

    def vol(self, K):
        return sabr_vol(K, self.F, self.tau, self.alpha, self.rho, self.nu, self.beta)


@dataclass(frozen=True)
class CalibrationResult:
    smile: SabrSmile
    residuals: np.ndarray  # (ATM, RR, strangle) in volatility units
    success: bool
    message: str
    strikes: dict  # K_atm, K_call, K_put under the smile; K_call_ms, K_put_ms for the market strangle
    jacobian_cond: float


def market_strangle_strikes(q: SmileQuotes):
    """Strikes of the market strangle: the ±Δ strikes at the flat volatility σ_ATM + BF."""
    s = q.sigma_atm + q.bf
    kc = strike_from_delta(q.delta, q.F, s, q.tau, CALL, q.convention, q.df_base)
    kp = strike_from_delta(-q.delta, q.F, s, q.tau, PUT, q.convention, q.df_base)
    return kc, kp


def _smile_delta_strikes(smile: SabrSmile, q: SmileQuotes):
    kc = strike_from_delta_smile(q.delta, q.F, q.tau, CALL, q.convention, smile.vol, q.df_base)
    kp = strike_from_delta_smile(-q.delta, q.F, q.tau, PUT, q.convention, smile.vol, q.df_base)
    return kc, kp


def _atm_strike_smile(smile: SabrSmile, q: SmileQuotes):
    """ATM DNS strike under the smile: fixed point K = K_DNS(σ(K))."""
    g = lambda lk: lk - np.log(atm_dns_strike(q.F, float(smile.vol(np.exp(lk))), q.tau, q.convention))
    lk0 = np.log(atm_dns_strike(q.F, q.sigma_atm, q.tau, q.convention))
    width = 0.5 * q.sigma_atm**2 * q.tau + 1e-6
    return float(np.exp(brentq(g, lk0 - 20 * width, lk0 + 20 * width, xtol=1e-15)))


RHO_MAX = 0.999


def _unpack(x, F, tau, beta):
    # |rho| is kept below 1: at rho = ±1 the Hagan expansion divides by zero.
    rho = float(np.clip(np.tanh(x[1]), -RHO_MAX, RHO_MAX))
    return SabrSmile(F, tau, float(np.exp(x[0])), rho, float(np.exp(x[2])), beta)


def _residuals(x, q: SmileQuotes, bf_reading: str, beta: float, k_ms):
    smile = _unpack(x, q.F, q.tau, beta)
    try:
        with np.errstate(all="ignore"):
            k_atm = _atm_strike_smile(smile, q)
            kc, kp = _smile_delta_strikes(smile, q)
    except (ValueError, FloatingPointError):
        return np.full(3, 1.0)  # outside the region where the fixed points exist
    if not np.all(np.isfinite([k_atm, kc, kp])):
        return np.full(3, 1.0)
    r_atm = float(smile.vol(k_atm)) - q.sigma_atm
    vc, vp = float(smile.vol(kc)), float(smile.vol(kp))
    r_rr = (vc - vp) - q.rr
    if bf_reading == "smile":
        r_bf = 0.5 * (vc + vp) - q.sigma_atm - q.bf
    else:
        kcm, kpm = k_ms
        s = q.sigma_atm + q.bf
        target = forward_premium(q.F, kcm, s, q.tau, CALL) + forward_premium(q.F, kpm, s, q.tau, PUT)
        value = forward_premium(q.F, kcm, smile.vol(kcm), q.tau, CALL) + forward_premium(q.F, kpm, smile.vol(kpm), q.tau, PUT)
        strangle_vega = vega(q.F, kcm, s, q.tau, 1.0) + vega(q.F, kpm, s, q.tau, 1.0)
        r_bf = float(value - target) / float(strangle_vega)  # volatility-equivalent units
    out = np.array([r_atm, r_rr, r_bf])
    # Non-finite residuals can stall MINPACK; treat them as an infeasible point.
    return out if np.all(np.isfinite(out)) else np.full(3, 1.0)


def calibrate_sabr(q: SmileQuotes, bf_reading: str = "market", beta: float = 1.0, x0=None, tol=1e-12, max_nfev=600) -> CalibrationResult:
    """Calibrate a SABR smile to ATM, RR and BF quotes (result R8)."""
    if bf_reading not in ("market", "smile"):
        raise ValueError("bf_reading must be 'market' or 'smile'")
    k_ms = market_strangle_strikes(q) if bf_reading == "market" else None
    starts = [x0] if x0 is not None and np.all(np.isfinite(x0)) else []
    rho_guess = np.clip(q.rr / max(4.0 * q.bf + 1e-4, 1e-3), -0.9, 0.9)
    for nu in (0.5, 1.5, 3.0):
        starts.append(np.array([np.log(q.sigma_atm), np.arctanh(rho_guess), np.log(nu)]))
    best = None
    for s in starts:
        sol = least_squares(_residuals, s, args=(q, bf_reading, beta, k_ms), method="lm", xtol=tol, ftol=tol, gtol=tol, max_nfev=max_nfev)
        if best is None or sol.cost < best.cost:
            best = sol
        if sol.cost < 1e-24:
            break
    smile = _unpack(best.x, q.F, q.tau, beta)
    k_atm = _atm_strike_smile(smile, q)
    kc, kp = _smile_delta_strikes(smile, q)
    strikes = {"K_atm": k_atm, "K_call": kc, "K_put": kp}
    if k_ms is not None:
        strikes.update({"K_call_ms": k_ms[0], "K_put_ms": k_ms[1]})
    if np.all(np.isfinite(best.jac)):
        sv = np.linalg.svd(best.jac, compute_uv=False)
        cond = float(sv[0] / sv[-1]) if sv[-1] > 0 else np.inf
    else:
        cond = np.nan
    res = _residuals(best.x, q, bf_reading, beta, k_ms)
    ok = bool(np.max(np.abs(res)) < 1e-8)
    return CalibrationResult(smile, res, ok, best.message, strikes, cond)


def quotes_from_smile(smile: SabrSmile, tau, convention: DeltaConvention, delta=0.25, df_base=1.0, bf_reading="market"):
    """Inverse map: the ATM, RR and BF quotes a given smile implies.

    Used for synthetic recovery tests and for E4 (the effect of misreading
    the butterfly convention).
    """
    F = smile.F
    probe = SmileQuotes(F, tau, float(smile.vol(F)), 0.0, 0.0, delta, df_base, convention)
    k_atm = _atm_strike_smile(smile, probe)
    s_atm = float(smile.vol(k_atm))
    q0 = SmileQuotes(F, tau, s_atm, 0.0, 0.0, delta, df_base, convention)
    kc, kp = _smile_delta_strikes(smile, q0)
    vc, vp = float(smile.vol(kc)), float(smile.vol(kp))
    rr = vc - vp
    if bf_reading == "smile":
        return s_atm, rr, 0.5 * (vc + vp) - s_atm

    def gap(bf):
        qq = SmileQuotes(F, tau, s_atm, rr, bf, delta, df_base, convention)
        kcm, kpm = market_strangle_strikes(qq)
        s = s_atm + bf
        flat = forward_premium(F, kcm, s, tau, CALL) + forward_premium(F, kpm, s, tau, PUT)
        smiley = forward_premium(F, kcm, smile.vol(kcm), tau, CALL) + forward_premium(F, kpm, smile.vol(kpm), tau, PUT)
        return float(flat - smiley)

    # Bracket the market-strangle butterfly around the smile-strangle value,
    # keeping the strangle volatility s_atm + bf positive.
    bf_ss = 0.5 * (vc + vp) - s_atm
    lo = max(bf_ss - 0.05, -0.5 * s_atm)
    hi = bf_ss + 0.05
    return s_atm, rr, brentq(gap, lo, hi, xtol=1e-14)
