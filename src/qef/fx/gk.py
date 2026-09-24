"""Garman–Kohlhagen pricing, delta conventions and delta-to-strike inversion.

Notation follows the research design. A pair is quoted as units of the quote
currency per unit of the base currency (EURUSD: USD per EUR; USDJPY: JPY per
USD). Premiums are in quote currency per unit of base-currency notional.

Discount factors: ``df_quote`` discounts quote-currency cash flows to the
valuation date. The base-currency discount factor is implied by the forward,
``df_base = F * df_quote / S``, so no separate covered-interest-parity
assumption enters (research design, section 3.3).

Sources: the pricing formula is Garman and Kohlhagen (1983). The spot,
forward and premium-adjusted delta definitions and the delta-neutral-straddle
ATM strike follow Reiswich and Wystup (2012). Full
references are in projects/01_fx_crash_insurance_carry/references.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

CALL = 1
PUT = -1
MAX_EXPANSIONS = 200  # bound on bracket expansions; a valid bracket needs far fewer


def _check_vol(sigma, tau):
    if not (np.isfinite(sigma) and np.isfinite(tau) and sigma > 0 and tau > 0):
        raise ValueError("volatility and time to expiry must be positive and finite")


@dataclass(frozen=True)
class DeltaConvention:
    """Delta quotation convention for one currency pair.

    ``spot``: spot delta if True, forward delta otherwise.
    ``premium_adjusted``: premium-adjusted delta if True (premium paid in the
    base currency), pips delta otherwise.
    """

    spot: bool = True
    premium_adjusted: bool = False


def d_plus_minus(F, K, sigma, tau):
    """Return (d+, d-) with d± = [ln(F/K) ± σ²τ/2] / (σ√τ)."""
    F, K, sigma, tau = np.broadcast_arrays(*map(np.asarray, (F, K, sigma, tau)))
    sd = sigma * np.sqrt(tau)
    d1 = (np.log(F / K) + 0.5 * sd**2) / sd
    return d1, d1 - sd


def forward_premium(F, K, sigma, tau, phi):
    """Undiscounted Garman–Kohlhagen premium φ[F Φ(φd+) − K Φ(φd−)]."""
    d1, d2 = d_plus_minus(F, K, sigma, tau)
    return phi * (F * norm.cdf(phi * d1) - K * norm.cdf(phi * d2))


def premium(F, K, sigma, tau, df_quote, phi):
    """Discounted premium in quote currency per unit of base notional."""
    return df_quote * forward_premium(F, K, sigma, tau, phi)


def vega(F, K, sigma, tau, df_quote):
    """∂premium/∂σ = D_q F φ(d+) √τ (identical for calls and puts)."""
    d1, _ = d_plus_minus(F, K, sigma, tau)
    return df_quote * F * norm.pdf(d1) * np.sqrt(tau)


def delta(F, K, sigma, tau, phi, convention: DeltaConvention, df_base=1.0):
    """Delta in the stated convention.

    Pips spot delta: φ D_b Φ(φd+). Premium-adjusted spot delta:
    φ D_b (K/F) Φ(φd−). Forward versions drop D_b. ``df_base`` is ignored for
    forward deltas.
    """
    d1, d2 = d_plus_minus(F, K, sigma, tau)
    scale = df_base if convention.spot else 1.0
    if convention.premium_adjusted:
        return phi * scale * (K / F) * norm.cdf(phi * d2)
    return phi * scale * norm.cdf(phi * d1)


def atm_dns_strike(F, sigma, tau, convention: DeltaConvention):
    """Delta-neutral-straddle strike: F e^{+σ²τ/2} (pips) or F e^{−σ²τ/2} (premium-adjusted)."""
    sign = -1.0 if convention.premium_adjusted else 1.0
    return F * np.exp(sign * 0.5 * sigma**2 * tau)


def pa_call_delta_maximiser(F, sigma, tau):
    """Strike K* maximising the premium-adjusted call delta (result R5(b)).

    K* solves σ√τ Φ(d−) = φ(d−). The ratio φ/Φ is strictly decreasing, so the
    root d* is unique; then K* = F exp(−d* σ√τ − σ²τ/2).
    """
    _check_vol(sigma, tau)
    sd = sigma * np.sqrt(tau)
    # g(d) = sd Φ(d) − φ(d) has the sign of sd − φ(d)/Φ(d). The inverse Mills
    # ratio φ/Φ decreases strictly from +∞ to 0, so g has exactly one zero and
    # is negative to its left, positive to its right (theory note, R5(b)).
    g = lambda d: sd * norm.cdf(d) - norm.pdf(d)
    lo, hi = -1.0, 1.0
    for _ in range(MAX_EXPANSIONS):
        if g(lo) <= 0:
            break
        lo *= 2.0
    for _ in range(MAX_EXPANSIONS):
        if g(hi) >= 0:
            break
        hi *= 2.0
    d_star = brentq(g, lo, hi, xtol=1e-14, rtol=1e-14)
    return F * np.exp(-d_star * sd - 0.5 * sd**2)


def strike_from_delta(target, F, sigma, tau, phi, convention: DeltaConvention, df_base=1.0):
    """Strike with the given delta under a flat volatility (result R5).

    Pips deltas have the closed form of R5(a). Premium-adjusted put deltas are
    strictly monotone in K. Premium-adjusted call deltas are unimodal; the
    market convention takes the root on [K*, ∞), which is returned here.
    ``target`` is signed (negative for puts).
    """
    if np.sign(target) != phi or target == 0:
        raise ValueError("target delta must be non-zero with the sign of the option type")
    _check_vol(sigma, tau)
    scale = df_base if convention.spot else 1.0
    sd = sigma * np.sqrt(tau)
    if not convention.premium_adjusted:
        a = abs(target) / scale
        if not 0 < a < 1:
            raise ValueError("|delta| must lie in (0, D_b) for pips spot delta")
        return F * np.exp(0.5 * sd**2 - phi * sd * norm.ppf(a))

    f = lambda K: delta(F, K, sigma, tau, phi, convention, df_base) - target
    if phi == PUT:
        # |Δ| increases from 0 to ∞ as K increases.
        lo, hi = F * np.exp(-10 * sd - 1e-8), F
        for _ in range(MAX_EXPANSIONS):
            if f(hi) <= 0:
                break
            hi *= 2.0
        return brentq(f, lo, hi, xtol=1e-14 * F, rtol=1e-14)
    k_star = pa_call_delta_maximiser(F, sigma, tau)
    if f(k_star) < 0:
        raise ValueError("target exceeds the maximal premium-adjusted call delta")
    hi = k_star * 2.0
    for _ in range(MAX_EXPANSIONS):
        if f(hi) <= 0:
            break
        hi *= 2.0
    return brentq(f, k_star, hi, xtol=1e-14 * F, rtol=1e-14)


def strike_from_delta_smile(
    target,
    F,
    tau,
    phi,
    convention: DeltaConvention,
    vol_of_strike: Callable[[float], float],
    df_base=1.0,
    n_scan=801,
    span=12.0,
):
    """Strike solving Δ(K, σ(K)) = target for a smile σ(·) (fixed point in R8).

    With a smile the delta map need not be monotone (result R5(c)). The
    root nearest the flat-volatility strike is found by expanding a bracket
    around that strike. For premium-adjusted calls the root must lie where
    delta decreases in K (the conventional branch); otherwise, or if no local
    bracket is found, the function scans a log-strike grid of ±``span`` ATM
    standard deviations and returns the conventional root closest to the
    flat-volatility strike. Raises if no root exists.
    """
    sigma_ref = float(vol_of_strike(F))
    sd = sigma_ref * np.sqrt(tau)
    g = lambda K: float(delta(F, K, float(vol_of_strike(K)), tau, phi, convention, df_base)) - target
    k0 = strike_from_delta(target, F, sigma_ref, tau, phi, convention, df_base)
    g0 = g(k0)
    if g0 == 0.0:
        return k0
    pa_call = convention.premium_adjusted and phi == CALL
    step = 0.05 * sd
    for _ in range(30):
        lo, hi = k0 * np.exp(-step), k0 * np.exp(step)
        glo, ghi = g(lo), g(hi)
        root = None
        if np.sign(glo) != np.sign(g0):
            root = brentq(g, lo, k0, xtol=1e-13 * F, rtol=1e-14)
        elif np.sign(ghi) != np.sign(g0):
            root = brentq(g, k0, hi, xtol=1e-13 * F, rtol=1e-14)
        if root is not None:
            if not pa_call or g(root * (1 + 1e-6)) < g(root * (1 - 1e-6)):
                return root
            break  # root on the non-conventional branch: fall back to the scan
        step *= 1.8
        if step > span * sd:
            break
    ks = F * np.exp(np.linspace(-span * sd, span * sd, n_scan))
    try:  # vectorised evaluation when the smile accepts arrays (SABR does)
        vols = np.asarray(vol_of_strike(ks), dtype=float)
        if vols.shape != ks.shape:
            raise ValueError
    except (TypeError, ValueError):
        vols = np.array([vol_of_strike(k) for k in ks])
    vals = delta(F, ks, vols, tau, phi, convention, df_base) - target
    if convention.premium_adjusted and phi == CALL:
        deltas = vals + target
        i_max = int(np.argmax(deltas))
        ks, vals = ks[i_max:], vals[i_max:]
    idx = np.where(np.sign(vals[:-1]) * np.sign(vals[1:]) <= 0)[0]
    if idx.size == 0:
        raise ValueError("no strike attains the target delta on the scanned range")
    flat_guess = strike_from_delta(target, F, sigma_ref, tau, phi, convention, df_base)
    i = idx[np.argmin(np.abs(np.log(ks[idx] / flat_guess)))]
    g = lambda K: delta(F, K, vol_of_strike(K), tau, phi, convention, df_base) - target
    if g(ks[i]) == 0:
        return float(ks[i])
    return brentq(g, ks[i], ks[i + 1], xtol=1e-13 * F, rtol=1e-14)
