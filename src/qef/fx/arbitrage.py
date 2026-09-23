"""Static-arbitrage checks for a smile (result R2).

Two tests are applied: discrete convexity of undiscounted call prices on a
strike grid (slopes in [−1, 0] and nondecreasing), and the Gatheral–Jacquier
(2014) density condition g(k) ≥ 0 on implied total variance w(k) = σ²(k)τ,
with k = ln(K/F).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .gk import CALL, forward_premium


@dataclass(frozen=True)
class ArbitrageReport:
    convex_ok: bool
    slope_bounds_ok: bool
    min_g: float
    density_ok: bool
    worst_k: float


def g_function(k, w, w1, w2):
    """g(k) = (1 − k w′/(2w))² − (w′²/4)(1/w + 1/4) + w″/2."""
    return (1.0 - k * w1 / (2.0 * w)) ** 2 - (w1**2 / 4.0) * (1.0 / w + 0.25) + 0.5 * w2


def check_smile(vol_of_strike, F, tau, k_min=-0.5, k_max=0.5, n=2001, tol=1e-10) -> ArbitrageReport:
    # tol absorbs floating-point error in second differences of call prices.
    k = np.linspace(k_min, k_max, n)
    K = F * np.exp(k)
    sig = np.asarray(vol_of_strike(K), dtype=float)
    c = forward_premium(F, K, sig, tau, CALL) / F  # normalised undiscounted calls
    slopes = np.diff(c) / np.diff(K / F)
    slope_bounds_ok = bool(np.all(slopes >= -1 - tol) and np.all(slopes <= tol))
    convex_ok = bool(np.all(np.diff(slopes) >= -tol))
    w = sig**2 * tau
    h = k[1] - k[0]
    w1 = np.gradient(w, h)
    w2 = np.gradient(w1, h)
    g = g_function(k, w, w1, w2)
    inner = slice(2, -2)  # one-sided differences at the ends are less accurate
    i = int(np.argmin(g[inner])) + 2
    return ArbitrageReport(convex_ok, slope_bounds_ok, float(g[i]), bool(g[i] >= -1e-8), float(k[i]))
