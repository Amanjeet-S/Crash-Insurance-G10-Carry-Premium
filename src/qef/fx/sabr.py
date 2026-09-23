"""Hagan et al. (2002) lognormal implied volatility for the SABR model.

The expansion is an asymptotic approximation, not an exact model price, so
SABR is used as an interpolator of quoted smiles, never as a closed-form
benchmark (research design, section 12).
"""

from __future__ import annotations

import numpy as np


def _z_over_x(z, rho):
    """z / x(z) with x(z) = ln[(√(1 − 2ρz + z²) + z − ρ) / (1 − ρ)], continuous at z = 0."""
    z = np.asarray(z, dtype=float)
    out = np.empty_like(z)
    small = np.abs(z) < 1e-7
    zs = z[small]
    out[small] = 1.0 - 0.5 * rho * zs + (2.0 - 3.0 * rho**2) * zs**2 / 12.0
    zl = z[~small]
    x = np.log((np.sqrt(1.0 - 2.0 * rho * zl + zl**2) + zl - rho) / (1.0 - rho))
    out[~small] = zl / x
    return out


def sabr_vol(K, F, tau, alpha, rho, nu, beta=1.0):
    """Black (lognormal) implied volatility from the Hagan SABR expansion.

    Parameters satisfy alpha > 0, −1 < rho < 1, nu ≥ 0, 0 ≤ beta ≤ 1.
    """
    K = np.asarray(K, dtype=float)
    one_b = 1.0 - beta
    log_fk = np.log(F / K)
    fk_pow = (F * K) ** (0.5 * one_b)
    denom = fk_pow * (1.0 + one_b**2 / 24.0 * log_fk**2 + one_b**4 / 1920.0 * log_fk**4)
    z = (nu / alpha) * fk_pow * log_fk if nu > 0 else np.zeros_like(K)
    ratio = _z_over_x(z, rho) if nu > 0 else np.ones_like(K)
    correction = 1.0 + (
        one_b**2 / 24.0 * alpha**2 / fk_pow**2
        + 0.25 * rho * beta * nu * alpha / fk_pow
        + (2.0 - 3.0 * rho**2) / 24.0 * nu**2
    ) * tau
    return alpha / denom * ratio * correction
