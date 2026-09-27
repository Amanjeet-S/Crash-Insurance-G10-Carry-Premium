"""Stationary bootstrap (Politis and Romano, 1994) with automatic block length.

The expected block length follows Politis and White (2004) as corrected by
Patton, Politis and White (2009): b = (2 G² / D)^{1/3} T^{1/3} with
G = Σ_{|k|≤M} λ(k/M) |k| γ̂k, D = 2 ĝ0², ĝ0 = Σ_{|k|≤M} λ(k/M) γ̂k and the
flat-top kernel λ(x) = min(1, 2(1 − |x|)). M = 2 max(m̂, 1), where m̂ is the
smallest m such that the autocorrelations at lags m + 1, ..., m + K lie inside
±2 √(log10 T / T), with K = max(5, ⌊log10 T⌋) (Politis and White, 2004). The
caps M ≤ ⌈√T⌉ + K and b ≤ ⌈min(3√T, T/3)⌉ follow the arch package (Sheppard,
function optimal_block_length), which takes them from Andrew Patton's MATLAB
implementation; that package served as the reference implementation. Its
choice of m̂ sits one lag above the definition used here, and it normalises
autocorrelations by lag-specific sums of squares rather than by γ̂0, so block
lengths can differ from arch. A third difference concerns a degenerate case:
when the kernel estimate ĝ0 of the long-run variance is not positive, the
formula is not valid and this module returns b = 1 (the i.i.d. bootstrap),
whereas arch applies the formula and caps b. With the first two choices
switched to arch's (``_block_length(x, arch_conventions=True)``), the block
lengths equal those of arch 8.0.0 to rounding whenever ĝ0 > 0
(tests/test_bootstrap_arch.py).
"""

from __future__ import annotations

import numpy as np


def optimal_block_length(x) -> float:
    """Expected block length for the stationary bootstrap."""
    return _block_length(x, arch_conventions=False)


def _block_length(x, arch_conventions: bool) -> float:
    """Block length with this module's conventions, or with those of arch.

    ``arch_conventions=True`` switches only the two choices named in the module
    docstring: autocorrelations normalised by lag-specific sums of squares, as in
    arch's ``_single_optimal_block``, and m̂ taken as the first lag of the run of
    K insignificant autocorrelations rather than the lag before it. It exists so
    that a test can check every other step against arch exactly.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    T = x.shape[0]
    e = x - x.mean()
    b_max = np.ceil(min(3.0 * np.sqrt(T), T / 3.0))
    K = max(5, int(np.log10(T)))
    m_max = int(np.ceil(np.sqrt(T))) + K
    band = 2.0 * np.sqrt(np.log10(T) / T)
    acv = np.array([float(e[k:] @ e[: T - k]) / T for k in range(m_max + 1)])
    if arch_conventions:
        ss = lambda v: float(v @ v)
        rho = np.array([abs(T * acv[k]) / np.sqrt(ss(e[k + 1:]) * ss(e[: T - k - 1])) for k in range(m_max + 1)])
        window = lambda i: rho[i - K : i]
    else:
        rho = np.abs(acv / acv[0]) if acv[0] > 0 else np.zeros_like(acv)
        window = lambda i: rho[i - K + 1 : i + 1]
    m_hat = None
    for i in range(K, m_max + 1):
        if np.all(window(i) < band):
            m_hat = i - K
            break
    M = min(2 * max(m_hat, 1), m_max) if m_hat is not None else m_max
    k = np.arange(1, M + 1)
    lam = np.where(k / M <= 0.5, 1.0, 2.0 * (1.0 - k / M))
    G = 2.0 * np.sum(lam * k * acv[1 : M + 1])
    g0 = acv[0] + 2.0 * np.sum(lam * acv[1 : M + 1])
    if g0 <= 0 or G == 0:
        return 1.0
    b = (2.0 * G**2 / (2.0 * g0**2)) ** (1.0 / 3.0) * T ** (1.0 / 3.0)
    return float(min(max(b, 1.0), b_max))


def stationary_indices(T: int, b: float, rng: np.random.Generator) -> np.ndarray:
    """One stationary-bootstrap resample of 0..T-1: geometric blocks with mean b, wrapped circularly."""
    p = 1.0 / b
    idx = np.empty(T, dtype=int)
    idx[0] = rng.integers(T)
    new_block = rng.random(T) < p
    starts = rng.integers(T, size=T)
    for t in range(1, T):
        idx[t] = starts[t] if new_block[t] else (idx[t - 1] + 1) % T
    return idx


def bootstrap_distribution(stat, series: list, B: int = 9999, seed: int = 20260924, block: list | None = None) -> np.ndarray:
    """Bootstrap distribution of stat(*resampled series).

    Each series is resampled independently with its own expected block length
    (used for statistics computed on separate regime subsamples).
    """
    rng = np.random.default_rng(seed)
    series = [np.asarray(s, dtype=float) for s in series]
    blocks = block or [optimal_block_length(s) for s in series]
    out = np.empty(B)
    for r in range(B):
        draws = [s[stationary_indices(s.shape[0], b, rng)] for s, b in zip(series, blocks)]
        out[r] = stat(*draws)
    return out
