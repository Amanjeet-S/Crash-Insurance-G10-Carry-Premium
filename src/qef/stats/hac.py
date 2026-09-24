"""Newey–West HAC variance with the automatic bandwidth of Newey and West (1994).

Bartlett kernel without prewhitening (research design, section 7):

- pilot truncation n = ⌊4 (T/100)^{2/9}⌋;
- s0 = σ̂0 + 2 Σ_{j=1}^{n} σ̂j and s1 = 2 Σ_{j=1}^{n} j σ̂j, with σ̂j = T⁻¹ Σ u_t u_{t−j};
- bandwidth γ̂ T^{1/3} with γ̂ = 1.1447 (s1/s0)^{2/3};
- lag L = ⌊bandwidth⌋ and weights 1 − j/(L + 1) (Newey and West, 1987).

The constants follow Newey and West (1994, Tables I and II) as implemented in
the R package sandwich (functions bwNeweyWest and NeweyWest, version 3.1-3),
which served as the reference implementation. Two guards depart from it: if
s0 ≤ 0 the bandwidth is set to zero (the iid variance), and the lag is capped
at T − 1. Neither binds on the project's series. Full references are in
projects/01_fx_crash_insurance_carry/references.md.
"""

from __future__ import annotations

import numpy as np


def _autocov(u: np.ndarray, j: int) -> float:
    T = u.shape[0]
    return float(u[j:] @ u[: T - j]) / T


def nw_bandwidth(u) -> float:
    """Automatic Bartlett bandwidth of Newey and West (1994), no prewhitening."""
    u = np.asarray(u, dtype=float)
    T = u.shape[0]
    n = int(np.floor(4.0 * (T / 100.0) ** (2.0 / 9.0)))
    sig = np.array([_autocov(u, j) for j in range(n + 1)])
    s0 = sig[0] + 2.0 * sig[1:].sum()
    s1 = 2.0 * (np.arange(1, n + 1) * sig[1:]).sum()
    if s0 <= 0:
        return 0.0
    gamma = 1.1447 * ((s1 / s0) ** 2) ** (1.0 / 3.0)
    return gamma * T ** (1.0 / 3.0)


def long_run_variance(u, lag: int | None = None) -> tuple[float, int]:
    """Bartlett-weighted long-run variance of a demeaned series; returns (variance, lag)."""
    u = np.asarray(u, dtype=float)
    if lag is None:
        lag = int(np.floor(nw_bandwidth(u)))
    lag = max(0, min(lag, u.shape[0] - 1))
    v = _autocov(u, 0)
    for j in range(1, lag + 1):
        v += 2.0 * (1.0 - j / (lag + 1.0)) * _autocov(u, j)
    return v, lag


def mean_and_se(x) -> dict:
    """Sample mean with its Newey–West standard error."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    lrv, lag = long_run_variance(x - x.mean())
    return {"mean": float(x.mean()), "se": float(np.sqrt(max(lrv, 0.0) / x.shape[0])), "lag": lag, "n": int(x.shape[0])}


def ols_hac(y, X) -> dict:
    """OLS coefficients with a Newey–West covariance matrix.

    The automatic bandwidth is computed on the equally weighted sum of the
    score columns excluding the intercept, as in sandwich's default weights.
    """
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X = y[ok], X[ok]
    T, k = X.shape
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    e = y - X @ beta
    scores = X * e[:, None]
    w = np.ones(k)
    const = np.all(X == 1.0, axis=0)
    if const.any() and k > 1:
        w[const] = 0.0
    lag = int(np.floor(nw_bandwidth(scores @ w)))
    lag = max(0, min(lag, T - 1))
    S = scores.T @ scores / T
    for j in range(1, lag + 1):
        G = scores[j:].T @ scores[: T - j] / T
        S += (1.0 - j / (lag + 1.0)) * (G + G.T)
    V = T * XtX_inv @ S @ XtX_inv
    return {"beta": beta, "se": np.sqrt(np.diag(V)), "cov": V, "lag": lag, "n": T, "resid": e}
