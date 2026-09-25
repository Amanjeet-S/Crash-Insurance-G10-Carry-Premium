"""Option-implied variance and skewness of currency returns from truncated quotes (R1, R3, R4).

The object is y = ln(X_T/F_X), with X the USD price of one unit of the non-USD
currency (research design, section 3) and T the option expiry. Moments are
taken under the USD T-forward measure Q^USD. P(K) and C(K) are the undiscounted
out-of-the-money Garman–Kohlhagen premia on the quoted pair at the smile
volatility; they are expectations under the quote-currency T-forward measure
Q^q, under which E^q S_T = F.

Change of numeraire (R1). For a USD-base pair (USDJPY) take the quote currency
as domestic and USD as foreign, so S is the domestic price of one unit of
foreign currency. With the T-maturity zero-coupon bonds as numeraires (over one
month with deterministic rates they coincide with the money-market accounts,
B^d_T = 1/D_q and B^f_T = 1/D_b), the density of R1 is
Z_T = S_T B^f_T/(S_0 B^d_T) = S_T D_q/(S_0 D_b) = S_T/F, since F = S_0 D_b/D_q.
R1(b) at t = 0 gives

    E^USD h(S_T) = E^q[h(S_T) S_T/F],

and in particular E^USD(1/S_T) = 1/F, the USD forward price of the currency
(R1(d)). Under a flat smile ln(S_T/F) ~ N(−s²/2, s²) under Q^q, with s = σ√τ.
The factor S_T/F = exp(ln(S_T/F)) shifts the mean of this normal law by its
variance, so ln(S_T/F) ~ N(s²/2, s²) under Q^USD and y = −ln(S_T/F) ~ N(−s²/2, s²).

Contracts (R3), adapted from Bakshi, Kapadia and Madan (2003, Theorem 1): they
are centred at the forward F and priced under the forward measure rather than
expanded around the spot with a constant interest rate. The raw moment
E^USD y^k, k = 1, 2, 3, equals E^q h_k(S_T), with

- EURUSD-type pairs (USD is the quote currency, Q^q = Q^USD, y = ln(x/F)):
  h_k = f_k, f_k(x) = (ln(x/F))^k;
- USD-base pairs (y = −ln(x/F)): h_k = g_k, g_k(x) = (−ln(x/F))^k x/F.

With y(x) as above, y' = ±1/x, f_k' = k y^{k−1}/x and g_k' = (y^k − k y^{k−1})/F, so

    f_k''(K) = [k(k−1) y^{k−2} − k y^{k−1}] / K²,
    g_k''(K) = [k(k−1) y^{k−2} − k y^{k−1}] / (F K).

h_k(F) = 0, so R3 gives E^USD y^k = ∫_0^F h_k'' P dK + ∫_F^∞ h_k'' C dK. For
k = 1 the bracket is −1; for k ≥ 2 it equals k y^{k−2}(k − 1 − y), which
vanishes only at y = k − 1 and, for k = 3, at y = 0 (K = F). In a tail a weight
therefore changes sign at most once: at K = F e^{k−1} in the upper tail of an
EURUSD-type pair, at K = F e^{−(k−1)} in the lower tail of a USD-base pair.

Middle part. Between K_min and K_max, the 10Δ put and call strikes, the
integrals are computed by composite Simpson's rule with F a node and each side
integrated separately (R4 remark). The number of subintervals is doubled until
halving the step changes every contract by less than 1e-10 relative.

Tails (R4). A weight w is split into its positive and negative parts. Under the
power bounds of the lemma each part's tail integral lies in [0, B], with B from
the general inequalities of the R4 proposition,

    B = P(K_min) K_min^{−(γ+1)} ∫_0^{K_min} w_±(K) K^{γ+1} dK      (lower tail),
    B = C(K_max) K_max^{η−1} ∫_{K_max}^∞ w_±(K) K^{1−η} dK          (upper tail),

computed by adaptive quadrature. Each raw moment is then an interval, and the
variance v = m2 − m1² and skewness (m3 − 3 m1 m2 + 2 m1³)/v^{3/2} are reported as
their exact ranges over the box of contract intervals, never as point values.

Exponents. Setting (i), boundary elasticities: γ = K_min g(K_min)/G(K_min) and
η = K_max g(K_max)/Ḡ(K_max), with G = P', Ḡ = −C' and g = P'' = C'' from the
smile. G(u)/u^γ is nondecreasing on (0, K_min] exactly when the elasticity
u g(u)/G(u) is at least γ there, and u^η Ḡ(u) is nonincreasing on [K_max, ∞)
exactly when u g(u)/Ḡ(u) is at least η there. Setting (i) therefore assumes
that the tails beyond the last quotes are no heavier than the smile implies at
those quotes. Setting (ii), heavy tails: γ = η = 2. The upper bound needs η > 1
(for USD-base pairs the upper-tail weight decays only like (ln K)^{k−1}/K). Full
references are in projects/01_fx_crash_insurance_carry/references.md.
"""

from __future__ import annotations

import itertools
from typing import Callable

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.stats import norm

from .gk import CALL, PUT, DeltaConvention, d_plus_minus, forward_premium, strike_from_delta_smile

CONTRACTS = (1, 2, 3)
SIMPSON_RTOL = 1e-10  # largest relative change of a contract when the step is halved
QUAD_RTOL = 1e-12  # relative tolerance of the adaptive tail quadrature
QUAD_LIMIT = 500
HEAVY_TAIL_EXPONENT = 2.0  # setting (ii)


def _vols(vol_of_strike: Callable, K) -> np.ndarray:
    """Smile volatilities at an array of strikes, vectorised when the smile accepts arrays."""
    K = np.asarray(K, dtype=float)
    try:
        return np.broadcast_to(np.asarray(vol_of_strike(K), dtype=float), K.shape).copy()
    except (TypeError, ValueError):
        return np.array([float(vol_of_strike(k)) for k in K.ravel()]).reshape(K.shape)


def log_return(x, F, usd_base: bool):
    """y = ln(X_T/F_X) as a function of the terminal rate x on the quoted pair."""
    u = np.log(np.asarray(x, dtype=float) / F)
    return -u if usd_base else u


def contract_payoff(k: int, x, F, usd_base: bool):
    """h_k(x): f_k(x) = y^k for EURUSD-type pairs, g_k(x) = y^k x/F for USD-base pairs."""
    x = np.asarray(x, dtype=float)
    y = log_return(x, F, usd_base)
    return y**k * x / F if usd_base else y**k


def contract_weight(k: int, K, F, usd_base: bool):
    """h_k''(K) = [k(k−1) y^{k−2} − k y^{k−1}] / K² (EURUSD-type) or / (F K) (USD-base)."""
    K = np.asarray(K, dtype=float)
    y = log_return(K, F, usd_base)
    bracket = (k * (k - 1) * y ** (k - 2) if k >= 2 else 0.0) - k * y ** (k - 1)
    return bracket / (F * K) if usd_base else bracket / K**2


def weight_sign_changes(k: int, F, usd_base: bool) -> tuple:
    """Strikes other than F at which h_k'' changes sign (y = k − 1)."""
    if k < 2:
        return ()
    return (F * np.exp(-(k - 1.0)),) if usd_base else (F * np.exp(k - 1.0),)


def otm_price(K, F, tau, vol_of_strike):
    """Undiscounted out-of-the-money premium: the put below F, the call at and above F."""
    K = np.asarray(K, dtype=float)
    return forward_premium(F, K, _vols(vol_of_strike, K), tau, np.where(K < F, PUT, CALL))


def _simpson(values, h):
    return h / 3.0 * (values[0] + values[-1] + 4.0 * values[1:-1:2].sum() + 2.0 * values[2:-1:2].sum())


def middle_contracts(F, K_min, K_max, tau, vol_of_strike, usd_base: bool, n0=32, rtol=SIMPSON_RTOL, n_max=2**16):
    """Middle parts ∫_{K_min}^F h_k'' P dK + ∫_F^{K_max} h_k'' C dK for k = 1, 2, 3.

    Composite Simpson's rule with n subintervals on each side and F a node.
    Starting from n0, n is doubled until halving the step changes every
    contract by less than ``rtol`` relative. Returns the values at the final n,
    n, the largest relative change at the last halving and whether it is below
    ``rtol``.
    """
    if not K_min < F < K_max:
        raise ValueError("the middle part needs K_min < F < K_max")

    def integrate(n):
        total = np.zeros(len(CONTRACTS))
        for a, b in ((K_min, F), (F, K_max)):
            K = np.linspace(a, b, n + 1)
            p = otm_price(K, F, tau, vol_of_strike)
            total += [_simpson(contract_weight(k, K, F, usd_base) * p, (b - a) / n) for k in CONTRACTS]
        return total

    n, prev = n0, integrate(n0)
    while True:
        n *= 2
        cur = integrate(n)
        change = float(np.max(np.abs(cur - prev) / np.maximum(np.abs(cur), np.finfo(float).tiny)))
        if change < rtol or n >= n_max:
            return {"values": cur, "n": n, "rel_change": change, "converged": change < rtol}
        prev = cur


def tail_bounds(weight: Callable, side: str, K_edge, price_edge, exponent, breaks=()):
    """R4 bounds (B₊, B₋) on the tail integrals of the positive and negative parts of a weight.

    ``side="lower"``: the tail (0, K_min], with K_edge = K_min, price_edge =
    P(K_min) and exponent γ > 0. ``side="upper"``: the tail [K_max, ∞), with
    K_edge = K_max, price_edge = C(K_max) and exponent η > 1. With K = K_edge e^x
    the integrals of the R4 proposition become P(K_min) K_min ∫_{−∞}^0
    w±(K) e^{(γ+2)x} dx and C(K_max) K_max ∫_0^∞ w±(K) e^{(2−η)x} dx. ``breaks``
    are the strikes at which w changes sign; the tail is split there and each
    piece, of constant sign, is integrated by adaptive quadrature. Returns
    (B₊, B₋, largest quadrature error estimate relative to its integral).
    """
    if side == "lower":
        if not exponent > 0:
            raise ValueError("the lower-tail exponent must be positive")
        power = exponent + 2.0
        edges = [-np.inf, *sorted(np.log(b / K_edge) for b in breaks if 0 < b < K_edge), 0.0]
    elif side == "upper":
        if not exponent > 1:
            raise ValueError("the upper-tail exponent must exceed 1")
        power = 2.0 - exponent
        edges = [0.0, *sorted(np.log(b / K_edge) for b in breaks if b > K_edge), np.inf]
    else:
        raise ValueError(f"unknown side {side!r}")

    def integrand(x):
        with np.errstate(all="ignore"):
            v = float(weight(K_edge * np.exp(x))) * np.exp(power * x)
        # Non-finite values arise only where exp under- or overflows (|x| of several
        # hundred), where the integrand is below any representable contribution.
        return v if np.isfinite(v) else 0.0

    pos = neg = err = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        val, abserr = quad(integrand, lo, hi, epsabs=0.0, epsrel=QUAD_RTOL, limit=QUAD_LIMIT)
        if val >= 0:
            pos += val
        else:
            neg -= val
        if val != 0:
            err = max(err, abserr / abs(val))
    return price_edge * K_edge * pos, price_edge * K_edge * neg, err


def closed_form_bounds(F, K_min, P_min, gamma, K_max, C_max, eta) -> dict:
    """Closed forms of the R4 proposition for w = K^{−2} and the logarithmic weights.

    ``lower_log`` bounds ∫_0^{K_min} 2(1 + ln(F/K)) K^{−2} P dK; ``upper_log``
    bounds the integral of |2(1 − ln(K/F))| K^{−2} C over [K_max, ∞). The
    proposition states the latter for F < K_max ≤ eF only; outside that range
    it is NaN.
    """
    ell = np.log(K_max / F)
    upper_log = (2.0 * C_max / K_max * ((1.0 - ell) / eta - 1.0 / eta**2 + 2.0 * np.exp(-eta * (1.0 - ell)) / eta**2)
                 if 0.0 < ell <= 1.0 else np.nan)
    return {
        "lower_inverse_square": P_min / (gamma * K_min),
        "lower_log": 2.0 * P_min / K_min * ((1.0 + np.log(F / K_min)) / gamma + 1.0 / gamma**2),
        "upper_inverse_square": C_max / (eta * K_max),
        "upper_log": upper_log,
    }


def closed_form_check(F, K_min, P_min, gamma, K_max, C_max, eta) -> float:
    """Largest relative difference between ``tail_bounds`` and the applicable ``closed_form_bounds``."""
    cf = closed_form_bounds(F, K_min, P_min, gamma, K_max, C_max, eta)
    inv_sq = lambda K: K**-2.0
    num = {
        "lower_inverse_square": lambda: tail_bounds(inv_sq, "lower", K_min, P_min, gamma),
        "lower_log": lambda: tail_bounds(lambda K: 2.0 * (1.0 + np.log(F / K)) / K**2, "lower", K_min, P_min, gamma),
        "upper_inverse_square": lambda: tail_bounds(inv_sq, "upper", K_max, C_max, eta),
        "upper_log": lambda: tail_bounds(lambda K: 2.0 * (1.0 - np.log(K / F)) / K**2, "upper", K_max, C_max, eta, (np.e * F,)),
    }
    return max(abs(sum(num[k]()[:2]) / cf[k] - 1.0) for k in cf if np.isfinite(cf[k]))


def smile_distribution(K, F, tau, vol_of_strike, h=3e-4) -> dict:
    """P, C, G = P', Ḡ = −C' and g = P'' = C'' at strike K under the smile σ(K).

    For the undiscounted Black put p(K, σ): p_K = Φ(−d−), p_KK = φ(d−)/(Kσ√τ),
    p_σ = F φ(d+)√τ, p_Kσ = φ(d−) d+/σ and p_σσ = p_σ d+ d−/σ; the call has
    c_K = −Φ(d−) and the same second derivatives. Along the smile
    G = p_K + p_σ σ', Ḡ = Φ(d−) − p_σ σ' and
    g = p_KK + 2 p_Kσ σ' + p_σσ σ'² + p_σ σ''. σ' and σ'' are five-point
    differences in log strike with step h.
    """
    K = float(K)
    s = _vols(vol_of_strike, K * np.exp(h * np.arange(-2.0, 3.0)))
    sx = (s[0] - 8.0 * s[1] + 8.0 * s[3] - s[4]) / (12.0 * h)
    sxx = (-s[0] + 16.0 * s[1] - 30.0 * s[2] + 16.0 * s[3] - s[4]) / (12.0 * h**2)
    sig = float(s[2])
    ds, d2s = sx / K, (sxx - sx) / K**2
    d1, d2 = (float(d) for d in d_plus_minus(F, K, sig, tau))
    vega = F * norm.pdf(d1) * np.sqrt(tau)
    p_kk = norm.pdf(d2) / (K * sig * np.sqrt(tau))
    p_ks = norm.pdf(d2) * d1 / sig
    p_ss = vega * d1 * d2 / sig
    return {
        "P": float(forward_premium(F, K, sig, tau, PUT)),
        "C": float(forward_premium(F, K, sig, tau, CALL)),
        "G": float(norm.cdf(-d2) + vega * ds),
        "Gbar": float(norm.cdf(d2) - vega * ds),
        "g": float(p_kk + 2.0 * p_ks * ds + p_ss * ds**2 + vega * d2s),
    }


def boundary_elasticities(F, K_min, K_max, tau, vol_of_strike):
    """Setting (i): γ = K_min g(K_min)/G(K_min) and η = K_max g(K_max)/Ḡ(K_max)."""
    lo = smile_distribution(K_min, F, tau, vol_of_strike)
    hi = smile_distribution(K_max, F, tau, vol_of_strike)
    return K_min * lo["g"] / lo["G"], K_max * hi["g"] / hi["Gbar"]


def contract_intervals(middle, F, K_min, P_min, K_max, C_max, usd_base: bool, gamma, eta):
    """Raw-moment intervals: the middle part plus the R4 range [−B₋, B₊] of each tail.

    Returns (lo, hi, lower-tail bounds, upper-tail bounds, largest relative
    quadrature error estimate); the tail bounds are arrays of (B₊, B₋) per contract.
    """
    lo, hi = np.array(middle, dtype=float), np.array(middle, dtype=float)
    lower, upper = np.zeros((len(CONTRACTS), 2)), np.zeros((len(CONTRACTS), 2))
    err = 0.0
    for i, k in enumerate(CONTRACTS):
        w = lambda K, k=k: contract_weight(k, K, F, usd_base)
        breaks = weight_sign_changes(k, F, usd_base)
        for store, side, K_e, p_e, e in ((lower, "lower", K_min, P_min, gamma), (upper, "upper", K_max, C_max, eta)):
            b_pos, b_neg, qe = tail_bounds(w, side, K_e, p_e, e, breaks)
            store[i] = b_pos, b_neg
            lo[i] -= b_neg
            hi[i] += b_pos
            err = max(err, qe)
    return lo, hi, lower, upper, err


def variance_skewness(m1, m2, m3):
    """Variance m2 − m1² and skewness (m3 − 3 m1 m2 + 2 m1³)/v^{3/2} from raw moments."""
    v = m2 - m1**2
    return v, (m3 - 3.0 * m1 * m2 + 2.0 * m1**3) / v**1.5


def _variance(m):
    return m[1] - m[0] ** 2


def _variance_grad(m):
    return np.array([-2.0 * m[0], 1.0, 0.0])


def _skewness(m):
    return variance_skewness(*m)[1]


def _skewness_grad(m):
    # ∂s/∂m1 = 3(m1 m3 − m2²)/v^{5/2}, ∂s/∂m2 = 1.5(m1 m2 − m3)/v^{5/2}, ∂s/∂m3 = v^{−3/2}.
    m1, m2, m3 = m
    v = m2 - m1**2
    return np.array([3.0 * (m1 * m3 - m2**2), 1.5 * (m1 * m2 - m3), v]) / v**2.5


def _vertices(lo, hi):
    return [np.array(p, dtype=float) for p in itertools.product(*zip(lo, hi))]


def _variance_candidates(lo, hi):
    """v increases in m2 and, in m1, is largest at m1 = 0: vertices plus m1 = 0 if it lies inside."""
    pts = _vertices(lo, hi)
    if lo[0] < 0.0 < hi[0]:
        pts += [np.array([0.0, m2, lo[2]]) for m2 in (lo[1], hi[1])]
    return pts


def _skewness_candidates(lo, hi):
    """Vertices and edge critical points of the skewness.

    ∂s/∂m3 = v^{−3/2} > 0, so the extremes are at the endpoints of m3. The m1
    and m2 derivatives are proportional to m1 m3 − m2² and m1 m2 − m3; both
    vanish only if m2(m1² − m2) = 0, which is excluded when v > 0 (then
    m2 > m1² ≥ 0), so there is no critical point inside the (m1, m2)
    rectangle. On an edge m1 = a the only critical point is m2 = m3/a; on an
    edge m2 = b it is m1 = b²/m3.
    """
    pts = _vertices(lo, hi)
    for m3 in (lo[2], hi[2]):
        for m1 in (lo[0], hi[0]):
            if m1 != 0.0 and lo[1] < m3 / m1 < hi[1]:
                pts.append(np.array([m1, m3 / m1, m3]))
        for m2 in (lo[1], hi[1]):
            if m3 != 0.0 and lo[0] < m2**2 / m3 < hi[0]:
                pts.append(np.array([m2**2 / m3, m2, m3]))
    return pts


def box_range(fun, grad, lo, hi, candidates, n_random=4, seed=20260924):
    """Range of ``fun`` over the box [lo, hi]: candidate points plus a bounded multi-start search.

    The candidates (vertices and analytic critical points on the edges) give
    the exact range when the function has no other critical point in the box.
    Bounded L-BFGS-B, run in unit-cube coordinates from every vertex, the
    centre and ``n_random`` fixed-seed points, minimises and maximises
    independently; its optima are included. Returns (min, max, gap), where gap
    is how far the search went beyond the candidates (zero up to rounding).
    """
    lo, hi = np.asarray(lo, dtype=float), np.asarray(hi, dtype=float)
    width = hi - lo
    values = [fun(c) for c in candidates]
    c_min, c_max = min(values), max(values)
    rng = np.random.default_rng(seed)
    starts = [np.array(t, dtype=float) for t in itertools.product((0.0, 1.0), repeat=len(lo))]
    starts += [np.full(len(lo), 0.5), *rng.random((n_random, len(lo)))]
    o_min, o_max = c_min, c_max
    for sign in (1.0, -1.0):
        obj = lambda t: sign * fun(lo + t * width)
        jac = lambda t: sign * grad(lo + t * width) * width
        for t0 in starts:
            res = minimize(obj, t0, jac=jac, method="L-BFGS-B", bounds=[(0.0, 1.0)] * len(lo))
            val = fun(lo + np.clip(res.x, 0.0, 1.0) * width)
            o_min, o_max = min(o_min, val), max(o_max, val)
    return o_min, o_max, max(c_min - o_min, o_max - c_max)


def moment_ranges(lo, hi) -> dict:
    """Ranges of variance and skewness over the box of raw-moment intervals [lo, hi].

    If the box contains points with v ≤ 0, the skewness is unbounded and its
    range is reported as (−∞, ∞).
    """
    v_lo, v_hi, v_gap = box_range(_variance, _variance_grad, lo, hi, _variance_candidates(lo, hi))
    if not v_lo > 0:
        return {"var_lo": v_lo, "var_hi": v_hi, "skew_lo": -np.inf, "skew_hi": np.inf, "range_gap": v_gap}
    s_lo, s_hi, s_gap = box_range(_skewness, _skewness_grad, lo, hi, _skewness_candidates(lo, hi))
    return {"var_lo": v_lo, "var_hi": v_hi, "skew_lo": s_lo, "skew_hi": s_hi, "range_gap": max(v_gap, s_gap)}


def truncated_moment_intervals(F, tau, vol_of_strike, usd_base: bool, K_min, K_max, heavy=HEAVY_TAIL_EXPONENT) -> dict:
    """Contract, variance and skewness intervals from prices observed on [K_min, K_max].

    Settings: (i) boundary elasticities of the smile; (ii) γ = η = ``heavy``.
    Returns a flat dict with the quadrature check, the exponents and, per
    setting, the raw-moment intervals, the variance and skewness ranges, the
    tail-quadrature error, the closed-form check of the tail quadrature and a status.
    """
    mid = middle_contracts(F, K_min, K_max, tau, vol_of_strike, usd_base)
    P_min, C_max = float(otm_price(K_min, F, tau, vol_of_strike)), float(otm_price(K_max, F, tau, vol_of_strike))
    gamma, eta = boundary_elasticities(F, K_min, K_max, tau, vol_of_strike)
    out = {"K_min": K_min, "K_max": K_max, "P_min": P_min, "C_max": C_max,
           "n_simpson": mid["n"], "simpson_rel_change": mid["rel_change"], "simpson_converged": mid["converged"],
           **{f"mid_m{k}": v for k, v in zip(CONTRACTS, mid["values"])},
           "gamma_i": gamma, "eta_i": eta, "gamma_ii": heavy, "eta_ii": heavy}
    for name, (g, e) in (("i", (gamma, eta)), ("ii", (heavy, heavy))):
        if not (np.isfinite(g) and np.isfinite(e) and g > 0 and e > 1):
            out[f"status_{name}"] = "exponent_invalid"
            continue
        lo, hi, _, _, qerr = contract_intervals(mid["values"], F, K_min, P_min, K_max, C_max, usd_base, g, e)
        ranges = moment_ranges(lo, hi)
        for i, k in enumerate(CONTRACTS):
            out[f"m{k}_lo_{name}"], out[f"m{k}_hi_{name}"] = lo[i], hi[i]
        out.update({f"{key}_{name}": v for key, v in ranges.items()})
        out[f"tail_quad_err_{name}"] = qerr
        out[f"closed_form_rel_diff_{name}"] = closed_form_check(F, K_min, P_min, g, K_max, C_max, e)
        out[f"status_{name}"] = "ok" if ranges["var_lo"] > 0 else "variance_not_positive"
    return out


def implied_moment_intervals(F, tau, vol_of_strike, usd_base: bool, convention: DeltaConvention, df_base=1.0,
                             delta=0.10, heavy=HEAVY_TAIL_EXPONENT) -> dict:
    """``truncated_moment_intervals`` with K_min, K_max the smile ±``delta`` strikes in the pair's convention."""
    K_min = strike_from_delta_smile(-delta, F, tau, PUT, convention, vol_of_strike, df_base)
    K_max = strike_from_delta_smile(delta, F, tau, CALL, convention, vol_of_strike, df_base)
    return truncated_moment_intervals(F, tau, vol_of_strike, usd_base, K_min, K_max, heavy)
