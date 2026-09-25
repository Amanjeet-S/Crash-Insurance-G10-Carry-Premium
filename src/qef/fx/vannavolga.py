"""Vanna–volga smiles from ATM, risk-reversal and butterfly quotes (robustness variant 9).

Research design, section 8, variant 9, with the pillar and butterfly choices
recorded in the research log (24 September 2026). Quotes are in volatility
units and premiums are undiscounted, as in qef.fx.gk.

Method. Pillar strikes K1 < K2 < K3 carry volatilities σ1, σ2, σ3. Every
Greek is taken under the flat volatility σ = σ2 = σ_ATM (Castagna and
Mercurio, 2007, Section 3, p. 4; Castagna and Mercurio, 2006, fn 5). Three
smile functions are selected by ``method``:

``"price"``         the vanna–volga price, eq. (7), with the closed-form weights
                    of Prop. 3.1, eq. (6); its implied volatility is found by
                    inverting the Garman–Kohlhagen premium (p. 6).
``"first_order"``   the first-order approximation ς1, eq. (11), p. 10.
``"second_order"``  the second-order approximation ς2, eq. (12), p. 10.

Primary method: ``"price"``. The properties Castagna and Mercurio (2007)
establish belong to the price (7): pillar reproduction (Prop. 3.1, pp. 5–6),
strike-set invariance (Prop. 6.1), consistency with static replication
(Prop. 6.2) and the boundary limits of p. 6. Eqs. (11) and (12) are expansions
of the price in σ (Castagna and Mercurio, 2006, App. A, proofs of Props. 5.1
and 5.2), so they add approximation error exactly at the 10Δ hedge strikes,
which lie outside [K1, K3]. There ς1 overvalues the wings and violates the
bound of Lee (2004), as Castagna and Mercurio (2007, p. 10) state, and it is
the worst 10Δ extrapolator in Reiswich (2010, Section 4.2, Table 4.2). ς2 adds
a failure mode of its own, a negative radicand (Castagna and Mercurio, 2007,
pp. 10–11; Reiswich, 2010, Section 3.8.2, Tables 3.3–3.4). The price fails
only where it leaves the no-arbitrage bounds and no implied volatility exists
(Reiswich, 2010, Section 3.3.2), which is an arbitrage flag to be counted.
Reiswich (2010, Section 3.8) left the price out of his robustness study
because most of its calibrations missed his precision, so its calibration
robustness is not documented in the sources read; calibration status is
therefore recorded for every month.

Undefined values. ``vol`` returns NaN, and never raises or clips, where

- K is not positive and finite (all methods);
- ``"price"``: the premium of the out-of-the-money option is not strictly
  inside (0, F) for a call or (0, K) for a put, so no implied volatility
  exists (Reiswich, 2010, Section 3.3.2; convexity of (7) is not guaranteed,
  Castagna and Mercurio, 2007, fn 10), or the premium is not resolved in
  floating point (total volatility σ√τ outside TOTAL_VOL_BRACKET);
- ``"first_order"``: ς1 ≤ 0, which (11) permits (Reiswich, 2010, Section 3.3.2);
- ``"second_order"``: the radicand of (12) is negative (Castagna and Mercurio,
  2007, pp. 10–11), or ς2 ≤ 0.

Derived here, not stated in the sources:

- Forward form. Every term of (7) carries the quote-currency discount factor,
  and V(K)/V(K_i) = φ(d+(K))/φ(d+(K_i)) with d+ = [ln(F/K) + σ²τ/2]/(σ√τ),
  where F = S0 e^{(rd − rf)T}. Dividing (7) by the discount factor gives the
  same formula for undiscounted premiums.
- Puts. By put–call parity C^MKT(K_i) − C^BS(K_i) is the same for puts, so (7)
  applied to puts gives the same implied volatility. The out-of-the-money
  option is inverted, for accuracy.
- Eq. (12) is evaluated as ς2 = σ2 + b/(σ2 + √(σ2² + a b)), with a = d+ d− and
  b = 2σ2 D1 + D2. This equals (12) wherever a ≠ 0. At a = 0, where (12) is
  0/0, it equals the limit σ2 + D1 + D2/(2σ2). The ATM delta-neutral-straddle
  strike K2 is such a point in either delta convention (d+ = 0 for pips delta,
  d− = 0 for premium-adjusted delta), and there ς2 = σ2.
- At K1 and K3 the radicand of (12) is (σ2 + a(K_j)(σ_j − σ2))², so ς2
  reproduces σ_j only if σ2 + a(K_j)(σ_j − σ2) ≥ 0. Calibration checks the
  pillars, so a violation is reported, not assumed away.

Pillars from quotes. K2 is the delta-neutral-straddle strike at σ_ATM and K1,
K3 are the ∓Δ put and call strikes at the pillar volatilities
σ1,3 = σ_ATM + BF_ss ∓ RR/2, in the pair's delta convention (Castagna and
Mercurio, 2006, Section 2, eqs. (1)–(7), for pips spot delta; Reiswich, 2010,
Section 3.3.3, eq. (3.8), in the pair's convention; Bossens et al., 2010,
Section 3.1, eqs. (2)–(4), Section 3.2, eq. (5), and Section 3.3). BF_ss is
the smile strangle. The premium-adjusted call strike is the conventional-branch
root of ``strike_from_delta``. ATM and RR are matched by construction
(Reiswich, 2010, pp. 53–54).

Butterfly readings. ``"smile"``: BF_ss = BF (Castagna and Mercurio, 2006,
eqs. (4)–(5)). ``"market"``: BF_ss is the root of a one-dimensional search,
started at the quoted strangle, at which the smile prices the market strangle
at the market-strangle strikes (Reiswich, 2010, Section 3.3.3, eqs. (3.8)–(3.9);
Bossens et al., 2010, Section 3.3, pseudo-algorithm 2). The residual is in
volatility units, as in qef.fx.smile. A unique root is observed but not proved
(Reiswich, 2010, p. 54).

References (versions read):

- Castagna, A. and Mercurio, F. (2007). The vanna-volga method for implied
  volatilities. Risk, January 2007, pp. 106–111. Equation and page numbers
  refer to the authors' preprint of 30 October 2006.
- Castagna, A. and Mercurio, F. (2006). Consistent pricing of FX options.
  Working paper, Banca IMI, revised version of 1 September 2006.
- Bossens, F., Rayee, G., Skantzos, N. S. and Deelstra, G. (2010). Vanna-volga
  methods applied to FX derivatives: from theory to market practice.
  International Journal of Theoretical and Applied Finance 13(8), 1293–1324.
  arXiv:0904.1074v3 read.
- Reiswich, D. (2010). The Foreign Exchange Volatility Surface. Doctoral
  dissertation, Frankfurt School of Finance & Management.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq, elementwise
from scipy.stats import norm

from .gk import CALL, PUT, DeltaConvention, atm_dns_strike, d_plus_minus, delta, forward_premium, strike_from_delta, strike_from_delta_smile, vega
from .smile import SmileQuotes, market_strangle_strikes

PRICE, FIRST_ORDER, SECOND_ORDER = "price", "first_order", "second_order"
METHODS = (PRICE, FIRST_ORDER, SECOND_ORDER)
PRIMARY_METHOD = PRICE
TOTAL_VOL_BRACKET = (1e-9, 40.0)  # σ√τ bracket of the implied-volatility inversion
SCALAR_LIMIT = 4  # up to this many premiums are inverted one by one with brentq
TOL = 1e-8  # calibration success: max |residual| in volatility units, as in calibrate_sabr
_RTOL = 4 * np.finfo(float).eps
# Root searches on a smile with undefined regions fail with these; they are reported as NaN.
_SEARCH_ERRORS = (ValueError, RuntimeError, FloatingPointError, ZeroDivisionError)


def gk_implied_vol(p, F, K, tau, phi):
    """Volatility at which the undiscounted Garman–Kohlhagen premium equals ``p``.

    NaN where no volatility exists: p must lie strictly inside the no-arbitrage
    bounds max(φ(F − K), 0) < p < F (calls) or K (puts). The root is bracketed
    on σ√τ ∈ TOTAL_VOL_BRACKET (premiums outside that range are NaN) and found by
    Brent's method for a few premiums, or by scipy's elementwise bracketed solver.
    """
    p, K, phi = np.broadcast_arrays(np.asarray(p, float), np.asarray(K, float), np.asarray(phi))
    out = np.full(p.shape, np.nan)
    lo, hi = (s / np.sqrt(tau) for s in TOTAL_VOL_BRACKET)
    f = lambda sig, p_, K_, phi_: forward_premium(F, K_, sig, tau, phi_) - p_
    with np.errstate(all="ignore"):
        upper = np.where(phi == CALL, F, K)
        ok = np.isfinite(p) & np.isfinite(K) & (K > 0) & (p > np.maximum(phi * (F - K), 0.0)) & (p < upper)
        ok &= (f(lo, p, K, phi) < 0) & (f(hi, p, K, phi) > 0)
        idx = np.flatnonzero(ok)
        if idx.size == 0:
            return out
        pf, Kf, phf = p.ravel()[idx], K.ravel()[idx], phi.ravel()[idx]
        if idx.size <= SCALAR_LIMIT:
            vals = [brentq(lambda s: float(f(s, a, b, c)), lo, hi, xtol=1e-15, rtol=_RTOL) for a, b, c in zip(pf, Kf, phf)]
        else:
            res = elementwise.find_root(f, (lo, hi), args=(pf, Kf, phf))
            vals = np.where(res.success, res.x, np.nan)
    np.put(out, idx, vals)
    return out


@dataclass(frozen=True)
class VannaVolgaSmile:
    """Vanna–volga smile through three pillars; ``vol`` accepts scalars and arrays.

    ``strikes`` = (K1, K2, K3) with 0 < K1 < K2 < K3 and ``vols`` = (σ1, σ2, σ3);
    σ2 is the flat volatility of the Greeks. The constructor raises ValueError
    for pillars outside that domain, where the weights (6) are not defined
    (Castagna and Mercurio, 2006, App. A, eq. (23), needs K1 < K2 < K3).
    """

    F: float
    tau: float
    strikes: tuple
    vols: tuple
    method: str = PRIMARY_METHOD
    _cost: np.ndarray = field(init=False, repr=False, compare=False)  # C^MKT(K_i) − C^BS(K_i), undiscounted
    _dd: np.ndarray = field(init=False, repr=False, compare=False)  # d+ d− at the pillars, under σ2

    def __post_init__(self):
        if self.method not in METHODS:
            raise ValueError(f"method must be one of {METHODS}")
        K = np.asarray(self.strikes, dtype=float)
        s = np.asarray(self.vols, dtype=float)
        if K.shape != (3,) or s.shape != (3,):
            raise ValueError("three pillar strikes and volatilities are required")
        if not (np.isfinite(self.F) and self.F > 0 and np.isfinite(self.tau) and self.tau > 0):
            raise ValueError("forward and time to expiry must be positive and finite")
        if not (np.all(np.isfinite(K)) and 0 < K[0] < K[1] < K[2]):
            raise ValueError("pillar strikes must satisfy 0 < K1 < K2 < K3")
        if not np.all(np.isfinite(s) & (s > 0)):
            raise ValueError("pillar volatilities must be positive and finite")
        object.__setattr__(self, "strikes", tuple(float(k) for k in K))
        object.__setattr__(self, "vols", tuple(float(v) for v in s))
        phi = np.where(K >= self.F, CALL, PUT)  # out-of-the-money option at each pillar
        cost = forward_premium(self.F, K, s, self.tau, phi) - forward_premium(self.F, K, s[1], self.tau, phi)
        cost[1] = 0.0  # σ = σ2: the K2 term of (7) vanishes (Castagna and Mercurio, 2007, p. 6)
        d1, d2 = d_plus_minus(self.F, K, s[1], self.tau)
        object.__setattr__(self, "_cost", cost)
        object.__setattr__(self, "_dd", d1 * d2)

    def log_weights(self, K):
        """(y1, y2, y3)(K): the log-strike factors of eq. (6), and the coefficients of eq. (11)."""
        K1, K2, K3 = self.strikes
        lk = np.log(np.asarray(K, dtype=float))
        l1, l2, l3 = np.log(K1), np.log(K2), np.log(K3)
        return np.stack([(l2 - lk) * (l3 - lk) / ((l2 - l1) * (l3 - l1)),
                         (lk - l1) * (l3 - lk) / ((l2 - l1) * (l3 - l2)),
                         (lk - l1) * (lk - l2) / ((l3 - l1) * (l3 - l2))])

    def weights(self, K):
        """(x1, x2, x3)(K) of Prop. 3.1, eq. (6), with V(K)/V(K_i) = φ(d+(K))/φ(d+(K_i)) under σ2."""
        s2 = self.vols[1]
        d_k, _ = d_plus_minus(self.F, np.asarray(K, dtype=float), s2, self.tau)
        d_i, _ = d_plus_minus(self.F, np.asarray(self.strikes), s2, self.tau)
        ratio = norm.pdf(d_k)[None, ...] / norm.pdf(d_i).reshape((3,) + (1,) * d_k.ndim)
        return ratio * self.log_weights(K)

    def vv_premium(self, K, phi=CALL):
        """Undiscounted vanna–volga premium, eq. (7), for calls or puts; NaN for K not positive and finite.

        It does not depend on ``method`` and is returned even where it leaves
        the no-arbitrage bounds (there ``vol`` is NaN under ``"price"``).
        """
        K = np.asarray(K, dtype=float)
        with np.errstate(all="ignore"):
            valid = np.isfinite(K) & (K > 0)
            Ks = np.where(valid, K, self.F)
            out = forward_premium(self.F, Ks, self.vols[1], self.tau, phi) + np.tensordot(self._cost, self.weights(Ks), axes=1)
        out = np.where(valid, out, np.nan)
        return out[()] if out.ndim == 0 else out

    def vol(self, K, method: str | None = None):
        """Implied volatility at K under ``method`` (default: the smile's own); NaN where undefined."""
        m = self.method if method is None else method
        if m not in METHODS:
            raise ValueError(f"method must be one of {METHODS}")
        K = np.asarray(K, dtype=float)
        with np.errstate(all="ignore"):
            valid = np.isfinite(K) & (K > 0)
            Ks = np.where(valid, K, self.F)  # placeholder; invalid strikes are set to NaN below
            out = {PRICE: self._price_vol, FIRST_ORDER: self._first_order, SECOND_ORDER: self._second_order}[m](Ks)
        out = np.where(valid, out, np.nan)
        return out[()] if out.ndim == 0 else out

    def _price_vol(self, K):
        phi = np.where(K >= self.F, CALL, PUT)
        p = forward_premium(self.F, K, self.vols[1], self.tau, phi) + np.tensordot(self._cost, self.weights(K), axes=1)
        return gk_implied_vol(p, self.F, K, self.tau, phi)

    def _first_order(self, K):
        s = np.tensordot(np.asarray(self.vols), self.log_weights(K), axes=1)
        return np.where(s > 0, s, np.nan)

    def _second_order(self, K):
        s1, s2, s3 = self.vols
        y = self.log_weights(K)
        d1, d2 = d_plus_minus(self.F, K, s2, self.tau)
        a = d1 * d2
        D1 = np.tensordot(np.asarray(self.vols), y, axes=1) - s2
        D2 = y[0] * self._dd[0] * (s1 - s2) ** 2 + y[2] * self._dd[2] * (s3 - s2) ** 2
        b = 2.0 * s2 * D1 + D2
        rad = s2**2 + a * b
        s = s2 + b / (s2 + np.sqrt(np.where(rad >= 0, rad, np.nan)))
        return np.where(s > 0, s, np.nan)


def pillars(q: SmileQuotes, bf_ss: float):
    """Pillar strikes (K1, K2, K3) and volatilities (σ1, σ2, σ3) for the smile strangle ``bf_ss``.

    Raises ValueError where a pillar strike does not exist (a non-positive
    pillar volatility, or a premium-adjusted call target above the maximal delta).
    """
    s1 = q.sigma_atm + bf_ss - 0.5 * q.rr
    s3 = q.sigma_atm + bf_ss + 0.5 * q.rr
    k2 = float(atm_dns_strike(q.F, q.sigma_atm, q.tau, q.convention))
    k1 = float(strike_from_delta(-q.delta, q.F, s1, q.tau, PUT, q.convention, q.df_base))
    k3 = float(strike_from_delta(q.delta, q.F, s3, q.tau, CALL, q.convention, q.df_base))
    return (k1, k2, k3), (s1, q.sigma_atm, s3)


def vv_smile(q: SmileQuotes, bf_ss: float, method: str = PRIMARY_METHOD) -> VannaVolgaSmile:
    """Vanna–volga smile through the pillars of ``q`` at smile strangle ``bf_ss``; raises ValueError where undefined."""
    strikes, vols = pillars(q, bf_ss)
    return VannaVolgaSmile(q.F, q.tau, strikes, vols, method)


@dataclass(frozen=True)
class VVCalibration:
    smile: VannaVolgaSmile | None
    bf_ss: float
    residuals: np.ndarray  # (ATM, RR, strangle) in volatility units
    success: bool
    status: str  # ok, not_converged, quotes_not_reproduced, no_root, pillars_undefined, strangle_undefined
    strikes: dict  # K_call_ms, K_put_ms for the market strangle


def _strangle_residual(smile: VannaVolgaSmile, q: SmileQuotes, k_ms):
    """(smile strangle value − flat strangle value)/strangle vega at the market-strangle strikes, as in qef.fx.smile."""
    kc, kp = k_ms
    s = q.sigma_atm + q.bf
    target = forward_premium(q.F, kc, s, q.tau, CALL) + forward_premium(q.F, kp, s, q.tau, PUT)
    vc, vp = smile.vol(np.array([kc, kp]))
    value = forward_premium(q.F, kc, vc, q.tau, CALL) + forward_premium(q.F, kp, vp, q.tau, PUT)
    return float(value - target) / float(vega(q.F, kc, s, q.tau, 1.0) + vega(q.F, kp, s, q.tau, 1.0))


def quote_residuals(smile: VannaVolgaSmile, q: SmileQuotes, bf_reading: str, k_ms=None) -> np.ndarray:
    """(ATM, RR, BF) residuals in volatility units at the pillars; NaN where the smile is undefined.

    K2 is the smile's ATM strike and K1, K3 its ±Δ strikes: each pillar has the
    target delta at its own volatility, which the smile reproduces.
    """
    v1, v2, v3 = smile.vol(np.array(smile.strikes))
    r_atm, r_rr = v2 - q.sigma_atm, (v3 - v1) - q.rr
    if bf_reading == "smile":
        r_bf = 0.5 * (v1 + v3) - q.sigma_atm - q.bf
    else:
        r_bf = _strangle_residual(smile, q, k_ms)
    return np.array([r_atm, r_rr, r_bf], dtype=float)


def _bracket(f, x0, x_min, step=5e-4, n_expand=10, span=0.05, n_grid=101):
    """Bracket (a, b) of a sign change of f on (x_min, ∞) near x0, or None.

    The search first moves from x0 in the direction in which an increasing f
    crosses zero, doubling the step; if that meets an undefined value or no
    crossing, a grid of ±``span`` around x0 is scanned and the crossing nearest
    x0 is returned.
    """
    f0 = f(x0)
    if f0 == 0.0:
        return x0, x0
    if np.isfinite(f0):
        direction = -1.0 if f0 > 0 else 1.0
        x_prev, h = x0, step
        for _ in range(n_expand):
            x = x0 + direction * h
            if x <= x_min:
                x = x_min + 0.5 * (x_prev - x_min)
            fx = f(x)
            if not np.isfinite(fx):
                break
            if np.sign(fx) != np.sign(f0):
                return (x, x_prev) if x < x_prev else (x_prev, x)
            x_prev, h = x, 2.0 * h
    lo, hi = max(x0 - span, x_min + 1e-6), x0 + span
    if not lo < hi:
        return None
    grid = np.linspace(lo, hi, n_grid)
    vals = np.array([f(x) for x in grid])
    ok = np.isfinite(vals[:-1]) & np.isfinite(vals[1:]) & (np.sign(vals[:-1]) * np.sign(vals[1:]) <= 0)
    idx = np.flatnonzero(ok)
    if idx.size == 0:
        return None
    i = idx[np.argmin(np.abs(0.5 * (grid[idx] + grid[idx + 1]) - x0))]
    return grid[i], grid[i + 1]


def _solve(f, x0, x_min):
    """Root of f near x0 on (x_min, ∞) by bracketing and Brent's method; NaN if none is found."""
    br = _bracket(f, x0, x_min)
    if br is None:
        return np.nan
    if br[0] == br[1]:
        return float(br[0])
    try:
        return float(brentq(f, br[0], br[1], xtol=1e-15, rtol=_RTOL))
    except _SEARCH_ERRORS:  # an undefined value inside the bracket
        return np.nan


def calibrate_vv(q: SmileQuotes, bf_reading: str = "market", method: str = PRIMARY_METHOD) -> VVCalibration:
    """Vanna–volga smile for ATM, RR and BF quotes under a butterfly reading. Never raises on market data.

    ``"smile"``: the quoted butterfly is the smile strangle. ``"market"``: the
    smile strangle solves the market-strangle condition (Reiswich, 2010,
    eq. (3.9)); the root is bracketed from the quoted butterfly and refined by
    Brent's method. Each ``method`` is calibrated with its own smile, as in
    Reiswich (2010, Table 3.2).
    """
    if bf_reading not in ("market", "smile"):
        raise ValueError("bf_reading must be 'market' or 'smile'")
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}")
    nan3 = np.full(3, np.nan)

    def build(x):
        try:
            with np.errstate(all="ignore"):
                return vv_smile(q, x, method)
        except _SEARCH_ERRORS:
            return None

    if bf_reading == "smile":
        smile = build(q.bf)
        if smile is None:
            return VVCalibration(None, q.bf, nan3, False, "pillars_undefined", {})
        res = quote_residuals(smile, q, "smile")
        ok = bool(np.all(np.isfinite(res)) and np.max(np.abs(res)) < TOL)
        return VVCalibration(smile, q.bf, res, ok, "ok" if ok else "quotes_not_reproduced", {})

    try:
        k_ms = market_strangle_strikes(q)
    except _SEARCH_ERRORS:
        return VVCalibration(None, np.nan, nan3, False, "strangle_undefined", {})
    strikes = {"K_call_ms": float(k_ms[0]), "K_put_ms": float(k_ms[1])}

    def f(x):
        smile = build(x)
        if smile is None:
            return np.nan
        with np.errstate(all="ignore"):
            return _strangle_residual(smile, q, k_ms)

    x_min = -q.sigma_atm + 0.5 * abs(q.rr)  # both wing volatilities positive
    x = _solve(f, q.bf, x_min)
    if not np.isfinite(x):
        return VVCalibration(None, np.nan, nan3, False, "no_root", strikes)
    smile = build(x)
    if smile is None:
        return VVCalibration(None, x, nan3, False, "not_converged", strikes)
    with np.errstate(all="ignore"):
        res = quote_residuals(smile, q, "market", k_ms)
    ok = bool(np.all(np.isfinite(res)) and np.max(np.abs(res)) < TOL)
    return VVCalibration(smile, x, res, ok, "ok" if ok else "not_converged", strikes)


def _defined_between(smile: VannaVolgaSmile, a, b, n=129):
    """True if the smile is finite on a log-strike grid from a to b, both included."""
    ks = np.exp(np.linspace(np.log(a), np.log(b), n))
    return bool(np.all(np.isfinite(smile.vol(ks))))


def smile_delta_strike(smile: VannaVolgaSmile, target, phi, convention: DeltaConvention, df_base=1.0, tol=1e-10,
                       fallback=True, span=12.0, n_scan=801):
    """Strike with Δ(K, σ(K)) = target on a vanna–volga smile; NaN where none is found.

    Only the strikes connected to the ATM pillar K2 by a range on which the
    smile is defined are admitted: across a strike where the method is
    undefined the delta map is not continuous, so a root beyond it is not a
    strike of this smile. ``strike_from_delta_smile`` (qef.fx.gk) is tried
    first; it assumes a smile defined everywhere, so its result is kept only if
    σ(K) is finite, |Δ − target| ≤ ``tol``, K is connected to K2 and, for
    premium-adjusted calls, K is on the conventional (decreasing) branch. If
    that fails and ``fallback`` is set, the log-strike grid of
    ``strike_from_delta_smile`` (±``span`` ATM standard deviations,
    ``n_scan`` points), cut to the defined range around K2, is scanned with
    the same selection rule: for premium-adjusted calls only roots right of the
    largest delta, and the root nearest the flat-volatility strike at σ2.
    """
    F, tau, k2, s2 = smile.F, smile.tau, smile.strikes[1], smile.vols[1]
    pa_call = convention.premium_adjusted and phi == CALL
    dlt = lambda k: float(delta(F, k, float(smile.vol(k)), tau, phi, convention, df_base))

    def accept(k):
        if not (np.isfinite(k) and k > 0 and np.isfinite(smile.vol(k))):
            return False
        if not abs(dlt(k) - target) <= tol:
            return False
        if pa_call and not dlt(k * (1 + 1e-6)) < dlt(k * (1 - 1e-6)):
            return False
        return _defined_between(smile, k2, k)

    with np.errstate(all="ignore"):
        try:
            k = float(strike_from_delta_smile(target, F, tau, phi, convention, smile.vol, df_base))
        except _SEARCH_ERRORS:
            k = np.nan
        if accept(k):
            return k
        if not fallback:
            return np.nan
        sd = s2 * np.sqrt(tau)
        ks = np.unique(np.append(F * np.exp(np.linspace(-span * sd, span * sd, n_scan)), k2))
        vols = np.asarray(smile.vol(ks), dtype=float)
        finite = np.isfinite(vols)
        i2 = int(np.searchsorted(ks, k2))
        lo, hi = i2, i2  # the run of finite grid values around K2
        while lo > 0 and finite[lo - 1]:
            lo -= 1
        while hi < ks.size - 1 and finite[hi + 1]:
            hi += 1
        ks = ks[lo:hi + 1]
        d = delta(F, ks, vols[lo:hi + 1], tau, phi, convention, df_base)
        if pa_call:
            i_max = int(np.argmax(d))
            ks, d = ks[i_max:], d[i_max:]
        vals = d - target
        idx = np.flatnonzero(np.sign(vals[:-1]) * np.sign(vals[1:]) <= 0)
        if idx.size == 0:
            return np.nan
        try:
            guess = strike_from_delta(target, F, s2, tau, phi, convention, df_base)
        except ValueError:
            guess = k2
        i = idx[np.argmin(np.abs(np.log(ks[idx] / guess)))]
        try:
            k = float(brentq(lambda x: dlt(x) - target, ks[i], ks[i + 1], xtol=1e-13 * F, rtol=1e-14))
        except _SEARCH_ERRORS:
            return np.nan
        return k if accept(k) else np.nan


def implied_quotes(smile: VannaVolgaSmile, convention: DeltaConvention, delta=0.25, df_base=1.0, bf_reading="market") -> dict:
    """ATM, risk-reversal and butterfly quotes a vanna–volga smile implies at ``delta``; NaN where undefined.

    The definitions are those of qef.fx.smile.quotes_from_smile, so vanna–volga
    and SABR predictions of held-out quotes are comparable: the risk reversal
    and smile strangle use the smile's own ±Δ strikes; the market strangle is
    the BF at which the flat volatility σ_ATM + BF, at its own ±Δ strikes,
    prices the strangle as the smile does. The ATM strike is K2, which is the
    delta-neutral-straddle strike at σ(K2) = σ2. The market-strangle root is
    bracketed from the smile strangle outwards, so an undefined region far in
    the wings does not prevent it.
    """
    F, tau = smile.F, smile.tau
    s_atm = float(smile.vol(smile.strikes[1]))
    kc = smile_delta_strike(smile, delta, CALL, convention, df_base)
    kp = smile_delta_strike(smile, -delta, PUT, convention, df_base)
    vc, vp = smile.vol(np.array([kc, kp]))
    out = {"atm": s_atm, "rr": float(vc - vp), "bf": np.nan, "K_call": kc, "K_put": kp}
    bf_ss = 0.5 * (vc + vp) - s_atm
    if bf_reading == "smile":
        out["bf"] = float(bf_ss)
        return out
    if not np.isfinite(bf_ss):
        return out

    def gap(bf):
        try:
            kcm, kpm = market_strangle_strikes(SmileQuotes(F, tau, s_atm, out["rr"], bf, delta, df_base, convention))
        except _SEARCH_ERRORS:
            return np.nan
        s = s_atm + bf
        flat = forward_premium(F, kcm, s, tau, CALL) + forward_premium(F, kpm, s, tau, PUT)
        vcm, vpm = smile.vol(np.array([kcm, kpm]))
        value = forward_premium(F, kcm, vcm, tau, CALL) + forward_premium(F, kpm, vpm, tau, PUT)
        return float(flat - value) / float(vega(F, kcm, s, tau, 1.0) + vega(F, kpm, s, tau, 1.0))

    with np.errstate(all="ignore"):
        root = _solve(gap, float(bf_ss), -s_atm)  # the strangle volatility s_atm + bf stays positive
        if np.isfinite(root) and abs(gap(root)) <= TOL:
            out["bf"] = root
    return out
