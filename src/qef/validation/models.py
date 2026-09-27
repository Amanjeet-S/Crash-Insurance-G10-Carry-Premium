"""Independent option-pricing models for validating the implied-moment code (R1, R3, R4).

Two models of the quoted pair supply exact moments, QuantLib prices, a smile and
the tail behaviour against which qef.fx.moments is checked. The object is
x = ln(S_T/F) under the quote-currency T-forward measure Q^q, under which
E^q e^x = 1. Its cumulant generating function K_q(z) = ln M(z), M(z) = E^q e^{zx},
is finite for real z in an interval (w₋, w₊) containing [0, 1] and analytic in
the strip w₋ < Re z < w₊. The market is a spot S_0 and flat continuously
compounded rates of the quote and base currencies, D_q = e^{−r_q τ},
D_b = e^{−r_b τ}, F = S_0 D_b/D_q, with an Actual/360 day count and a 30-day
expiry, so that τ = 30/360 = 1/12 in QuantLib as well.

Reference moments (R1). y = ln(X_T/F_X), with X the USD price of the non-USD
currency. For EURUSD-type pairs y = x and Q^USD = Q^q, so the cumulants of y are
K_q^(n)(0). For USD-base pairs y = −x, and R1(b) with the density S_T/F gives
E^USD e^{vy} = E^q e^{(1−v)x}, so the cumulants of y are (−1)^n K_q^(n)(1). The raw
moments are m1 = κ1, m2 = κ2 + κ1² and m3 = κ3 + 3κ1κ2 + κ1³.

Merton (1976). Eqs. (2)–(3): S(t)/S = exp[(α − σ²/2 − λk)t + σZ(t)] Y(n), with
Poisson arrivals of intensity λ and k = E(Y − 1). In the lognormal case (p. 135,
before eq. (18)) ln Y ~ N(γ − δ²/2, δ²) with γ = ln(1 + k); μ = γ − δ²/2 is the
mean log jump, so k = e^{μ+δ²/2} − 1. The forward has no drift under Q^q, so

    x = −σ²τ/2 − λkτ + σW_τ + Σ_{i≤N} ln Y_i,
    K_q(z) = z(−σ²τ/2 − λkτ) + σ²τz²/2 + λτ(e^{zμ+z²δ²/2} − 1),

whose derivatives are in closed form. Given N = n, x ~ N(m_n, s_n²) with
m_n = −σ²τ/2 − λkτ + nμ and s_n² = σ²τ + nδ²; the law of x is the Poisson(λτ)
mixture of these normals, and log f, log F and log(1 − F) are log-sum-exps over n.
Eq. (19) prices a call as Σ_n e^{−λ'τ}(λ'τ)^n/n! f_n with λ' = λ(1 + k) and f_n the
Black–Scholes price at variance rate σ² + nδ²/τ and interest rate r − λk + nγ/τ.
Since e^{−λ'τ}(λ'τ)^n/n! e^{λkτ − nγ} = e^{−λτ}(λτ)^n/n!, the undiscounted form is

    C(K) = Σ_n e^{−λτ}(λτ)^n/n! B(F_n, K, s_n²),  F_n = F e^{−λkτ + nγ},

with B the undiscounted Black call; the put follows term by term, since
Σ_n e^{−λτ}(λτ)^n/n! F_n = F. QuantLib 1.43 does not expose its jump-diffusion
engine in Python, so QuantLib prices come from BatesModel and BatesEngine with
constant variance: v0 = θ = σ², vol of variance 1e−8 and ρ = 0. QuantLib's Bates
jump has ln(1 + J) ~ N(ν, δ²), so ν = μ; the tests confirm this against eq. (19).

Heston (1993). Eqs. (1) and (4) give dS = μS dt + √v S dz1 and
dv = κ(θ − v)dt + σ√v dz2 with correlation ρ; for currency options, Section 2
(eq. (25)) works with x = ln[S F(t;T)/P(t;T)], the log forward. With risk-neutral
parameters (no price of volatility risk) and deterministic rates, under Q^q

    dx = −v/2 dt + √v dW1,  dv = κ(θ − v)dt + ξ√v dW2,  d⟨W1, W2⟩ = ρ dt.

The characteristic function is eq. (17). I use the second representation φ2 of
Albrecher, Mayer, Schoutens and Tistaert (2006, Section 3, with g2 as in eq. (2)),
which they prove does not cross the branch cut of the logarithm on the lines
used for Fourier pricing (Theorem 3). With u = −iz,

    b = κ − ρξz,  d = √(b² − ξ²(z² − z)),  g = (b − d)/(b + d),
    K_q(z) = κθ/ξ² [(b − d)τ − 2 ln((1 − g e^{−dτ})/(1 − g))]
             + v0 (b − d)(1 − e^{−dτ}) / (ξ²(1 − g e^{−dτ})).

M(w) is finite for real w exactly when τ < T*(w), where T*(w) = ∞ if
Δ(w) = (κ − ρξw)² − ξ²(w² − w) ≥ 0 and otherwise, with ω = √(−Δ(w)) and
χ(w) = ρξw − κ, T*(w) = 2 arctan2(ω, χ)/ω (Keller-Ressel 2011, Section 6.1,
eq. (6.2), after Andersen and Piterbarg 2007; arctan2 combines the two branches
of eq. (6.2), which assumes κ > ρξ; Lord and Kahl 2006, Section 3.1.1, note that
the formula also holds for w < 0). The critical moments w₋ < 0 and w₊ > 1 solve
T*(w) = τ. The derivatives of K_q come from Cauchy's integral formula on a circle
of radius r about c, discretised by the trapezoidal rule with m nodes,
K^(n)(c) ≈ n!/(m rⁿ) Σ_j K_q(c + r e^{iθ_j}) e^{−inθ_j} (Bornemann 2011, eqs.
(1.2)–(1.3), after Lyness 1967), with error at most 2E_{m−1}/rⁿ, where E_{m−1} is
the error of the best polynomial approximation of degree m − 1 on the disc
(Theorem 2.1, eq. (2.2)). Here r = 1/2 and m = 64.

Fourier inversion. On a line w = c + iu inside the strip, with k = ln(K/F),

    f(k)          =  (1/π) ∫_0^∞ Re[e^{−wk} M(w)] du                      (any c),
    Q^q(x > k)    =  (1/π) ∫_0^∞ Re[e^{−wk} M(w)/w] du                    (c > 0),
    Q^q(x ≤ k)    = −(1/π) ∫_0^∞ Re[e^{−wk} M(w)/w] du                    (c < 0),
    price(K)/F    =  (1/π) ∫_0^∞ Re[e^{−(w−1)k} M(w)/(w(w − 1))] du        (call if c > 1, put if c < 0).

The first is Fourier inversion of e^{ck} f(k). The tail probabilities are
Lee (2004b, Theorem 5.1, eq. (5.4)) for the payoff G3 with b0 = 1 and b1 = 0,
whose residue term is 1 when the line passes to the left of the pole at w = 0.
The call is the Carr–Madan formula as written by Albrecher et al. (2006,
eqs. (3)–(5)) with damping α = c − 1; the put is Lee's eq. (5.4) for G1 with
α < −1, where the residue term F − K turns the call into the put. c minimises
the modulus of the integrand at u = 0 over the admissible side (Lord and Kahl
2006, Section 3.2, eq. (61)); for the density this is the saddle point
K_q'(c) = k. On that line the integrand has no large oscillating part, so prices
and probabilities far below 1e−16 keep their relative accuracy. The integral is
taken by the trapezoidal rule in u, truncated where the modulus of the integrand
falls below 1e−18 of its value at u = 0, and the step is halved until the value
changes by less than 1e−12 relative.

QuantLib prices are NPV/D_q, undiscounted. Heston prices use AnalyticHestonEngine
with the angled-contour complex logarithm and exp-sinh quadrature at relative
tolerance 1e−12; Merton prices use BatesEngine with its default 144-point
Gauss–Laguerre quadrature. Measured against the independent prices above for the
parameters of the tests, the absolute error is below 1e−14 F for |ln(K/F)| ≤ 0.8
(below 3e−16 F for Heston); BatesEngine degrades beyond, to about 1e−11 F at
|ln(K/F)| = 1 and 3e−7 F at 1.3, so the smile does not use nodes there.

Smile. Nodes x_j = j h, h = 1e−3, |x_j| ≤ 0.8, kept outwards from the money while
the QuantLib out-of-the-money price exceeds 1e−13 F. At each node the Black
volatility solves forward_premium(F, K, σ, τ, φ) = price (the undiscounted
Garman–Kohlhagen premium of qef.fx.gk, which the moment code uses to turn
volatilities back into prices) by bisection in ln σ; vol_of_strike is the
not-a-knot cubic spline of these volatilities in x = ln(K/F), held flat beyond the
outermost nodes. It reproduces the model prices exactly at the nodes; the tests
measure its error at the midpoints between nodes and propagate it to the moments.

Hazard rates (R4). With G(u) = Q^q(S_T ≤ u), Ḡ = 1 − G and g = G', u g(u)/G(u) = r(x)
and u g(u)/Ḡ(u) = h(x) at x = ln(u/F), where r = f/F is the reversed hazard rate
and h = f/(1 − F) the hazard rate of x. As noted in qef.fx.moments, the R4 lemma
holds with exponent γ on (0, K_min] exactly when γ ≤ inf{r(x): x ≤ ln(K_min/F)},
and with η on [K_max, ∞) exactly when η ≤ inf{h(x): x ≥ ln(K_max/F)}.
tail_exponent evaluates these infima on a grid that is fine (step 0.002) within
0.5 of the boundary and geometric out to 40 beyond it, and refines the smallest
grid value by bounded minimisation. Merton's r and h tend to ∞ in the tails; for
Heston, r → −w₋ and h → w₊ (the tails of x are exponential with the critical
moments as rates), so the grid reaches far enough to show the approach.

Price bounds. Lee (2004a, Theorem 2.1, eqs. (2.1)–(2.2)), undiscounted: for all
K > 0, q > 0 and p > 0,

    P(K) ≤ E S_T^{−q} q^q (1 + q)^{−(1+q)} K^{1+q},  C(K) ≤ E S_T^{p+1} p^p (p + 1)^{−(p+1)} K^{−p},

with E S_T^w = F^w M(w). price_bound minimises the bound at the edge strike over
the order and returns it with the exponent in the form of the R4 lemma, so that
qef.fx.moments.tail_bounds bounds the tails of a spanning integral beyond a
strike range.

References (version consulted):
- Albrecher, H., Mayer, P., Schoutens, W. and Tistaert, J. (2006). The little
  Heston trap. Working paper, version of 11 September 2006 (first version
  6 December 2005); published in Wilmott Magazine (2007). Full text: Sections
  2–5, eqs. (1)–(8), Theorem 3.
- Andersen, L. B. G. and Piterbarg, V. V. (2007). Moment explosions in stochastic
  volatility models. Finance and Stochastics 11(1), 29–50. Not consulted; cited as
  the source of the explosion time through Keller-Ressel (2011) and Lord and Kahl
  (2006).
- Bornemann, F. (2011). Accuracy and stability of computing high-order
  derivatives of analytic functions by Cauchy integrals. Foundations of
  Computational Mathematics 11(1), 1–63. arXiv:0910.1841v2 (20 November 2009),
  Sections 1–2.
- Heston, S. L. (1993). A closed-form solution for options with stochastic
  volatility with applications to bond and currency options. Review of Financial
  Studies 6(2), 327–343. Published article: eqs. (1), (4), (10), (17), (18), (25).
- Keller-Ressel, M. (2011). Moment explosions and long-term behavior of affine
  stochastic volatility models. Mathematical Finance 21(1), 73–98.
  arXiv:0802.1823v2 (13 October 2008), Section 6.1, eqs. (6.1)–(6.2).
- Lee, R. W. (2004a). The moment formula for implied volatility at extreme
  strikes. Mathematical Finance 14(3), 469–480. Preprint of 11 August 2003,
  Section 2.1, Theorem 2.1.
- Lee, R. W. (2004b). Option pricing by transform methods: extensions,
  unification, and error control. Journal of Computational Finance 7(3), 51–86.
  Author's version of 2 February 2005, Sections 2, 4 and 5, Theorems 4.2, 4.3
  and 5.1.
- Lord, R. and Kahl, C. (2007). Optimal Fourier inversion in semi-analytical
  option pricing. Journal of Computational Finance 10(4), 1–30. Tinbergen
  Institute Discussion Paper TI 2006-066/2 (2006), Sections 3.1–3.2,
  eqs. (43)–(45) and (60)–(62).
- Lyness, J. N. (1967). Numerical algorithms based on the theory of complex
  variable. Not consulted; cited through Bornemann (2011).
- Merton, R. C. (1976). Option pricing when underlying stock returns are
  discontinuous. Journal of Financial Economics 3(1–2), 125–144. Published
  article: eqs. (2), (3), (16)–(19).
"""

from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from typing import Callable

import numpy as np
import QuantLib as ql
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar
from scipy.special import log_ndtr, logsumexp, ndtr
from scipy.stats import poisson

from qef.fx.gk import CALL, PUT, forward_premium

EVALUATION_DATE = ql.Date(30, ql.September, 2026)
DAY_COUNT = ql.Actual360()
BATES_VOL_OF_VARIANCE = 1e-8  # constant variance: Merton through the Bates engine
CAUCHY_RADIUS = 0.5
CAUCHY_NODES = 64
LINE_RTOL = 1e-12  # relative change of a line integral at the last halving of the step
LINE_CUTOFF = 1e-18  # truncate a line integral where the integrand falls below this fraction of its value at u = 0
MAX_HALVINGS = 14
SMILE_SPACING = 1e-3  # node spacing in ln(K/F)
SMILE_REACH = 0.8  # largest |ln(K/F)| of a node
SMILE_FLOOR = 1e-13  # a node needs an out-of-the-money price above SMILE_FLOOR·F
TAIL_NEAR, TAIL_NEAR_STEP, TAIL_FAR, TAIL_FAR_POINTS = 0.5, 0.002, 40.0, 60


@dataclass(frozen=True)
class Market:
    """A currency pair: spot and flat continuously compounded quote- and base-currency rates."""

    spot: float
    rate_quote: float
    rate_base: float
    days: int = 30

    @property
    def tau(self) -> float:
        return self.days / 360.0  # Actual/360, as in the QuantLib curves

    @property
    def df_quote(self) -> float:
        return math.exp(-self.rate_quote * self.tau)

    @property
    def df_base(self) -> float:
        return math.exp(-self.rate_base * self.tau)

    @property
    def forward(self) -> float:
        return self.spot * self.df_base / self.df_quote

    def quantlib(self):
        """(expiry, quote-currency curve, base-currency curve, spot) for QuantLib; sets the evaluation date."""
        ql.Settings.instance().evaluationDate = EVALUATION_DATE
        curve = lambda r: ql.YieldTermStructureHandle(ql.FlatForward(EVALUATION_DATE, r, DAY_COUNT))
        return (EVALUATION_DATE + self.days, curve(self.rate_quote), curve(self.rate_base),
                ql.QuoteHandle(ql.SimpleQuote(self.spot)))


@dataclass(frozen=True)
class Merton:
    """Lognormal jump diffusion: diffusive volatility, jump intensity and mean and s.d. of the log jump."""

    sigma: float
    lam: float
    mu: float
    delta: float

    @property
    def k(self) -> float:
        return math.expm1(self.mu + 0.5 * self.delta**2)

    def strip(self, tau):
        return -math.inf, math.inf

    def cgf(self, z, tau):
        z = np.asarray(z)
        s2 = self.sigma**2 * tau
        jump = np.expm1(z * self.mu + 0.5 * z * z * self.delta**2)
        return z * (-0.5 * s2 - self.lam * tau * self.k) + 0.5 * s2 * z * z + self.lam * tau * jump

    def cgf_derivatives(self, c, tau):
        """K_q', K_q'' and K_q''' at real c in closed form."""
        s2, lt, b = self.sigma**2 * tau, self.lam * tau, self.delta**2
        p, e = self.mu + c * b, math.exp(c * self.mu + 0.5 * c * c * b)
        return np.array([-0.5 * s2 - lt * self.k + s2 * c + lt * p * e,
                         s2 + lt * (p * p + b) * e,
                         lt * (p**3 + 3.0 * p * b) * e])

    def _mixture(self, tau, n):
        m = -0.5 * self.sigma**2 * tau - self.lam * tau * self.k + n * self.mu
        return poisson.logpmf(n, self.lam * tau), m, np.sqrt(self.sigma**2 * tau + n * self.delta**2)

    def log_distribution(self, x, tau, n0=64):
        """log f, log F and log(1 − F) of x under Q^q from the Poisson mixture of normals.

        The number of terms is doubled until the last one is below e^{−60} of each total.
        """
        x = np.atleast_1d(np.asarray(x, dtype=float))
        n_terms = n0
        while True:
            lw, m, s = self._mixture(tau, np.arange(n_terms)[:, None])
            z = (x - m) / s
            parts = (lw - 0.5 * z * z - np.log(s) - 0.5 * math.log(2.0 * math.pi), lw + log_ndtr(z), lw + log_ndtr(-z))
            totals = [logsumexp(p, axis=0) for p in parts]
            if max(float(np.max(p[-1] - t)) for p, t in zip(parts, totals)) < -60.0:
                return tuple(totals)
            n_terms *= 2

    def log_density_and_tail(self, x, tau, side):
        lf, lF, lS = self.log_distribution(x, tau)
        return lf, (lF if side == "lower" else lS)

    def series_price(self, K, F, tau, n0=64):
        """Undiscounted out-of-the-money prices from Merton's eq. (19) in forward form.

        Terms are added in blocks of n0 until a block changes no price by more than 1e−17 relative.
        """
        K = np.atleast_1d(np.asarray(K, dtype=float))
        phi = np.where(K < F, PUT, CALL)
        gam = math.log1p(self.k)
        total, start = np.zeros_like(K), 0
        while True:
            n = np.arange(start, start + n0)[:, None]
            w = poisson.pmf(n, self.lam * tau)
            s = np.sqrt(self.sigma**2 * tau + n * self.delta**2)
            F_n = F * np.exp(-self.lam * self.k * tau + n * gam)
            d1 = (np.log(F_n / K) + 0.5 * s * s) / s
            block = (w * phi * (F_n * ndtr(phi * d1) - K * ndtr(phi * (d1 - s)))).sum(axis=0)
            total += block
            start += n0
            if np.all(np.abs(block) <= 1e-17 * np.abs(total)):
                return total

    def engine(self, market: Market):
        _, quote, base, spot = market.quantlib()
        v = self.sigma**2
        process = ql.BatesProcess(quote, base, spot, v, 1.0, v, BATES_VOL_OF_VARIANCE, 0.0, self.lam, self.mu, self.delta)
        return ql.BatesEngine(ql.BatesModel(process))


@dataclass(frozen=True)
class Heston:
    """Heston model of x: initial variance, mean reversion, long-run variance, vol of variance, correlation."""

    v0: float
    kappa: float
    theta: float
    xi: float
    rho: float

    def cgf(self, z, tau):
        z = np.asarray(z, dtype=complex)
        b = self.kappa - self.rho * self.xi * z
        d = np.sqrt(b * b - self.xi**2 * (z * z - z))
        g = (b - d) / (b + d)
        e = np.exp(-d * tau)
        A = self.kappa * self.theta / self.xi**2 * ((b - d) * tau - 2.0 * np.log((1.0 - g * e) / (1.0 - g)))
        B = (b - d) / self.xi**2 * (1.0 - e) / (1.0 - g * e)
        return A + B * self.v0

    def explosion_time(self, w):
        """T*(w) of Keller-Ressel (2011, eq. (6.2)): the moment E^q e^{wx} is finite exactly for τ < T*(w)."""
        w = np.asarray(w, dtype=float)
        disc = (self.kappa - self.rho * self.xi * w) ** 2 - self.xi**2 * (w * w - w)
        chi = self.rho * self.xi * w - self.kappa
        omega = np.sqrt(np.maximum(-disc, 0.0))
        with np.errstate(divide="ignore", invalid="ignore"):
            t = 2.0 * np.arctan2(omega, chi) / omega
        return np.where(disc < 0, t, np.inf)

    def strip(self, tau):
        return _heston_strip(self, tau)

    def cgf_derivatives(self, c, tau):
        lo, hi = self.strip(tau)
        radius = min(CAUCHY_RADIUS, 0.5 * (c - lo), 0.5 * (hi - c))
        return cauchy_derivatives(lambda z: self.cgf(z, tau), c, radius)

    def log_density_and_tail(self, x, tau, side):
        return fourier_log_density(self, x, tau), fourier_log_tail(self, x, tau, side)

    def log_distribution(self, x, tau):
        """log f, log F and log(1 − F) of x under Q^q by Fourier inversion."""
        return (fourier_log_density(self, x, tau), fourier_log_tail(self, x, tau, "lower"),
                fourier_log_tail(self, x, tau, "upper"))

    def engine(self, market: Market):
        _, quote, base, spot = market.quantlib()
        process = ql.HestonProcess(quote, base, spot, self.v0, self.kappa, self.theta, self.xi, self.rho)
        return ql.AnalyticHestonEngine(ql.HestonModel(process), ql.AnalyticHestonEngine.AngledContour,
                                       ql.AnalyticHestonEngine_Integration.expSinh(1e-12))


@functools.lru_cache(maxsize=None)
def _heston_strip(model: Heston, tau: float):
    """Critical moments (w₋, w₊): the roots of T*(w) = τ beyond the zeros of Δ(w), a concave quadratic."""
    a = -model.xi**2 * (1.0 - model.rho**2)
    b = model.xi**2 - 2.0 * model.kappa * model.rho * model.xi
    root = math.sqrt(b * b - 4.0 * a * model.kappa**2)
    zeros = sorted(((-b + root) / (2.0 * a), (-b - root) / (2.0 * a)))
    edges = []
    for w0, direction in ((zeros[0], -1.0), (zeros[1], 1.0)):
        # T* decreases from ∞ just beyond the zero of Δ to 0 as |w| → ∞.
        f = lambda w: float(model.explosion_time(w)) - tau
        near, step = w0 + direction * 1e-9 * max(1.0, abs(w0)), 1.0
        while f(w0 + direction * step) > 0:
            step *= 2.0
        edges.append(brentq(f, near, w0 + direction * step, xtol=1e-12, rtol=1e-14))
    return edges[0], edges[1]


def cauchy_derivatives(fun: Callable, centre, radius=CAUCHY_RADIUS, nodes=CAUCHY_NODES, orders=(1, 2, 3)):
    """Derivatives of an analytic function, real on the real axis, by the trapezoidal Cauchy integral."""
    theta = 2.0 * np.pi * np.arange(nodes) / nodes
    values = fun(centre + radius * np.exp(1j * theta))
    return np.array([math.factorial(n) * float(np.mean(values * np.exp(-1j * n * theta)).real) / radius**n
                     for n in orders])


def cumulants_to_raw(k1, k2, k3):
    return np.array([k1, k2 + k1**2, k3 + 3.0 * k1 * k2 + k1**3])


def reference_moments(model, tau, usd_base: bool):
    """E^USD y^k, k = 1, 2, 3: cumulants K_q^(n)(0) (EURUSD-type) or (−1)^n K_q^(n)(1) (USD-base)."""
    kappa = model.cgf_derivatives(1.0 if usd_base else 0.0, tau)
    if usd_base:
        kappa = kappa * np.array([-1.0, 1.0, -1.0])
    return cumulants_to_raw(*kappa)


def _admissible(model, tau, side):
    """Open interval of real c for the line of each integral: side 'put' (c < 0), 'call' (c > 1),
    'lower' (c < 0), 'upper' (c > 0) or 'any'."""
    lo, hi = model.strip(tau)
    lo = lo * (1.0 - 1e-12) if math.isfinite(lo) else lo
    hi = hi * (1.0 - 1e-12) if math.isfinite(hi) else hi
    return {"put": (lo, 0.0), "lower": (lo, 0.0), "call": (1.0, hi), "upper": (0.0, hi), "any": (lo, hi)}[side]


def _convex_argmin(psi, lo, hi, start):
    """Minimiser of a strictly convex psi on the open interval (lo, hi) that contains start.

    The bracket is expanded in doubling steps on each side until psi increases,
    then Brent's bounded method finds the minimum.
    """
    def safe(c):
        v = psi(c)
        return v if np.isfinite(v) else 1e300

    f0 = safe(start)

    def expand(direction):
        edge = hi if direction > 0 else lo
        a, fa, step = start, f0, 1.0
        while True:
            b = a + direction * step
            if not lo < b < hi:
                return a + (edge - a) * (1.0 - 1e-9)
            fb = safe(b)
            if fb > fa:
                return b
            a, fa, step = b, fb, 2.0 * step

    left, right = expand(-1.0), expand(1.0)
    res = minimize_scalar(safe, bounds=(left, right), method="bounded",
                          options={"xatol": 1e-8 * max(1.0, abs(left), abs(right))})
    return float(res.x) if res.fun <= f0 else start


def _line_integral(model, tau, k, c, kernel):
    """(1/π) ∫_0^∞ Re[e^{−wk + K_q(w) − s} κ(w)] du on w = c + iu, with s = −ck + K_q(c); returns (value, s)."""
    K_c = float(np.real(model.cgf(c, tau)))
    scale = -c * k + K_c
    lo, hi = model.strip(tau)
    h = min(1e-2 * max(1.0, abs(c)), 0.1 * (c - lo), 0.1 * (hi - c))
    curvature = (float(np.real(model.cgf(c + h, tau))) - 2.0 * K_c + float(np.real(model.cgf(c - h, tau)))) / h**2
    width = 1.0 / math.sqrt(curvature)  # the integrand is close to exp(−K_q''(c) u²/2) near u = 0
    k0 = abs(kernel(c))
    modulus = lambda u: math.exp(float(np.real(model.cgf(c + 1j * u, tau))) - K_c) * abs(kernel(c + 1j * u)) / k0
    top = 8.0 * width
    while modulus(top) > LINE_CUTOFF or modulus(1.5 * top) > LINE_CUTOFF:
        top *= 1.5
    integrand = lambda u: np.real(np.exp(-(c + 1j * u) * k + model.cgf(c + 1j * u, tau) - scale) * kernel(c + 1j * u))
    step = 0.5 * width
    u = np.arange(0.0, top + step, step)
    v = integrand(u)
    total = step * (v.sum() - 0.5 * v[0])
    for _ in range(MAX_HALVINGS):
        step *= 0.5
        new = 0.5 * total + step * integrand(np.arange(step, top + step, 2.0 * step)).sum()
        if abs(new - total) <= LINE_RTOL * abs(new):
            return new / np.pi, scale
        total = new
    raise RuntimeError("line integral did not converge")


def fourier_log_density(model, x, tau):
    """log f(x) under Q^q, integrating on the line through the saddle point K_q'(c) = x."""
    lo, hi = _admissible(model, tau, "any")
    out = []
    for xx in np.atleast_1d(np.asarray(x, dtype=float)):
        psi = lambda c: -c * xx + float(np.real(model.cgf(c, tau)))
        c = _convex_argmin(psi, lo, hi, 0.0)
        val, scale = _line_integral(model, tau, xx, c, lambda w: np.ones_like(w))
        out.append(math.log(val) + scale)
    return np.array(out)


def fourier_log_tail(model, x, tau, side):
    """log Q^q(x ≤ ·) (side 'lower') or log Q^q(x > ·) ('upper') by the tail formulas of Lee (2004b)."""
    sign = -1.0 if side == "lower" else 1.0
    lo, hi = _admissible(model, tau, side)
    out = []
    for xx in np.atleast_1d(np.asarray(x, dtype=float)):
        psi = lambda c: -c * xx + float(np.real(model.cgf(c, tau))) - math.log(abs(c))
        c = _convex_argmin(psi, lo, hi, sign * 1e-3)
        val, scale = _line_integral(model, tau, xx, c, lambda w: 1.0 / w)
        out.append(math.log(sign * val) + scale)
    return np.array(out)


def fourier_price(model, K, F, tau):
    """Undiscounted out-of-the-money prices (put below F, call at and above) by Fourier inversion."""
    out = []
    for kk in np.log(np.atleast_1d(np.asarray(K, dtype=float)) / F):
        side, start = ("put", -1.0) if kk < 0 else ("call", 2.0)
        lo, hi = _admissible(model, tau, side)
        psi = lambda c: -(c - 1.0) * kk + float(np.real(model.cgf(c, tau))) - math.log(abs(c * (c - 1.0)))
        c = _convex_argmin(psi, lo, hi, start)
        val, scale = _line_integral(model, tau, kk, c, lambda w: 1.0 / (w * (w - 1.0)))
        # e^{−(w−1)k} = e^{k} e^{−wk}: the line integral uses e^{−wk}, so rescale by e^{k}.
        out.append(F * val * math.exp(scale + kk))
    return np.array(out)


def quantlib_prices(model, market: Market, K, phi):
    """Undiscounted QuantLib prices NPV/D_q of calls (φ = 1) and puts (φ = −1)."""
    expiry = market.quantlib()[0]
    engine = model.engine(market)
    K = np.atleast_1d(np.asarray(K, dtype=float))
    phi = np.broadcast_to(phi, K.shape)
    out = np.empty(K.shape)
    for i, (k, p) in enumerate(zip(K, phi)):
        payoff = ql.PlainVanillaPayoff(ql.Option.Call if p == CALL else ql.Option.Put, float(k))
        option = ql.VanillaOption(payoff, ql.EuropeanExercise(expiry))
        option.setPricingEngine(engine)
        out[i] = option.NPV() / market.df_quote
    return out


def quantlib_otm_prices(model, market: Market, K):
    """Undiscounted out-of-the-money QuantLib prices: the put below F, the call at and above F."""
    K = np.atleast_1d(np.asarray(K, dtype=float))
    return quantlib_prices(model, market, K, np.where(K < market.forward, PUT, CALL))


def black_vols(F, K, tau, prices, lo=1e-4, hi=4.0, iterations=64):
    """Volatilities σ with forward_premium(F, K, σ, τ, φ) = price, φ = −1 below F, by bisection in ln σ."""
    K = np.asarray(K, dtype=float)
    phi = np.where(K < F, PUT, CALL)
    target = np.log(prices)
    a, b = np.full(K.shape, math.log(lo)), np.full(K.shape, math.log(hi))
    for _ in range(iterations):
        mid = 0.5 * (a + b)
        with np.errstate(divide="ignore"):
            above = np.log(forward_premium(F, K, np.exp(mid), tau, phi)) > target
        a, b = np.where(above, a, mid), np.where(above, mid, b)
    return np.exp(0.5 * (a + b))


class ModelSmile:
    """vol_of_strike from a model: a not-a-knot cubic spline of Black volatilities in ln(K/F), flat beyond."""

    def __init__(self, F, tau, x, vols, prices):
        self.F, self.tau = F, tau
        self.x, self.vols, self.prices = x, vols, prices
        self.spline = CubicSpline(x, vols, bc_type="not-a-knot")
        self.K_lo, self.K_hi = F * math.exp(x[0]), F * math.exp(x[-1])

    def __call__(self, K):
        return self.spline(np.clip(np.log(np.asarray(K, dtype=float) / self.F), self.x[0], self.x[-1]))


def model_smile(model, market: Market, spacing=SMILE_SPACING, reach=SMILE_REACH, floor=SMILE_FLOOR) -> ModelSmile:
    """Smile from QuantLib prices at nodes x_j = j·spacing, |x_j| ≤ reach, kept outwards while price > floor·F."""
    F, tau = market.forward, market.tau
    n = int(round(reach / spacing))
    x = spacing * np.arange(-n, n + 1)
    prices = quantlib_otm_prices(model, market, F * np.exp(x))
    ok = prices > floor * F
    centre = n
    lo = centre - int(np.argmin(np.append(ok[centre::-1], False)))
    hi = centre + int(np.argmin(np.append(ok[centre:], False)))
    x, prices = x[lo + 1:hi], prices[lo + 1:hi]
    return ModelSmile(F, tau, x, black_vols(F, F * np.exp(x), tau, prices), prices)


def hazard_rates(model, x, tau):
    """Reversed hazard rate r = f/F and hazard rate h = f/(1 − F) of x under Q^q."""
    lf, lF, lS = model.log_distribution(x, tau)
    return np.exp(lf - lF), np.exp(lf - lS)


def _tail_rate(model, x, tau, side):
    lf, lt = model.log_density_and_tail(x, tau, side)
    return np.exp(lf - lt)


def tail_exponent(model, tau, x_edge, side) -> dict:
    """inf of r over x ≤ x_edge (side 'lower') or of h over x ≥ x_edge ('upper'): the largest valid R4 exponent.

    Returns the infimum, where it is attained, the rate at the edge, the rate at
    the far end of the grid and whether the rate is nondecreasing away from the
    edge on the whole grid.
    """
    sign = -1.0 if side == "lower" else 1.0
    near = np.arange(0.0, TAIL_NEAR + 0.5 * TAIL_NEAR_STEP, TAIL_NEAR_STEP)
    far = np.geomspace(TAIL_NEAR, TAIL_FAR, TAIL_FAR_POINTS + 1)[1:]
    offsets = np.concatenate([near, far])
    rates = _tail_rate(model, x_edge + sign * offsets, tau, side)
    i = int(np.argmin(rates))
    best, where = float(rates[i]), float(offsets[i])
    lo, hi = offsets[max(i - 1, 0)], offsets[min(i + 1, len(offsets) - 1)]
    if hi > lo:
        res = minimize_scalar(lambda d: float(_tail_rate(model, x_edge + sign * d, tau, side)[0]),
                              bounds=(lo, hi), method="bounded", options={"xatol": 1e-9})
        if res.fun < best:
            best, where = float(res.fun), float(res.x)
    return {"inf": best, "at": x_edge + sign * where, "edge": float(rates[0]), "far": float(rates[-1]),
            "monotone": bool(np.all(np.diff(rates) >= 0.0))}


def price_bound(model, tau, F, K_edge, side):
    """Lee (2004a, Theorem 2.1) bound in R4 form: (B, e) with P(K) ≤ B (K/K_edge)^{e+1} on (0, K_edge]
    (side 'lower') or C(K) ≤ B (K/K_edge)^{1−e} on [K_edge, ∞) ('upper'), B minimised over the order."""
    lo, hi = model.strip(tau)
    ln_k = math.log(K_edge / F)
    if side == "lower":
        # ln[E S^{−q} q^q (1+q)^{−(1+q)} K^{1+q}] = K_q(−q) + q ln(K/F) + ln K + q ln q − (1+q) ln(1+q)
        psi = lambda q: (float(np.real(model.cgf(-q, tau))) + q * ln_k + math.log(K_edge)
                         + q * math.log(q) - (1.0 + q) * math.log1p(q))
        order = _convex_argmin(psi, 0.0, -lo, 1.0)
        return math.exp(psi(order)), order
    # ln[E S^{p+1} p^p (p+1)^{−(p+1)} K^{−p}] = K_q(p+1) − p ln(K/F) + ln F + p ln p − (p+1) ln(p+1)
    psi = lambda p: (float(np.real(model.cgf(p + 1.0, tau))) - p * ln_k + math.log(F)
                     + p * math.log(p) - (p + 1.0) * math.log1p(p))
    order = _convex_argmin(psi, 0.0, hi - 1.0, 1.0)
    return math.exp(psi(order)), order + 1.0
