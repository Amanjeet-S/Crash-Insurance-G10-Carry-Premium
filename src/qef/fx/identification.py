"""Sharp identification of option-implied moments from truncated quotes (R10).

Setting and notation are those of ``qef.fx.moments``: X = S_T > 0 under the
quote-currency forward measure Q^q, G(u) = Q^q(X ≤ u), Ḡ(u) = Q^q(X > u),
P(K) = ∫_0^K G and C(K) = ∫_K^∞ Ḡ, observed only on [K_min, K_max]. The raw
moments of y = ln(X_T/F_X) under Q^USD are m_k = mid_k + T_k^L + T_k^U, where
mid_k is the observed middle part and T_k^L = ∫_0^{K_min} w_k P,
T_k^U = ∫_{K_max}^∞ w_k C are the tail integrals of the contract weights
w_k = h_k'' of R3.

Observed at the boundaries are P0 = P(K_min), G0 = G(K_min) = P'(K_min+),
C0 = C(K_max) and Ḡ0 = Q(X ≥ K_max) = −C'(K_max−), the slopes from inside
the quoted range (theory/notes.tex, R10). R4 assumes
H_L(γ): G(u)/u^γ nondecreasing on (0, K_min], and H_U(η): u^η Ḡ(u)
nonincreasing on [K_max, ∞), η > 1.

Representation (lower tail). Under H_L(γ), G(K_min z) = G0 z^γ μ([0, z]) for
0 < z < 1 and a probability measure μ on [0, 1]: z^{−γ}G(K_min z)/G0 is
nondecreasing with values in [0, 1], so it is the distribution function of μ,
with the mass not reached below z = 1 placed at 1. Every such μ gives a law
compatible with the observations. By Tonelli's theorem a step at z_c
contributes

    a(z_c)   = (1 − z_c^{γ+1})/(γ+1)                                   to P0/(G0 K_min),
    b_k(z_c) = (K_min/(γ+1)) ∫_{z_c}^1 w_k(K_min z)(z^{γ+1} − z_c^{γ+1}) dz  to T_k^L/(G0 K_min),

so the identified set of (T_1^L, T_2^L, T_3^L) is
G0 K_min {∫ b dμ : ∫ a dμ = P0/(G0 K_min)}. Upper tail, u = K_max z, z ≥ 1:
Ḡ(K_max z) = Ḡ0 z^{−η} ν((z, ∞]) with ν a probability on [1, ∞],

    a(z_c)   = (1 − z_c^{1−η})/(η−1)                                   to C0/(Ḡ0 K_max),
    b_k(z_c) = (K_max/(η−1)) ∫_1^{z_c} w_k(K_max z)(z^{1−η} − z_c^{1−η}) dz  to T_k^U/(Ḡ0 K_max).

The two tails are constrained separately (the mean F and the total mass are
fixed by the observed prices), so the identified set of (m1, m2, m3) is the
convex set mid + S_L + S_U. Since a is monotone in z_c, a feasible measure
exists if and only if P0/(G0 K_min) ≤ a(0) = 1/(γ+1), that is
γ ≤ γ̄ = K_min G0/P0 − 1, and likewise η ≤ η̄ = 1 + K_max Ḡ0/C0; at the
maximal exponents μ = δ_0 and ν = δ_∞, a pure power tail, so the moments are
point-identified.

Computation. The measures are restricted to a grid of step points z_c (graded
towards the boundary), which gives an inner approximation that converges as the
grid is refined. The range of m1 is a linear programme; for a fixed mean t the
ranges of m2 and of the third central moment m3 − 3t m2 + 2t³ are linear
programmes, and their extremes over t are found on a grid of t with local
refinement. The variance is m2 − t² and the sign of the skewness is that of
the third central moment. The proofs are in the theory notes (R10).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import highspy
import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.optimize import minimize_scalar
from scipy.sparse import csc_matrix

from .moments import CONTRACTS, contract_weight, middle_contracts, smile_distribution

N_X = 20001  # nodes of the graded integration grid in log strike
N_C = 1000  # step points per tail; the inner approximation error is of order 1e-5 relative
GRADE = 3.0  # x = ∓L u^GRADE, dense near the boundary
L_MAX = 300.0  # largest |log(K/K_boundary)| integrated; beyond it the integrands are below e^{-L_MAX (exponent)}
THETA_MIN = 0.01  # the breakdown search does not go below this fraction of the maximal exponents


@dataclass(frozen=True)
class Boundary:
    """Observed boundary quantities and the exponents they imply."""

    K_min: float
    K_max: float
    P0: float
    G0: float
    C0: float
    Gbar0: float
    gamma_i: float  # setting (i): K g/G at K_min
    eta_i: float  # setting (i): K g/Ḡ at K_max

    @property
    def gamma_bar(self) -> float:
        """Largest γ compatible with (P0, G0): the strike elasticity of the put price at K_min, minus one."""
        return self.K_min * self.G0 / self.P0 - 1.0

    @property
    def eta_bar(self) -> float:
        """Largest η compatible with (C0, Ḡ0): one plus the absolute strike elasticity of the call price at K_max."""
        return 1.0 + self.K_max * self.Gbar0 / self.C0

    def exponents(self, theta: float) -> tuple[float, float]:
        """γ = θ γ̄ and η = 1 + θ(η̄ − 1), θ in (0, 1]."""
        return theta * self.gamma_bar, 1.0 + theta * (self.eta_bar - 1.0)


def boundary(F, tau, vol_of_strike, K_min, K_max) -> Boundary:
    lo = smile_distribution(K_min, F, tau, vol_of_strike)
    hi = smile_distribution(K_max, F, tau, vol_of_strike)
    return Boundary(K_min, K_max, lo["P"], lo["G"], hi["C"], hi["Gbar"],
                    K_min * lo["g"] / lo["G"], K_max * hi["g"] / hi["Gbar"])


def _step_points(n_x: int, n_c: int) -> np.ndarray:
    """Indices of the step points on the integration grid, spread evenly in the grid parameter and near the boundary."""
    even = np.round(np.linspace(0, n_x - 1, n_c)).astype(int)
    near = np.round(np.geomspace(1, n_x - 1, n_c)).astype(int)
    return np.unique(np.concatenate([even, near, [0, n_x - 1]]))


def lower_curve(F, K_min, gamma, usd_base, n_c=N_C, n_x=N_X):
    """a(z_c) and b_k(z_c), k = 1, 2, 3, on step points z_c in [0, 1] (the first is z_c = 0)."""
    L = min(L_MAX, max(40.0, 60.0 / gamma))  # contributions below K_min e^{-L} are of order e^{-γL}
    x = -L * np.linspace(1.0, 0.0, n_x) ** GRADE  # increasing from −L to 0, dense near 0
    z = np.exp(x)
    idx = n_x - 1 - _step_points(n_x, n_c)[::-1]
    zc = z[idx]
    A = np.concatenate([[1.0], 1.0 - zc ** (gamma + 1)]) / (gamma + 1)
    # Integrals from z to 1, accumulated from the boundary outwards on the reflected grid;
    # differencing cumulative integrals taken from the left end loses all precision, because
    # w(K_min z) z grows like 1/z as z -> 0.
    back = lambda f: cumulative_simpson(f[::-1], x=-x[::-1], initial=0.0)[::-1]
    B = []
    for k in CONTRACTS:
        w = contract_weight(k, K_min * z, F, usd_base) * K_min
        I1, I0 = back(w * z ** (gamma + 2)), back(w * z)  # dz = z dx
        B.append(np.concatenate([[I1[0]], I1[idx] - zc ** (gamma + 1) * I0[idx]]) / (gamma + 1))
    return A, np.array(B)


def upper_curve(F, K_max, eta, usd_base, n_c=N_C, n_x=N_X):
    """a(z_c) and b_k(z_c) on step points z_c in [1, ∞] (the last is z_c = ∞)."""
    L = min(L_MAX, max(40.0, 60.0 / (eta - 1.0)))
    x = L * np.linspace(0.0, 1.0, n_x) ** GRADE
    z = np.exp(x)
    idx = _step_points(n_x, n_c)
    zc = z[idx]
    A = np.concatenate([1.0 - zc ** (1 - eta), [1.0]]) / (eta - 1)
    B = []
    for k in CONTRACTS:
        w = contract_weight(k, K_max * z, F, usd_base) * K_max
        I1 = cumulative_simpson(w * z ** (2.0 - eta), x=x, initial=0.0)
        I0 = cumulative_simpson(w * z, x=x, initial=0.0)
        B.append(np.concatenate([I1[idx] - zc ** (1 - eta) * I0[idx], [I1[-1]]]) / (eta - 1))
    return A, np.array(B)


class SharpMomentSet:
    """The identified set mid + S_L + S_U of (m1, m2, m3) under H_L(γ) and H_U(η)."""

    def __init__(self, F, tau, vol_of_strike, usd_base, bd: Boundary, gamma, eta, n_c=N_C, mid=None,
                 put_beyond=None, call_beyond=None):
        """``put_beyond`` = (K, P(K)) with K < K_min and ``call_beyond`` = (K, C(K)) with K > K_max add
        observed prices beyond the boundaries as equality constraints (see tail_price_interval)."""
        if not (0 < gamma <= bd.gamma_bar and 1 < eta <= bd.eta_bar):
            raise ValueError("exponents outside the admissible range (0, γ̄] × (1, η̄]")
        self.mid = np.asarray(mid if mid is not None else
                              middle_contracts(F, bd.K_min, bd.K_max, tau, vol_of_strike, usd_base)["values"], dtype=float)
        AL, BL = lower_curve(F, bd.K_min, gamma, usd_base, n_c)
        AU, BU = upper_curve(F, bd.K_max, eta, usd_base, n_c)
        aL, aU = bd.P0 / (bd.G0 * bd.K_min), bd.C0 / (bd.Gbar0 * bd.K_max)
        self.M = np.hstack([bd.G0 * bd.K_min * BL, bd.Gbar0 * bd.K_max * BU])
        nL, n = len(AL), self.M.shape[1]
        rows = np.zeros((5, n))
        rows[0, :nL], rows[1, :nL], rows[2, nL:], rows[3, nL:] = 1.0, AL, 1.0, AU
        self._mean_scale = np.abs(self.M[0]).max()
        rows[4] = self.M[0] / self._mean_scale  # mean row, free unless a mean is imposed
        # At the maximal exponent the constraint forces the step at z_c = 0 (lower) or ∞ (upper);
        # clipping by a relative 1e-12 keeps the linear programme feasible in floating point.
        rhs = [1.0, min(aL, AL.max() * (1 - 1e-12)), 1.0, min(aU, AU.max() * (1 - 1e-12))]
        # An observed price beyond a boundary is linear in the mixing measure: a step at z_c
        # contributes max(0, a(z_c) − a_k) (lower tail, k = K/K_min) or max(0, a_U(z_c) − a_k)
        # (upper tail, k = K/K_max), in the units of the constraint on the boundary price.
        extra = []
        if put_beyond is not None:
            K, price = put_beyond
            row = np.zeros(n)
            row[:nL] = np.maximum(0.0, AL - _a_lower(K / bd.K_min, gamma))
            extra.append((row, price / (bd.G0 * bd.K_min)))
        if call_beyond is not None:
            K, price = call_beyond
            row = np.zeros(n)
            row[nL:] = np.maximum(0.0, AU - _a_upper(K / bd.K_max, eta))
            extra.append((row, price / (bd.Gbar0 * bd.K_max)))
        if extra:
            rows = np.vstack([rows] + [r for r, _ in extra])
        rhs = np.array(rhs)
        # One HiGHS model per set; each problem changes only the costs and the mean row, and
        # the simplex solver starts from the previous basis (Huangfu and Hall, 2018).
        lp = highspy.HighsLp()
        lp.num_col_, lp.num_row_ = n, rows.shape[0]
        lp.col_cost_, lp.col_lower_, lp.col_upper_ = np.zeros(n), np.zeros(n), np.full(n, highspy.kHighsInf)
        extra_rhs = np.array([v for _, v in extra])
        lp.row_lower_ = np.concatenate([rhs, [-highspy.kHighsInf], extra_rhs])
        lp.row_upper_ = np.concatenate([rhs, [highspy.kHighsInf], extra_rhs])
        A = csc_matrix(rows)
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        lp.a_matrix_.start_, lp.a_matrix_.index_, lp.a_matrix_.value_ = A.indptr, A.indices, A.data
        self._h = highspy.Highs()
        self._h.setOptionValue("output_flag", False)
        self._h.setOptionValue("threads", 1)  # callers parallelise across currency-months
        self._h.passModel(lp)
        self._cols = np.arange(n, dtype=np.int32)

    def _lp(self, cost, mean=None):
        """Minimum of cost · (μ, ν) over the feasible measures, with the mean fixed at ``mean`` if given."""
        h = self._h
        h.changeColsCost(len(self._cols), self._cols, np.asarray(cost, dtype=float))
        if mean is None:
            h.changeRowBounds(4, -highspy.kHighsInf, highspy.kHighsInf)
        else:
            v = (mean - self.mid[0]) / self._mean_scale
            h.changeRowBounds(4, v, v)
        h.run()
        return h.getInfo().objective_function_value if h.getModelStatus() == highspy.HighsModelStatus.kOptimal else np.nan

    def m1_range(self):
        c = self.M[0]
        return self.mid[0] + self._lp(c), self.mid[0] - self._lp(-c)

    def third_central(self, t, sign):
        """sign = +1: minimum, −1: maximum of m3 − 3t m2 + 2t³ among feasible laws with mean t."""
        c = self.M[2] - 3.0 * t * self.M[1]
        base = self.mid[2] - 3.0 * t * self.mid[1] + 2.0 * t ** 3
        return base + sign * self._lp(sign * c, mean=t)

    def second_central(self, t, sign):
        return self.mid[1] + sign * self._lp(sign * self.M[1], mean=t) - t * t

    def _extreme(self, f, sign, n_t, n_refine=2):
        """Extreme over the feasible means t: a grid of n_t means, refined around the n_refine best."""
        if not hasattr(self, "_m1"):
            self._m1 = self.m1_range()
        lo, hi = self._m1
        if not hi - lo > 1e-12 * max(abs(lo), abs(hi), 1e-300):
            return f(0.5 * (lo + hi), sign)
        eps = 1e-7 * (hi - lo)
        ts = np.linspace(lo + eps, hi - eps, n_t)
        vals = np.array([sign * f(t, sign) for t in ts])
        best = np.nanmin(vals)
        for i in np.argsort(np.where(np.isnan(vals), np.inf, vals))[:n_refine]:
            a, b = ts[max(i - 1, 0)], ts[min(i + 1, n_t - 1)]
            r = minimize_scalar(lambda t: sign * f(t, sign), bounds=(a, b), method="bounded",
                                options={"xatol": 1e-6 * (hi - lo)})
            if np.isfinite(r.fun):
                best = min(best, r.fun)
        return sign * best

    def variance_range(self, n_t=25):
        return self._extreme(self.second_central, 1.0, n_t), self._extreme(self.second_central, -1.0, n_t)

    def third_central_range(self, n_t=25):
        return self._extreme(self.third_central, 1.0, n_t), self._extreme(self.third_central, -1.0, n_t)


def _a_lower(k, gamma):
    """a at the step z_c = k in the lower tail: (1 − k^{γ+1})/(γ+1)."""
    return (1.0 - k ** (gamma + 1.0)) / (gamma + 1.0)


def _a_upper(k, eta):
    """a_U at the step z_c = k in the upper tail: (1 − k^{1−η})/(η − 1), tending to ln k as η → 1."""
    return math.log(k) if eta == 1.0 else (1.0 - k ** (1.0 - eta)) / (eta - 1.0)


def tail_price_interval(bd: Boundary, K, exponent, side):
    """Identified interval of the put price at K < K_min under H_L(γ) (side 'lower') or of the call
    price at K > K_max under H_U(η) (side 'upper'), given the boundary price and slope.

    A step at z_c contributes max(0, a(z_c) − a_k) to the price, a convex function of a, so the
    lower end is the step itself and the upper end the chord to the extreme step (R10(b)):
    lower tail [max(0, P0 − G0 K_min a_k), P0 k^{γ+1}], upper tail [max(0, C0 − Ḡ0 K_max a_k),
    C0 k^{1−η}]. The upper ends are R4's bounds. The exponent may be 0 (lower) or 1 (upper),
    the limits of the weakest hypotheses of the class.
    """
    if side == "lower":
        k = K / bd.K_min
        if not 0 < k < 1 or not 0 <= exponent <= bd.gamma_bar:
            raise ValueError("need K < K_min and 0 ≤ γ ≤ γ̄")
        return max(0.0, bd.P0 - bd.G0 * bd.K_min * _a_lower(k, exponent)), bd.P0 * k ** (exponent + 1.0)
    k = K / bd.K_max
    if not k > 1 or not 1 <= exponent <= bd.eta_bar:
        raise ValueError("need K > K_max and 1 ≤ η ≤ η̄")
    return max(0.0, bd.C0 - bd.Gbar0 * bd.K_max * _a_upper(k, exponent)), bd.C0 * k ** (1.0 - exponent)


def max_consistent_exponent(bd: Boundary, K, price, side, xtol=1e-10):
    """Largest exponent of the class at which an observed price beyond the boundary lies in the
    identified interval, and whether the price is consistent with the class at all.

    The upper end of the interval decreases and the lower end increases with the exponent, so the
    consistent exponents form an interval starting at the weakest hypothesis (γ → 0 or η → 1).
    Returns (exponent, θ, consistent), θ = γ/γ̄ or (η − 1)/(η̄ − 1); (nan, nan, False) if the price
    lies outside the interval even for the weakest hypothesis.
    """
    lo_e, hi_e = (0.0, bd.gamma_bar) if side == "lower" else (1.0, bd.eta_bar)
    inside = lambda e: (lambda iv: iv[0] <= price <= iv[1])(tail_price_interval(bd, K, e, side))
    if not inside(lo_e):
        return math.nan, math.nan, False
    if inside(hi_e):
        e = hi_e
    else:
        a, b = lo_e, hi_e
        while b - a > xtol * max(1.0, b):
            m = 0.5 * (a + b)
            a, b = (m, b) if inside(m) else (a, m)
        e = a
    theta = e / bd.gamma_bar if side == "lower" else (e - 1.0) / (bd.eta_bar - 1.0)
    return e, theta, True


def pareto_moments(F, tau, vol_of_strike, usd_base, bd: Boundary, mid=None) -> np.ndarray:
    """Raw moments at the maximal exponents (γ̄, η̄), where the tails are pure power laws."""
    mid = np.asarray(mid if mid is not None else
                     middle_contracts(F, bd.K_min, bd.K_max, tau, vol_of_strike, usd_base)["values"], dtype=float)
    _, BL = lower_curve(F, bd.K_min, bd.gamma_bar, usd_base, n_c=2)
    _, BU = upper_curve(F, bd.K_max, bd.eta_bar, usd_base, n_c=2)
    return mid + bd.G0 * bd.K_min * BL[:, 0] + bd.Gbar0 * bd.K_max * BU[:, -1]


def skewness_sign_breakdown(F, tau, vol_of_strike, usd_base, bd: Boundary, tol=1e-2, n_c=N_C, n_t=15,
                            theta_min=THETA_MIN):
    """Smallest θ at which the sign of the skewness is identified, and that sign.

    At θ = 1 the moments are point-identified; their third central moment gives
    the only sign that can be identified. The identified sets shrink as θ rises
    (R10), so the set of θ at which the third-central-moment range excludes zero
    is an interval ending at 1, and bisection finds its left end to within ``tol``.
    Returns (θ*, sign); θ* = ``theta_min`` if the sign is already identified
    there. If the Pareto third central moment is exactly zero, R10(e) has no
    breakdown point and the sign is identified at no fraction; the function then
    returns (1.0, 0) as a convention, and the sign 0 marks the case, which does
    not occur in the sample (research log, 3 October 2026).
    """
    mid = middle_contracts(F, bd.K_min, bd.K_max, tau, vol_of_strike, usd_base)["values"]
    m = pareto_moments(F, tau, vol_of_strike, usd_base, bd, mid)
    k3 = m[2] - 3 * m[0] * m[1] + 2 * m[0] ** 3
    if k3 == 0:
        return 1.0, 0
    sign = -1 if k3 < 0 else 1

    def identified(theta):
        S = SharpMomentSet(F, tau, vol_of_strike, usd_base, bd, *bd.exponents(theta), n_c=n_c, mid=mid)
        # sign −1 needs the maximum below zero; sign +1 needs the minimum above zero
        ext = S._extreme(S.third_central, 1.0 if sign > 0 else -1.0, n_t)
        return ext > 0 if sign > 0 else ext < 0

    lo, hi = theta_min, 1.0
    if identified(theta_min):
        return theta_min, sign
    while hi - lo > tol:
        mid_t = 0.5 * (lo + hi)
        lo, hi = (lo, mid_t) if identified(mid_t) else (mid_t, hi)
    return hi, sign
