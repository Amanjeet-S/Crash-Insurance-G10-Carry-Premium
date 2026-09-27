import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm

from qef.fx.conventions import G10
from qef.fx.gk import CALL, PUT, strike_from_delta_smile
from qef.fx.identification import (
    SharpMomentSet,
    boundary,
    lower_curve,
    pareto_moments,
    skewness_sign_breakdown,
    upper_curve,
)
from qef.fx.moments import contract_intervals, contract_weight, middle_contracts, moment_ranges, otm_price, tail_bounds
from qef.fx.sabr import sabr_vol

TAU = 1 / 12
F_EUR = 1.10
FAST = dict(n_c=300)


def flat(sig):
    return lambda K: np.full_like(np.asarray(K, dtype=float), sig)


def strikes(ccy, F, vol, delta=0.10):
    conv = G10[ccy].delta
    return (strike_from_delta_smile(-delta, F, TAU, PUT, conv, vol, 0.998),
            strike_from_delta_smile(delta, F, TAU, CALL, conv, vol, 0.998))


def lp_range(A, b, a):
    from scipy.optimize import linprog
    Aeq, beq = np.vstack([np.ones(len(A)), A]), [1.0, a]
    lo = linprog(b, A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs").fun
    hi = -linprog(-b, A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs").fun
    return lo, hi


@pytest.mark.parametrize("delta", [0.05, 0.10, 0.25])
def test_lognormal_satisfies_setting_i_below_the_maximal_exponent(delta):
    # ln X normal: log-concave distribution and survival functions, so the boundary
    # elasticities are admissible (γ_i ≤ γ̄, η_i ≤ η̄) at every boundary.
    vol = flat(0.10)
    bd = boundary(F_EUR, TAU, vol, *strikes("EUR", F_EUR, vol, delta))
    assert 0 < bd.gamma_i <= bd.gamma_bar and 1 < bd.eta_i <= bd.eta_bar
    # γ̄ + 1 is the strike elasticity of the put price at K_min
    h = 1e-6 * bd.K_min
    P = lambda K: float(otm_price(K, F_EUR, TAU, vol))
    elasticity = bd.K_min * (P(bd.K_min + h) - P(bd.K_min - h)) / (2 * h) / P(bd.K_min)
    assert bd.gamma_bar + 1 == pytest.approx(elasticity, rel=1e-7)


def test_strongly_skewed_smile_can_refute_setting_i():
    F = 0.65
    vol = lambda K: sabr_vol(K, F, TAU, 0.12, -0.5, 2.5, 1.0)
    bd = boundary(F, TAU, vol, *strikes("AUD", F, vol))
    assert bd.gamma_i > bd.gamma_bar  # H_L(γ_i) is incompatible with (P0, G0)
    with pytest.raises(ValueError):
        SharpMomentSet(F, TAU, vol, False, bd, bd.gamma_i, 2.0, **FAST)


def test_maximal_exponent_gives_the_pareto_tail():
    # At γ̄ the law below K_min is G0 (u/K_min)^γ̄, so P(K) = P0 (K/K_min)^{γ̄+1} there.
    vol = flat(0.10)
    bd = boundary(F_EUR, TAU, vol, *strikes("EUR", F_EUR, vol))
    g = bd.gamma_bar
    assert bd.G0 * bd.K_min / (g + 1) == pytest.approx(bd.P0, rel=1e-12)
    _, BL = lower_curve(F_EUR, bd.K_min, g, False, n_c=2)
    _, BU = upper_curve(F_EUR, bd.K_max, bd.eta_bar, False, n_c=2)
    for k in (1, 2, 3):
        w = lambda K: contract_weight(k, K, F_EUR, False)
        lo = quad(lambda K: w(K) * bd.P0 * (K / bd.K_min) ** (g + 1), 0, bd.K_min, epsabs=0, epsrel=1e-12, limit=400)[0]
        e = bd.eta_bar
        hi = quad(lambda K: w(K) * bd.C0 * (K / bd.K_max) ** (1 - e), bd.K_max, np.inf, epsabs=0, epsrel=1e-12, limit=400)[0]
        assert bd.G0 * bd.K_min * BL[k - 1, 0] == pytest.approx(lo, rel=1e-8)
        assert bd.Gbar0 * bd.K_max * BU[k - 1, -1] == pytest.approx(hi, rel=1e-8)
    m = pareto_moments(F_EUR, TAU, vol, False, bd)
    assert np.all(np.isfinite(m))


@pytest.mark.parametrize("theta", [0.05, 0.4, 0.8])
def test_single_integral_bounds_and_sharpness_of_r4(theta):
    # For a nonnegative weight the sharp upper bound is R4's bound and the sharp lower bound
    # is the single step at z* with a(z*) = P0/(G0 K_min), not zero.
    vol = flat(0.10)
    bd = boundary(F_EUR, TAU, vol, *strikes("EUR", F_EUR, vol))
    g, e = bd.exponents(theta)
    A, B = lower_curve(F_EUR, bd.K_min, g, False, n_c=1000)
    a = bd.P0 / (bd.G0 * bd.K_min)
    k = 2  # w_2 = 2(1 − y)/K² > 0 in the lower tail
    lo, hi = (bd.G0 * bd.K_min * v for v in lp_range(A, B[k - 1], a))
    r4 = tail_bounds(lambda K: contract_weight(k, K, F_EUR, False), "lower", bd.K_min, bd.P0, g)[0]
    assert hi == pytest.approx(r4, rel=1e-6)
    z_star = brentq(lambda z: (1 - z ** (g + 1)) / (g + 1) - a, 0.0, 1.0, xtol=1e-15)
    w = lambda z: contract_weight(k, bd.K_min * z, F_EUR, False) * bd.K_min
    step = quad(lambda z: w(z) * (z ** (g + 1) - z_star ** (g + 1)) / (g + 1), z_star, 1, epsabs=0, epsrel=1e-12)[0]
    assert lo == pytest.approx(bd.G0 * bd.K_min * step, rel=5e-4) and lo > 0  # inner approximation on the step grid


def test_joint_set_contains_the_lognormal_truth_and_is_sharp():
    vol = flat(0.10)
    s2 = 0.01 * TAU
    bd = boundary(F_EUR, TAU, vol, *strikes("EUR", F_EUR, vol))
    theta_i = min(bd.gamma_i / bd.gamma_bar, (bd.eta_i - 1) / (bd.eta_bar - 1))
    mid = middle_contracts(F_EUR, bd.K_min, bd.K_max, TAU, vol, False)["values"]
    prev = None
    for theta in (0.2, 0.5, theta_i):
        S = SharpMomentSet(F_EUR, TAU, vol, False, bd, *bd.exponents(theta), mid=mid, **FAST)
        v, k3 = S.variance_range(), S.third_central_range()
        assert v[0] <= s2 <= v[1] and k3[0] <= 0.0 <= k3[1]
        # inside R4's box, and nested in θ
        lo, hi, *_ = contract_intervals(mid, F_EUR, bd.K_min, bd.P0, bd.K_max, bd.C0, False, *bd.exponents(theta))
        box = moment_ranges(lo, hi)
        assert box["var_lo"] <= v[0] * (1 + 1e-9) and v[1] <= box["var_hi"] * (1 + 1e-9)
        if prev is not None:
            assert prev[0][0] <= v[0] * (1 + 1e-9) and v[1] <= prev[0][1] * (1 + 1e-9)
            assert prev[1][0] <= k3[0] + 1e-15 and k3[1] <= prev[1][1] + 1e-15
        prev = (v, k3)
    # Beyond the lognormal's own tail exponent the hypothesis is false and the truth can fall outside.
    S = SharpMomentSet(F_EUR, TAU, vol, False, bd, *bd.exponents(0.97), mid=mid, **FAST)
    assert S.variance_range()[0] > s2


def test_step_grid_refinement_changes_little():
    F = 0.65
    vol = lambda K: sabr_vol(K, F, TAU, 0.12, -0.5, 2.5, 1.0)
    bd = boundary(F, TAU, vol, *strikes("AUD", F, vol))
    ex = bd.exponents(0.3)
    coarse = SharpMomentSet(F, TAU, vol, False, bd, *ex, n_c=300).third_central_range()
    fine = SharpMomentSet(F, TAU, vol, False, bd, *ex, n_c=1000).third_central_range()
    width = fine[1] - fine[0]
    assert abs(coarse[0] - fine[0]) < 1e-3 * width and abs(coarse[1] - fine[1]) < 1e-3 * width


def test_skewness_sign_breakdown_is_consistent():
    F = 0.65
    vol = lambda K: sabr_vol(K, F, TAU, 0.12, -0.5, 2.5, 1.0)
    bd = boundary(F, TAU, vol, *strikes("AUD", F, vol))
    theta, sign = skewness_sign_breakdown(F, TAU, vol, False, bd, tol=0.01, **FAST)
    assert sign == -1 and 0.01 < theta < 1
    mid = middle_contracts(F, bd.K_min, bd.K_max, TAU, vol, False)["values"]
    above = SharpMomentSet(F, TAU, vol, False, bd, *bd.exponents(min(theta + 0.01, 1.0)), mid=mid, **FAST)
    below = SharpMomentSet(F, TAU, vol, False, bd, *bd.exponents(theta - 0.02), mid=mid, **FAST)
    assert above.third_central_range()[1] < 0 <= below.third_central_range()[1]


def test_tail_price_interval_matches_the_linear_programme_and_r4():
    # A price beyond the boundary is linear in the mixing measure; its identified interval is closed-form.
    from qef.fx.identification import tail_price_interval
    vol = flat(0.10)
    bd = boundary(F_EUR, TAU, vol, *strikes("EUR", F_EUR, vol))
    K_lo, K_hi = 0.97 * bd.K_min, 1.03 * bd.K_max
    for theta in (0.1, 0.5, 0.9):
        g, e = bd.exponents(theta)
        A, _ = lower_curve(F_EUR, bd.K_min, g, False, n_c=1000)
        k = K_lo / bd.K_min
        row = np.maximum(0.0, A - (1 - k ** (g + 1)) / (g + 1))
        lo, hi = (bd.G0 * bd.K_min * v for v in lp_range(A, row, bd.P0 / (bd.G0 * bd.K_min)))
        L, U = tail_price_interval(bd, K_lo, g, "lower")
        assert (lo, hi) == (pytest.approx(L, rel=1e-4, abs=1e-12 * bd.P0), pytest.approx(U, rel=1e-6))  # linprog feasibility tolerance
        assert U == pytest.approx(bd.P0 * k ** (g + 1), rel=1e-14)  # R4 lemma (a)
        AU, _ = upper_curve(F_EUR, bd.K_max, e, False, n_c=1000)
        k = K_hi / bd.K_max
        row = np.maximum(0.0, AU - (1 - k ** (1 - e)) / (e - 1))
        lo, hi = (bd.Gbar0 * bd.K_max * v for v in lp_range(AU, row, bd.C0 / (bd.Gbar0 * bd.K_max)))
        L, U = tail_price_interval(bd, K_hi, e, "upper")
        assert (lo, hi) == (pytest.approx(L, rel=1e-4, abs=1e-12 * bd.C0), pytest.approx(U, rel=1e-6))  # linprog feasibility tolerance
    # At the maximal exponent the interval is the Pareto price.
    L, U = tail_price_interval(bd, K_lo, bd.gamma_bar, "lower")
    assert L == pytest.approx(U, rel=1e-12)


def test_consistent_exponents_of_an_observed_price():
    from qef.fx.identification import max_consistent_exponent
    vol = flat(0.10)
    bd = boundary(F_EUR, TAU, vol, *strikes("EUR", F_EUR, vol))
    K5p, K5c = strikes("EUR", F_EUR, vol, 0.05)
    p5, c5 = float(otm_price(K5p, F_EUR, TAU, vol)), float(otm_price(K5c, F_EUR, TAU, vol))
    # The lognormal satisfies H_L(γ_i) and H_U(η_i), so its 5Δ prices are consistent at least up to those.
    g, th, ok = max_consistent_exponent(bd, K5p, p5, "lower")
    assert ok and g >= bd.gamma_i * (1 - 1e-6) and th == pytest.approx(g / bd.gamma_bar)
    e, th, ok = max_consistent_exponent(bd, K5c, c5, "upper")
    assert ok and e >= bd.eta_i * (1 - 1e-6)
    # A put price above the convexity chord P0 K/K_min is inconsistent with every hypothesis of the class.
    assert max_consistent_exponent(bd, K5p, 1.01 * bd.P0 * K5p / bd.K_min, "lower")[2] is False


def test_observed_prices_beyond_the_boundary_shrink_the_identified_set():
    vol = flat(0.10)
    s2 = 0.01 * TAU
    bd = boundary(F_EUR, TAU, vol, *strikes("EUR", F_EUR, vol))
    K5p, K5c = strikes("EUR", F_EUR, vol, 0.05)
    p5, c5 = float(otm_price(K5p, F_EUR, TAU, vol)), float(otm_price(K5c, F_EUR, TAU, vol))
    mid = middle_contracts(F_EUR, bd.K_min, bd.K_max, TAU, vol, False)["values"]
    ex = bd.exponents(0.5)
    free = SharpMomentSet(F_EUR, TAU, vol, False, bd, *ex, mid=mid, **FAST)
    tied = SharpMomentSet(F_EUR, TAU, vol, False, bd, *ex, mid=mid, put_beyond=(K5p, p5), call_beyond=(K5c, c5), **FAST)
    v0, v1 = free.variance_range(), tied.variance_range()
    k0, k1 = free.third_central_range(), tied.third_central_range()
    assert v0[0] <= v1[0] <= s2 <= v1[1] <= v0[1] and k0[0] <= k1[0] <= 0.0 <= k1[1] <= k0[1]
    assert (v1[1] - v1[0]) < 0.9 * (v0[1] - v0[0])
