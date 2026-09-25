import itertools

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import norm

from qef.fx.conventions import G10
from qef.fx.gk import CALL, PUT, d_plus_minus, forward_premium, strike_from_delta_smile
from qef.fx.moments import (
    CONTRACTS,
    _simpson,
    boundary_elasticities,
    closed_form_bounds,
    closed_form_check,
    contract_payoff,
    contract_weight,
    implied_moment_intervals,
    middle_contracts,
    moment_ranges,
    otm_price,
    smile_distribution,
    tail_bounds,
    variance_skewness,
    weight_sign_changes,
)
from qef.fx.sabr import sabr_vol

TAU = 1 / 12
PAIRS = [("EUR", 1.10), ("JPY", 150.0)]  # an EURUSD-type pair and a USD-base pair


def flat(sig):
    return lambda K: np.full_like(np.asarray(K, dtype=float), sig)


def normal_raw_moments(mu, s):
    return np.array([mu, mu**2 + s**2, mu**3 + 3 * mu * s**2])


def ten_delta_strikes(ccy, F, vol):
    conv = G10[ccy].delta
    return (strike_from_delta_smile(-0.10, F, TAU, PUT, conv, vol, 0.998),
            strike_from_delta_smile(0.10, F, TAU, CALL, conv, vol, 0.998))


def true_tail_parts(k, F, usd_base, side, K_edge, vol):
    """Tail integrals of the positive and negative parts of h_k'' against the smile's prices."""
    f = lambda x: float(contract_weight(k, K_edge * np.exp(x), F, usd_base) * otm_price(K_edge * np.exp(x), F, TAU, vol)) * K_edge * np.exp(x)
    cuts = [np.log(b / K_edge) for b in weight_sign_changes(k, F, usd_base)]
    if side == "lower":
        edges = [-40.0, *sorted(c for c in cuts if c < 0), 0.0]
    else:
        edges = [0.0, *sorted(c for c in cuts if c > 0), 40.0]
    pos = neg = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        v = quad(f, lo, hi, epsabs=0.0, epsrel=1e-12, limit=500)[0]
        pos, neg = (pos + v, neg) if v >= 0 else (pos, neg - v)
    return pos, neg


@pytest.mark.parametrize("usd_base", [False, True])
@pytest.mark.parametrize("k", CONTRACTS)
def test_contract_second_derivative_matches_finite_differences(k, usd_base):
    F = 150.0 if usd_base else 1.10
    for K in F * np.array([0.3, 0.8, 0.97, 1.02, 1.3, 4.0]):
        h = 1e-4 * K
        fd = (contract_payoff(k, K + h, F, usd_base) - 2 * contract_payoff(k, K, F, usd_base) + contract_payoff(k, K - h, F, usd_base)) / h**2
        assert contract_weight(k, K, F, usd_base) == pytest.approx(fd, rel=1e-6, abs=1e-12 / F**2)


@pytest.mark.parametrize("usd_base", [False, True])
def test_weights_change_sign_only_where_stated(usd_base):
    F = 150.0 if usd_base else 1.10
    K = F * np.exp(np.linspace(-4, 4, 4000))  # no node at F, F e^{±1} or F e^{±2}
    for k in CONTRACTS:
        w = contract_weight(k, K, F, usd_base)
        flips = K[1:][np.sign(w[1:]) != np.sign(w[:-1])]
        expected = list(weight_sign_changes(k, F, usd_base)) + ([F] if k == 3 else [])
        assert len(flips) == len(expected)
        for b in expected:
            assert np.min(np.abs(np.log(flips / b))) < 3e-3


@pytest.mark.parametrize("a", [-1.0, 1.0, 2.0])
def test_usd_law_is_the_tilted_lognormal(a):
    # E^USD (S_T/F)^a = E^q[(S_T/F)^{a+1}] under N(−s²/2, s²) equals E e^{a z} under z ~ N(s²/2, s²).
    s = 0.2
    dens_q = lambda z: norm.pdf(z, -s**2 / 2, s)
    lhs = quad(lambda z: np.exp((a + 1) * z) * dens_q(z), -12 * s, 12 * s, epsabs=0, epsrel=1e-12)[0]
    rhs = quad(lambda z: np.exp(a * z) * norm.pdf(z, s**2 / 2, s), -12 * s, 12 * s, epsabs=0, epsrel=1e-12)[0]
    assert lhs == pytest.approx(rhs, rel=1e-10)
    assert rhs == pytest.approx(np.exp(a * (a + 1) * s**2 / 2), rel=1e-10)  # a = −1: E^USD(1/S_T) = 1/F


@pytest.mark.parametrize("ccy,F", PAIRS)
def test_flat_smile_untruncated_reproduces_lognormal_moments(ccy, F):
    sig = 0.10
    s = sig * np.sqrt(TAU)
    usd_base = G10[ccy].usd_base
    wide = (F * np.exp(-15 * s), F * np.exp(15 * s))  # prices beyond 15 s are below 1e-45
    # Q^USD: y ~ N(−s²/2, s²) for both pair types.
    r = middle_contracts(F, *wide, TAU, flat(sig), usd_base)
    assert r["converged"]
    np.testing.assert_allclose(r["values"], normal_raw_moments(-s**2 / 2, s), rtol=1e-9)
    # Q^q: ln(S_T/F) ~ N(−s²/2, s²); for a USD-base pair y = −ln(S_T/F) ~ N(+s²/2, s²) under Q^q.
    q = middle_contracts(F, *wide, TAU, flat(sig), False)["values"]
    np.testing.assert_allclose(q, normal_raw_moments(-s**2 / 2, s), rtol=1e-9)
    sign = np.array([-1.0, 1.0, -1.0]) if usd_base else np.ones(3)
    np.testing.assert_allclose(sign * q, normal_raw_moments((1 if usd_base else -1) * s**2 / 2, s), rtol=1e-9)
    v, sk = variance_skewness(*r["values"])
    assert v == pytest.approx(s**2, rel=1e-9) and abs(sk) < 1e-6


def test_change_of_numeraire_identity():
    # E^USD(1/S_T) = E^q[(1/S_T)(S_T/F)] = 1/F, with the Q^q density g = C'' of a skewed smile;
    # E^q(1/S_T) exceeds 1/F (Jensen, R1 remark).
    F, a, rho, nu = 150.0, 0.10, 0.3, 1.0
    vol = lambda K: sabr_vol(K, F, TAU, a, rho, nu, 1.0)
    s = a * np.sqrt(TAU)
    K = np.linspace(F * np.exp(-10 * s), F * np.exp(10 * s), 8001)
    q = np.array([smile_distribution(k, F, TAU, vol)["g"] for k in K])
    h = K[1] - K[0]
    assert q.min() > 0
    assert _simpson(K * q, h) == pytest.approx(F, rel=1e-7)
    assert _simpson((1 / K) * (K / F) * q, h) == pytest.approx(1 / F, rel=1e-7)
    assert _simpson(q / K, h) > 1 / F * (1 + 1e-4)
    # Spanned USD moments agree with integration against the density (K/F) g.
    r = middle_contracts(F, K[0], K[-1], TAU, vol, True)["values"]
    y = -np.log(K / F)
    np.testing.assert_allclose(r, [_simpson(y**k * (K / F) * q, h) for k in CONTRACTS], rtol=1e-4)


def test_usd_moments_equal_spanning_with_options_on_the_inverted_rate():
    # Options on X = 1/S in USD: P_X(k) = (k/F) C(1/k) and C_X(k) = (k/F) P(1/k) by R1(b).
    # Spanning y = ln(X/F_X) with f_k under Q^USD must equal spanning g_k under Q^q.
    F, a, rho, nu = 150.0, 0.10, 0.3, 1.5
    vol = lambda K: sabr_vol(K, F, TAU, a, rho, nu, 1.0)
    K_min, K_max = ten_delta_strikes("JPY", F, vol)
    direct = middle_contracts(F, K_min, K_max, TAU, vol, True)["values"]
    FX, n, via_x = 1 / F, 8192, np.zeros(3)
    for lo, hi in ((1 / K_max, FX), (FX, 1 / K_min)):
        k = np.linspace(lo, hi, n + 1)
        price_x = (k / F) * otm_price(1 / k, F, TAU, vol)
        via_x += [_simpson(contract_weight(j, k, FX, False) * price_x, (hi - lo) / n) for j in CONTRACTS]
    np.testing.assert_allclose(direct, via_x, rtol=1e-10)


def test_smile_distribution_matches_finite_differences():
    F = 0.65
    vol = lambda K: sabr_vol(K, F, TAU, 0.12, -0.5, 2.5, 1.0)
    P = lambda x: float(forward_premium(F, x, float(vol(x)), TAU, PUT))
    C = lambda x: float(forward_premium(F, x, float(vol(x)), TAU, CALL))
    d1 = lambda f, K, h: (f(K + h) - f(K - h)) / (2 * h)
    d2 = lambda f, K, h: (f(K + h) - 2 * f(K) + f(K - h)) / h**2
    rich = lambda d, f, K, h: (4 * d(f, K, h / 2) - d(f, K, h)) / 3  # Richardson: O(h⁴) error
    for K in (0.60, 0.64, 0.66, 0.70):
        d = smile_distribution(K, F, TAU, vol)
        h = 2e-4 * K
        assert d["G"] == pytest.approx(rich(d1, P, K, h), rel=1e-8)
        assert d["Gbar"] == pytest.approx(-rich(d1, C, K, h), rel=1e-8)
        assert d["g"] == pytest.approx(rich(d2, P, K, h), rel=1e-6)


def test_flat_smile_boundary_elasticities_are_lognormal():
    F, sig = 1.10, 0.10
    s = sig * np.sqrt(TAU)
    K_min, K_max = ten_delta_strikes("EUR", F, flat(sig))
    gamma, eta = boundary_elasticities(F, K_min, K_max, TAU, flat(sig))
    d_lo, d_hi = float(d_plus_minus(F, K_min, sig, TAU)[1]), float(d_plus_minus(F, K_max, sig, TAU)[1])
    assert gamma == pytest.approx(norm.pdf(d_lo) / (s * norm.cdf(-d_lo)), rel=1e-12)
    assert eta == pytest.approx(norm.pdf(d_hi) / (s * norm.cdf(d_hi)), rel=1e-12)


@pytest.mark.parametrize("ccy,F", PAIRS)
def test_bounds_contain_true_lognormal_tails(ccy, F):
    vol = flat(0.10)
    usd_base = G10[ccy].usd_base
    K_min, K_max = ten_delta_strikes(ccy, F, vol)
    gamma, eta = boundary_elasticities(F, K_min, K_max, TAU, vol)
    P_min, C_max = float(otm_price(K_min, F, TAU, vol)), float(otm_price(K_max, F, TAU, vol))
    for k in CONTRACTS:
        w = lambda K, k=k: contract_weight(k, K, F, usd_base)
        breaks = weight_sign_changes(k, F, usd_base)
        for side, K_e, p_e, e in (("lower", K_min, P_min, gamma), ("upper", K_max, C_max, eta)):
            b_pos, b_neg, _ = tail_bounds(w, side, K_e, p_e, e, breaks)
            t_pos, t_neg = true_tail_parts(k, F, usd_base, side, K_e, vol)
            assert 0 <= t_pos <= b_pos * (1 + 1e-10) and 0 <= t_neg <= b_neg * (1 + 1e-10)


@pytest.mark.parametrize("usd_base", [False, True])
def test_bounds_decrease_in_the_exponent(usd_base):
    F = 150.0 if usd_base else 1.10
    K_min, K_max, P_min, C_max = 0.95 * F, 1.05 * F, 1e-3 * F, 1e-3 * F
    for k in CONTRACTS:
        w = lambda K, k=k: contract_weight(k, K, F, usd_base)
        breaks = weight_sign_changes(k, F, usd_base)
        lower = np.array([tail_bounds(w, "lower", K_min, P_min, g, breaks)[:2] for g in (0.5, 1, 2, 5, 20, 60)])
        upper = np.array([tail_bounds(w, "upper", K_max, C_max, e, breaks)[:2] for e in (1.2, 1.5, 2, 5, 20, 60)])
        for b in (*lower.T, *upper.T):
            assert np.all(b >= 0)
            assert np.all(np.diff(b) < 0) or np.all(b == 0)


@pytest.mark.parametrize("ccy,F", PAIRS)
def test_truncated_lognormal_intervals_contain_true_moments(ccy, F):
    sig = 0.10
    s = sig * np.sqrt(TAU)
    out = implied_moment_intervals(F, TAU, flat(sig), G10[ccy].usd_base, G10[ccy].delta, 0.998)
    true = normal_raw_moments(-s**2 / 2, s)
    assert out["simpson_converged"] and out["simpson_rel_change"] < 1e-10
    for name in ("i", "ii"):
        assert out[f"status_{name}"] == "ok"
        for j, k in enumerate(CONTRACTS):
            assert out[f"m{k}_lo_{name}"] <= true[j] <= out[f"m{k}_hi_{name}"]
        assert out[f"var_lo_{name}"] <= s**2 <= out[f"var_hi_{name}"]
        assert out[f"skew_lo_{name}"] <= 0 <= out[f"skew_hi_{name}"]
        assert out[f"closed_form_rel_diff_{name}"] < 1e-10
    for key in ("var", "skew"):  # heavier tails give wider intervals
        assert out[f"{key}_lo_ii"] <= out[f"{key}_lo_i"] and out[f"{key}_hi_i"] <= out[f"{key}_hi_ii"]


@pytest.mark.parametrize("gamma,eta", [(0.7, 1.3), (2.0, 2.0), (12.0, 30.0), (60.0, 60.0)])
@pytest.mark.parametrize("F,K_min,K_max", [(1.10, 1.05, 1.16), (150.0, 120.0, 400.0)])
def test_numeric_tail_bounds_equal_closed_forms(gamma, eta, F, K_min, K_max):
    P_min, C_max = 0.003 * F, 0.002 * F
    cf = closed_form_bounds(F, K_min, P_min, gamma, K_max, C_max, eta)
    inv_sq = lambda K: K**-2.0
    assert tail_bounds(inv_sq, "lower", K_min, P_min, gamma)[0] == pytest.approx(cf["lower_inverse_square"], rel=1e-10)
    assert tail_bounds(inv_sq, "upper", K_max, C_max, eta)[0] == pytest.approx(cf["upper_inverse_square"], rel=1e-10)
    lower_log = tail_bounds(lambda K: 2 * (1 + np.log(F / K)) / K**2, "lower", K_min, P_min, gamma)
    assert lower_log[0] == pytest.approx(cf["lower_log"], rel=1e-10) and lower_log[1] == 0
    upper_log = tail_bounds(lambda K: 2 * (1 - np.log(K / F)) / K**2, "upper", K_max, C_max, eta, (np.e * F,))
    assert upper_log[0] + upper_log[1] == pytest.approx(cf["upper_log"], rel=1e-10)
    assert closed_form_check(F, K_min, P_min, gamma, K_max, C_max, eta) < 1e-10


def test_upper_log_closed_form_only_applies_up_to_eF():
    # The R4 closed form for the upper logarithmic weight is stated for F < K_max ≤ eF.
    F, K_min, P_min, C_max = 1.0, 0.9, 1e-3, 1e-3
    assert np.isfinite(closed_form_bounds(F, K_min, P_min, 2.0, 2.5, C_max, 2.0)["upper_log"])
    cf = closed_form_bounds(F, K_min, P_min, 2.0, 3.0, C_max, 2.0)
    assert np.isnan(cf["upper_log"]) and np.isfinite(cf["upper_inverse_square"])
    assert closed_form_check(F, K_min, P_min, 2.0, 3.0, C_max, 2.0) < 1e-10


def test_upper_tail_needs_eta_above_one():
    with pytest.raises(ValueError):
        tail_bounds(lambda K: 1 / K, "upper", 1.2, 0.01, 1.0)
    with pytest.raises(ValueError):
        tail_bounds(lambda K: 1 / K, "lower", 1.0, 0.01, 0.0)


@pytest.mark.parametrize("lo,hi,off_vertex", [
    ([-4.5e-4, 7.5e-4, -3e-5], [-3.0e-4, 1.1e-3, 1e-5], None),
    ([-2e-3, 6e-4, -4e-3], [-4e-4, 4e-3, 1e-3], None),
    ([0.08, 0.015, 0.0015], [0.10, 0.025, 0.0020], None),  # m2 = m3/m1 inside an edge but not extremal
    ([0.0756, 0.0091, 0.00082], [0.0789, 0.0131, 0.0009], "skew_lo"),  # minimum at the edge point m2 = m3/m1
    ([-0.0978, 0.0151, -0.0021], [-0.0896, 0.0276, -0.00199], "skew_hi"),  # maximum at the edge point m2 = m3/m1
    ([0.05, 0.0045, 3.5e-4], [0.06, 0.006, 4.0e-4], "skew_lo"),  # minimum at the edge point m1 = m2²/m3
    ([-1e-3, 5e-4, -1e-5], [1e-3, 9e-4, 1e-5], "var_hi"),  # m1 = 0 inside: variance maximum on an edge
])
def test_moment_ranges_match_grid_search(lo, hi, off_vertex):
    lo, hi = np.array(lo), np.array(hi)
    r = moment_ranges(lo, hi)
    g = [np.linspace(a, b, 61) for a, b in zip(lo, hi)]
    m1, m2, m3 = np.meshgrid(*g, indexing="ij")
    v, sk = variance_skewness(m1, m2, m3)
    assert r["range_gap"] < 1e-12
    for key, grid in (("var", v), ("skew", sk)):
        span = grid.max() - grid.min()
        assert r[f"{key}_lo"] <= grid.min() + 1e-12 * span and r[f"{key}_hi"] >= grid.max() - 1e-12 * span
        assert r[f"{key}_lo"] >= grid.min() - 1e-3 * span and r[f"{key}_hi"] <= grid.max() + 1e-3 * span
    if off_vertex:  # the extremum lies off the vertices, so the edge candidates are needed
        vertices = np.array(list(itertools.product(*zip(lo, hi))))
        at_vertices = variance_skewness(*vertices.T)[0 if off_vertex.startswith("var") else 1]
        span = at_vertices.max() - at_vertices.min()
        gain = at_vertices.min() - r[off_vertex] if off_vertex.endswith("lo") else r[off_vertex] - at_vertices.max()
        assert gain > 1e-3 * span


def test_nonpositive_variance_gives_unbounded_skewness():
    r = moment_ranges(np.array([-0.03, 5e-4, -1e-5]), np.array([0.0, 9e-4, 1e-5]))
    assert r["var_lo"] < 0 and r["skew_lo"] == -np.inf and r["skew_hi"] == np.inf


def test_simpson_halving_check_is_recorded():
    F = 0.65
    vol = lambda K: sabr_vol(K, F, TAU, 0.12, -0.5, 2.5, 1.0)
    K_min, K_max = ten_delta_strikes("AUD", F, vol)
    r = middle_contracts(F, K_min, K_max, TAU, vol, False, n0=8)
    assert r["converged"] and r["rel_change"] < 1e-10 and r["n"] % 16 == 0
    coarse = middle_contracts(F, K_min, K_max, TAU, vol, False, n0=8, rtol=1e-4)
    assert coarse["n"] < r["n"]
