"""Tests for the vanna–volga smile (qef.fx.vannavolga).

Equation numbers refer to Castagna and Mercurio (2007), authors' preprint.
QuantLib 1.43 exposes vanna–volga in its Python bindings only as barrier-option
engines (VannaVolgaBarrierEngine, VannaVolgaIKDoubleBarrierEngine,
VannaVolgaWODoubleBarrierEngine), not as a vanilla smile, so there is no
directly comparable function and no QuantLib cross-check here.
"""

import numpy as np
import pytest
from scipy.stats import norm

from qef.fx.gk import CALL, PUT, DeltaConvention, atm_dns_strike, d_plus_minus, delta, forward_premium, pa_call_delta_maximiser, vega
from qef.fx.smile import SabrSmile, SmileQuotes, market_strangle_strikes, quotes_from_smile
from qef.fx.vannavolga import (FIRST_ORDER, METHODS, PRICE, SECOND_ORDER, VannaVolgaSmile, calibrate_vv, gk_implied_vol,
                               implied_quotes, smile_delta_strike, vv_smile)

PIPS = DeltaConvention(spot=True, premium_adjusted=False)
PA = DeltaConvention(spot=True, premium_adjusted=True)

# One-month quotes of G10 magnitude: F, tau, ATM, 25Δ RR, 25Δ BF, base-currency DF, convention.
QUOTES = {
    "eurusd_like": SmileQuotes(1.10, 1 / 12, 0.09, -0.012, 0.0025, 0.25, 0.998, PIPS),
    "audusd_like": SmileQuotes(0.70, 1 / 12, 0.12, -0.025, 0.004, 0.25, 0.997, PIPS),
    "usdjpy_like": SmileQuotes(110.0, 1 / 12, 0.10, -0.02, 0.0035, 0.25, 1.002, PA),
    "usdcad_like": SmileQuotes(1.30, 1 / 12, 0.07, 0.008, 0.002, 0.25, 0.999, PA),
}
# Skewed quotes where the price is negative on the low-volatility call wing.
SKEWED = SmileQuotes(110.0, 1 / 12, 0.10, -0.04, 0.004, 0.25, 0.998, PA)


def _y(sm, K):
    """Coefficients of eq. (11), written out independently of the module."""
    K1, K2, K3 = sm.strikes
    y1 = np.log(K2 / K) * np.log(K3 / K) / (np.log(K2 / K1) * np.log(K3 / K1))
    y2 = np.log(K / K1) * np.log(K3 / K) / (np.log(K2 / K1) * np.log(K3 / K2))
    y3 = np.log(K / K1) * np.log(K / K2) / (np.log(K3 / K1) * np.log(K3 / K2))
    return y1, y2, y3


def _greeks(S0, rd, rf, T, sigma, K):
    """Vega, volga and vanna of a call under GK, eq. (5)."""
    d1 = (np.log(S0 / K) + (rd - rf + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    V = S0 * np.exp(-rf * T) * np.sqrt(T) * norm.pdf(d1)
    return np.array([V, V * d1 * d2 / sigma, -V * d2 / (S0 * sigma * np.sqrt(T))])


def _market():
    S0, rd, rf, T = 1.10, 0.03, 0.01, 0.25
    F = S0 * np.exp((rd - rf) * T)
    s2 = 0.095
    sm = VannaVolgaSmile(F, T, (1.04, F * np.exp(0.5 * s2**2 * T), 1.17), (0.112, s2, 0.090))
    return S0, rd, rf, T, F, sm


def test_weights_solve_hedging_system():
    # System (4): the weighted pillar calls match the vega, volga and vanna of the call at K.
    S0, rd, rf, T, F, sm = _market()
    K = F * np.exp(np.linspace(-0.3, 0.3, 41))
    x = sm.weights(K)
    pillar = _greeks(S0, rd, rf, T, sm.vols[1], np.array(sm.strikes))  # (greek, pillar)
    np.testing.assert_allclose(pillar @ x, _greeks(S0, rd, rf, T, sm.vols[1], K), rtol=1e-10, atol=1e-14)


def test_weights_are_kronecker_at_pillars():
    *_, sm = _market()
    np.testing.assert_allclose(sm.weights(np.array(sm.strikes)), np.eye(3), atol=1e-15)


def test_price_equals_strike_independent_greek_cost_form():
    # Remark 4.1: C(K) = C^BS(K) + y·(vega, volga, vanna)(K) with y = (A')^{-1} c, on discounted prices.
    S0, rd, rf, T, F, sm = _market()
    D = np.exp(-rd * T)
    Ki = np.array(sm.strikes)
    c = D * (forward_premium(F, Ki, np.array(sm.vols), T, CALL) - forward_premium(F, Ki, sm.vols[1], T, CALL))
    y = np.linalg.solve(_greeks(S0, rd, rf, T, sm.vols[1], Ki).T, c)
    K = F * np.exp(np.linspace(-0.25, 0.25, 31))
    expected = D * forward_premium(F, K, sm.vols[1], T, CALL) + y @ _greeks(S0, rd, rf, T, sm.vols[1], K)
    np.testing.assert_allclose(D * sm.vv_premium(K, CALL), expected, rtol=1e-10)


def test_put_call_parity_and_boundary_limits():
    *_, F, sm = _market()
    K = F * np.exp(np.linspace(-0.4, 0.4, 41))
    np.testing.assert_allclose(sm.vv_premium(K, CALL) - sm.vv_premium(K, PUT), F - K, atol=1e-14)
    # Boundary limits (preprint p. 6), undiscounted: C → F as K → 0 and C → 0 as K → ∞.
    assert sm.vv_premium(F * 1e-3, CALL) == pytest.approx(F - F * 1e-3, rel=1e-12)
    assert abs(sm.vv_premium(F * 3.0, CALL)) < 1e-12


def test_price_is_invariant_to_the_pillar_set():
    # Prop. 6.1: pillars moved to other strikes, priced by the smile itself, give the same price. The
    # Greeks' volatility is held at σ2 by keeping K2 as the middle pillar.
    *_, F, sm = _market()
    L = (0.98, sm.strikes[1], 1.21)
    moved = VannaVolgaSmile(F, sm.tau, L, tuple(sm.vol(np.array(L))))
    K = F * np.exp(np.linspace(-0.4, 0.4, 41))
    np.testing.assert_allclose(moved.vv_premium(K), sm.vv_premium(K), rtol=0, atol=1e-14)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("reading", ["market", "smile"])
@pytest.mark.parametrize("name", QUOTES)
def test_pillar_reproduction(method, reading, name):
    cal = calibrate_vv(QUOTES[name], reading, method)
    assert cal.success, cal.status
    sm = cal.smile
    np.testing.assert_allclose(sm.vol(np.array(sm.strikes)), sm.vols, rtol=0, atol=1e-12)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("conv", [PIPS, PA])
def test_flat_pillars_give_flat_smile(method, conv):
    q = SmileQuotes(1.10, 0.25, 0.10, 0.0, 0.0, 0.25, 0.995, conv)
    sm = vv_smile(q, 0.0, method)
    sd = 0.10 * np.sqrt(0.25)
    K = 1.10 * np.exp(np.linspace(-5 * sd, 5 * sd, 201))
    np.testing.assert_allclose(sm.vol(K), 0.10, rtol=0, atol=1e-12)
    np.testing.assert_allclose(sm.vv_premium(K), forward_premium(1.10, K, 0.10, 0.25, CALL), rtol=0, atol=1e-15)


@pytest.mark.parametrize("name", QUOTES)
def test_price_vol_reprices_vv_premium(name):
    q = QUOTES[name]
    sm = vv_smile(q, q.bf, PRICE)
    sd = q.sigma_atm * np.sqrt(q.tau)
    K = q.F * np.exp(np.linspace(-2.5 * sd, 2.5 * sd, 101))
    phi = np.where(K >= q.F, CALL, PUT)
    v = sm.vol(K)
    assert np.all(np.isfinite(v))
    np.testing.assert_allclose(forward_premium(q.F, K, v, q.tau, phi), sm.vv_premium(K, phi), rtol=1e-11)


@pytest.mark.parametrize("name", QUOTES)
def test_approximations_close_to_price_between_pillars(name):
    # Castagna and Mercurio (2007, p. 10) describe both approximations as accurate inside [K1, K3].
    # Tolerances set here for one-month G10-size quotes: 0.05 vol points for ς1 and 0.01 for ς2;
    # that ς2 is the closer of the two is checked, not taken from the source.
    q = QUOTES[name]
    sm = vv_smile(q, q.bf, PRICE)
    K = np.linspace(sm.strikes[0], sm.strikes[2], 201)
    exact = sm.vol(K)
    e1 = np.max(np.abs(sm.vol(K, FIRST_ORDER) - exact))
    e2 = np.max(np.abs(sm.vol(K, SECOND_ORDER) - exact))
    assert e1 < 5e-4 and e2 < 1e-4 and e2 < e1


def test_second_order_matches_eq12_away_from_zero_denominator():
    q = QUOTES["audusd_like"]
    sm = vv_smile(q, q.bf, SECOND_ORDER)
    s1, s2, s3 = sm.vols
    sd = q.sigma_atm * np.sqrt(q.tau)
    K = q.F * np.exp(np.linspace(-3 * sd, 3 * sd, 301))
    y1, y2, y3 = _y(sm, K)
    d1, d2 = d_plus_minus(q.F, K, s2, q.tau)
    dd = [np.prod(d_plus_minus(q.F, k, s2, q.tau)) for k in sm.strikes]
    D1 = y1 * s1 + y2 * s2 + y3 * s3 - s2
    D2 = y1 * dd[0] * (s1 - s2) ** 2 + y3 * dd[2] * (s3 - s2) ** 2
    rad = s2**2 + d1 * d2 * (2 * s2 * D1 + D2)
    keep = (np.abs(d1 * d2) > 1e-3) & (rad >= 0)
    direct = s2 + (-s2 + np.sqrt(rad[keep])) / (d1 * d2)[keep]
    np.testing.assert_allclose(sm.vol(K[keep]), direct, rtol=1e-12)
    assert keep.sum() > 250


@pytest.mark.parametrize("conv", [PIPS, PA])
def test_second_order_limit_where_d_plus_d_minus_vanishes(conv):
    q = SmileQuotes(1.10, 0.25, 0.10, -0.02, 0.004, 0.25, 0.995, conv)
    sm = vv_smile(q, q.bf, SECOND_ORDER)
    s1, s2, s3 = sm.vols
    K2 = sm.strikes[1]
    assert sm.vol(K2) == pytest.approx(s2, abs=1e-15)  # d+ d− = 0 at the DNS strike in either convention
    # The other zero of d+ d−: d− = 0 (pips) or d+ = 0 (premium-adjusted).
    Kb = q.F * np.exp((0.5 if conv.premium_adjusted else -0.5) * s2**2 * q.tau)
    assert abs(np.prod(d_plus_minus(q.F, Kb, s2, q.tau))) < 1e-12
    y1, y2, y3 = _y(sm, Kb)
    dd = [np.prod(d_plus_minus(q.F, k, s2, q.tau)) for k in sm.strikes]
    D1 = y1 * s1 + y2 * s2 + y3 * s3 - s2
    D2 = y1 * dd[0] * (s1 - s2) ** 2 + y3 * dd[2] * (s3 - s2) ** 2
    assert sm.vol(Kb) == pytest.approx(s2 + D1 + D2 / (2 * s2), abs=1e-14)
    # Continuity: the symmetric average of neighbouring values agrees to second order in the step.
    near = sm.vol(Kb * np.exp(np.array([-1e-6, 1e-6])))
    assert near.mean() == pytest.approx(sm.vol(Kb), abs=1e-10)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("name", QUOTES)
def test_market_reading_reprices_strangle_and_matches_quote_definitions(method, name):
    q = QUOTES[name]
    cal = calibrate_vv(q, "market", method)
    assert cal.success, cal.status
    sm = cal.smile
    # Market strangle: the smile prices it at the market-strangle strikes as the flat volatility does.
    kc, kp = market_strangle_strikes(q)
    assert (cal.strikes["K_call_ms"], cal.strikes["K_put_ms"]) == (kc, kp)
    s = q.sigma_atm + q.bf
    flat = forward_premium(q.F, kc, s, q.tau, CALL) + forward_premium(q.F, kp, s, q.tau, PUT)
    value = forward_premium(q.F, kc, sm.vol(kc), q.tau, CALL) + forward_premium(q.F, kp, sm.vol(kp), q.tau, PUT)
    assert abs(value - flat) / (vega(q.F, kc, s, q.tau, 1.0) + vega(q.F, kp, s, q.tau, 1.0)) < 1e-10
    # Risk reversal at the smile's own ±25Δ strikes, which are the pillars K3 and K1.
    k25c = smile_delta_strike(sm, 0.25, CALL, q.convention, q.df_base)
    k25p = smile_delta_strike(sm, -0.25, PUT, q.convention, q.df_base)
    assert k25c == pytest.approx(sm.strikes[2], rel=1e-10) and k25p == pytest.approx(sm.strikes[0], rel=1e-10)
    assert sm.vol(k25c) - sm.vol(k25p) == pytest.approx(q.rr, abs=1e-10)
    # ATM: K2 is a delta-neutral-straddle fixed point of the smile.
    assert atm_dns_strike(q.F, sm.vol(sm.strikes[1]), q.tau, q.convention) == pytest.approx(sm.strikes[1], rel=1e-13)
    assert sm.vol(sm.strikes[1]) == pytest.approx(q.sigma_atm, abs=1e-12)


@pytest.mark.parametrize("method", METHODS)
def test_smile_reading_uses_quoted_butterfly(method):
    q = QUOTES["usdjpy_like"]
    cal = calibrate_vv(q, "smile", method)
    assert cal.success and cal.bf_ss == q.bf
    v1, _, v3 = cal.smile.vol(np.array(cal.smile.strikes))
    assert 0.5 * (v1 + v3) - q.sigma_atm == pytest.approx(q.bf, abs=1e-12)


@pytest.mark.parametrize("reading", ["market", "smile"])
@pytest.mark.parametrize("conv", [PIPS, PA])
def test_implied_quotes_recover_quotes_generated_by_another_smile(reading, conv):
    # Quotes are generated by a SABR smile; the calibrated vanna–volga smile must imply them back
    # under the definitions of quotes_from_smile.
    F, tau, db = 1.1, 1 / 12, 0.998
    atm, rr, bf = quotes_from_smile(SabrSmile(F, tau, 0.09, -0.3, 1.5), tau, conv, 0.25, db, reading)
    q = SmileQuotes(F, tau, atm, rr, bf, 0.25, db, conv)
    cal = calibrate_vv(q, reading)
    assert cal.success, cal.status
    iq = implied_quotes(cal.smile, conv, 0.25, db, reading)
    np.testing.assert_allclose([iq["atm"], iq["rr"], iq["bf"]], [atm, rr, bf], rtol=0, atol=1e-9)


def test_market_and_smile_readings_differ_for_skewed_quotes():
    cal = calibrate_vv(QUOTES["usdjpy_like"], "market")
    assert cal.success and cal.bf_ss - QUOTES["usdjpy_like"].bf > 1e-4


def test_premium_adjusted_pillars():
    q = QUOTES["usdjpy_like"]
    sm = calibrate_vv(q, "smile").smile
    (K1, K2, K3), (s1, s2, s3) = sm.strikes, sm.vols
    assert K2 == pytest.approx(q.F * np.exp(-0.5 * q.sigma_atm**2 * q.tau), rel=1e-14)
    assert delta(q.F, K1, s1, q.tau, PUT, PA, q.df_base) == pytest.approx(-0.25, abs=1e-12)
    assert delta(q.F, K3, s3, q.tau, CALL, PA, q.df_base) == pytest.approx(0.25, abs=1e-12)
    # Conventional branch: right of the delta maximiser, where the premium-adjusted call delta decreases.
    assert K3 > pa_call_delta_maximiser(q.F, s3, q.tau)
    assert delta(q.F, K3 * 1.001, s3, q.tau, CALL, PA, q.df_base) < delta(q.F, K3 * 0.999, s3, q.tau, CALL, PA, q.df_base)


@pytest.mark.parametrize("method", METHODS)
def test_vol_accepts_scalars_and_arrays(method):
    sm = vv_smile(QUOTES["eurusd_like"], 0.0025, method)
    assert np.ndim(sm.vol(1.1)) == 0 and np.isfinite(sm.vol(1.1))
    K = np.linspace(1.05, 1.15, 6)
    assert sm.vol(K).shape == (6,)
    np.testing.assert_array_equal(sm.vol(K.reshape(2, 3)), sm.vol(K).reshape(2, 3))


@pytest.mark.parametrize("method", METHODS)
def test_invalid_strikes_give_nan(method):
    sm = vv_smile(QUOTES["eurusd_like"], 0.0025, method)
    for k in (0.0, -1.0, np.nan, np.inf):
        assert np.isnan(sm.vol(k))
    out = sm.vol(np.array([1.1, 0.0, np.nan]))
    assert np.isfinite(out[0]) and np.all(np.isnan(out[1:]))
    assert np.isnan(sm.vv_premium(-1.0))


def test_price_is_nan_where_the_premium_leaves_the_bounds():
    sm = vv_smile(SKEWED, SKEWED.bf, PRICE)
    sd = SKEWED.sigma_atm * np.sqrt(SKEWED.tau)
    K = SKEWED.F * np.exp(np.linspace(0.1 * sd, 4 * sd, 200))
    p = sm.vv_premium(K, CALL)
    v = sm.vol(K)
    assert np.any(p <= 0)  # negative vanna–volga call price on the low-volatility wing
    assert np.all(np.isnan(v[p <= 0])) and np.all(np.isfinite(v[p > 1e-10 * SKEWED.F]))


def test_second_order_is_nan_where_the_radicand_is_negative():
    sm = vv_smile(SKEWED, SKEWED.bf, SECOND_ORDER)
    s1, s2, s3 = sm.vols
    sd = SKEWED.sigma_atm * np.sqrt(SKEWED.tau)
    K = SKEWED.F * np.exp(np.linspace(-4 * sd, 4 * sd, 400))
    y1, y2, y3 = _y(sm, K)
    d1, d2 = d_plus_minus(SKEWED.F, K, s2, SKEWED.tau)
    dd = [np.prod(d_plus_minus(SKEWED.F, k, s2, SKEWED.tau)) for k in sm.strikes]
    rad = s2**2 + d1 * d2 * (2 * s2 * (y1 * s1 + y2 * s2 + y3 * s3 - s2) + y1 * dd[0] * (s1 - s2) ** 2 + y3 * dd[2] * (s3 - s2) ** 2)
    assert np.any(rad < 0)
    v = sm.vol(K)
    assert np.all(np.isnan(v[rad < 0]))
    assert np.all(np.isfinite(v[(rad > 0) & (np.abs(np.log(K / SKEWED.F)) < sd)]))


def test_first_order_is_nan_where_negative():
    # Concave pillars: the parabola in log strike of eq. (11) turns negative in the wings.
    F, tau = 1.1, 0.25
    sm = VannaVolgaSmile(F, tau, (1.0, F, 1.2), (0.05, 0.10, 0.05), FIRST_ORDER)
    K = F * np.exp(np.linspace(-0.6, 0.6, 121))
    y1, y2, y3 = _y(sm, K)
    s = 0.05 * y1 + 0.10 * y2 + 0.05 * y3
    assert np.any(s <= 0)
    v = sm.vol(K)
    assert np.all(np.isnan(v[s <= 0]))
    np.testing.assert_allclose(v[s > 0], s[s > 0], rtol=1e-13)


def test_calibration_returns_status_instead_of_raising():
    bad = SmileQuotes(1.1, 1 / 12, 0.09, -0.01, -0.2, 0.25, 0.998, PIPS)  # negative pillar and strangle volatilities
    for method in METHODS:
        assert calibrate_vv(bad, "smile", method).status == "pillars_undefined"
        cal = calibrate_vv(bad, "market", method)
        assert cal.status == "strangle_undefined"
        assert cal.smile is None and not cal.success and np.all(np.isnan(cal.residuals))


def test_delta_strike_is_nan_when_no_strike_attains_the_delta():
    # The first-order smile's premium-adjusted put delta stays above 0.10 in absolute value, so no
    # 10Δ put strike exists (the case of Reiswich, 2010, Fig. 3.6).
    sm = calibrate_vv(SKEWED, "market", FIRST_ORDER).smile
    assert np.isnan(smile_delta_strike(sm, -0.10, PUT, PA, SKEWED.df_base))
    iq = implied_quotes(sm, PA, 0.10, SKEWED.df_base, "market")
    assert np.isnan(iq["rr"]) and np.isnan(iq["bf"]) and np.isnan(iq["K_put"])


def test_delta_strike_found_before_an_undefined_region():
    # The flat-volatility 10Δ call strike lies where the price is undefined, so the search of
    # strike_from_delta_smile, which starts there, fails; the fallback finds the root between K2 and
    # the undefined region.
    q = SmileQuotes(110.0, 1 / 12, 0.10, -0.05, 0.004, 0.25, 0.998, PIPS)
    sm = calibrate_vv(q, "market", PRICE).smile
    assert np.isnan(smile_delta_strike(sm, 0.10, CALL, PIPS, q.df_base, fallback=False))
    k = smile_delta_strike(sm, 0.10, CALL, PIPS, q.df_base)
    assert sm.strikes[2] < k
    assert delta(q.F, k, sm.vol(k), q.tau, CALL, PIPS, q.df_base) == pytest.approx(0.10, abs=1e-10)
    assert np.all(np.isfinite(sm.vol(np.linspace(sm.strikes[1], k, 50))))
    assert np.any(np.isnan(sm.vol(np.linspace(k, 1.1 * k, 50))))


@pytest.mark.parametrize("phi", [CALL, PUT])
def test_gk_implied_vol_round_trip_and_bounds(phi):
    F, tau = 1.1, 0.25
    x = np.linspace(0.0, 0.3, 13)
    K = F * np.exp(phi * x)  # out of the money, as the smile inverts
    sig = np.linspace(0.05, 0.4, 13)
    p = forward_premium(F, K, sig, tau, phi)
    np.testing.assert_allclose(gk_implied_vol(p, F, K, tau, phi), sig, rtol=1e-12)
    np.testing.assert_allclose(gk_implied_vol(p[:3], F, K[:3], tau, phi), sig[:3], rtol=1e-12)  # scalar path
    K_itm = F * np.exp(-phi * np.array([0.02, 0.05]))
    p_itm = forward_premium(F, K_itm, 0.12, tau, phi)
    np.testing.assert_allclose(gk_implied_vol(p_itm, F, K_itm, tau, phi), 0.12, rtol=1e-10)
    intrinsic = np.maximum(phi * (F - K_itm), 0.0)
    assert np.all(np.isnan(gk_implied_vol(intrinsic, F, K_itm, tau, phi)))
    upper = F if phi == CALL else K
    assert np.all(np.isnan(gk_implied_vol(upper + 0 * K, F, K, tau, phi)))
    assert np.isnan(gk_implied_vol(-1e-3, F, F, tau, phi))


def test_invalid_pillars_raise():
    with pytest.raises(ValueError):
        VannaVolgaSmile(1.1, 0.25, (1.1, 1.05, 1.2), (0.1, 0.1, 0.1))
    with pytest.raises(ValueError):
        VannaVolgaSmile(1.1, 0.25, (1.0, 1.1, 1.2), (0.1, -0.1, 0.1))
    with pytest.raises(ValueError):
        VannaVolgaSmile(1.1, 0.25, (1.0, 1.1, 1.2), (0.1, 0.1, 0.1), "cubic")
