import numpy as np
import pytest

from qef.fx.gk import (
    CALL,
    PUT,
    DeltaConvention,
    atm_dns_strike,
    delta,
    forward_premium,
    pa_call_delta_maximiser,
    premium,
    strike_from_delta,
)

CONVENTIONS = [
    DeltaConvention(spot=True, premium_adjusted=False),
    DeltaConvention(spot=True, premium_adjusted=True),
    DeltaConvention(spot=False, premium_adjusted=False),
    DeltaConvention(spot=False, premium_adjusted=True),
]
S, DQ, DB, TAU = 1.10, 0.996, 0.998, 1 / 12
F = S * DB / DQ


def test_put_call_parity():
    K = np.linspace(0.9, 1.3, 9)
    c = premium(F, K, 0.08, TAU, DQ, CALL)
    p = premium(F, K, 0.08, TAU, DQ, PUT)
    np.testing.assert_allclose(c - p, DQ * (F - K), atol=1e-15)


def test_spot_delta_matches_finite_difference():
    # V(S) with forward F = S D_b / D_q; pips spot delta is dV/dS.
    K, sig, h = 1.12, 0.09, 1e-6
    v = lambda s: premium(s * DB / DQ, K, sig, TAU, DQ, CALL)
    fd = (v(S + h) - v(S - h)) / (2 * h)
    assert delta(F, K, sig, TAU, CALL, CONVENTIONS[0], DB) == pytest.approx(fd, rel=1e-8)


def test_premium_adjusted_delta_is_delta_less_premium_in_base_units():
    K, sig = 1.08, 0.11
    for phi in (CALL, PUT):
        pips = delta(F, K, sig, TAU, phi, CONVENTIONS[0], DB)
        pa = delta(F, K, sig, TAU, phi, CONVENTIONS[1], DB)
        assert pa == pytest.approx(pips - premium(F, K, sig, TAU, DQ, phi) / S, abs=1e-14)


@pytest.mark.parametrize("conv", CONVENTIONS)
def test_atm_dns_strike_is_delta_neutral(conv):
    for sig in (0.05, 0.12, 0.3):
        K = atm_dns_strike(F, sig, TAU, conv)
        total = delta(F, K, sig, TAU, CALL, conv, DB) + delta(F, K, sig, TAU, PUT, conv, DB)
        assert total == pytest.approx(0.0, abs=1e-14)


@pytest.mark.parametrize("conv", CONVENTIONS)
@pytest.mark.parametrize("target", [0.10, 0.25, 0.40])
@pytest.mark.parametrize("sig,tau", [(0.07, 1 / 52), (0.1, 1 / 12), (0.25, 1.0)])
def test_strike_from_delta_round_trip(conv, target, sig, tau):
    for phi in (CALL, PUT):
        K = strike_from_delta(phi * target, F, sig, tau, phi, conv, DB)
        assert delta(F, K, sig, tau, phi, conv, DB) == pytest.approx(phi * target, abs=1e-12)


def test_premium_adjusted_call_delta_unimodal_and_right_branch():
    conv = CONVENTIONS[3]  # forward, premium adjusted
    sig, tau = 0.3, 5.0
    k_star = pa_call_delta_maximiser(F, sig, tau)
    ks = F * np.exp(np.linspace(-3, 3, 4001))
    d = delta(F, ks, sig, tau, CALL, conv)
    assert ks[np.argmax(d)] == pytest.approx(k_star, rel=2e-3)
    inc = np.diff(d[ks < k_star]); dec = np.diff(d[ks > k_star])
    assert np.all(inc > 0) and np.all(dec < 0)
    target = 0.5 * d.max()
    K = strike_from_delta(target, F, sig, tau, CALL, conv)
    assert K > k_star
    with pytest.raises(ValueError):
        strike_from_delta(1.01 * d.max(), F, sig, tau, CALL, conv)


def test_forward_premium_limits():
    assert forward_premium(F, 1e-9, 0.1, TAU, CALL) == pytest.approx(F, rel=1e-9)
    assert forward_premium(F, 1e6, 0.1, TAU, CALL) == pytest.approx(0.0, abs=1e-15)


def test_diffusive_null_bound_r6():
    # Theory note R6: theta_0 = Phi(-d_+(xi)) and |theta_0 - Phi(-d1)| <= phi(0)|lambda| sqrt(tau)/sigma.
    from scipy.stats import norm

    F0, K, sig, tau = 1.1, 1.06, 0.09, 1 / 12
    p = lambda G: forward_premium(G, K, sig, tau, PUT)
    d1 = (np.log(F0 / K) + 0.5 * sig**2 * tau) / (sig * np.sqrt(tau))
    for lam in (-0.2, 0.05, 0.3):
        mu_u = F0 * np.expm1(lam * tau)
        theta0 = -(p(F0 * np.exp(lam * tau)) - p(F0)) / mu_u
        assert abs(theta0 - norm.cdf(-d1)) <= norm.pdf(0) * abs(lam) * np.sqrt(tau) / sig


@pytest.mark.parametrize("sigma", [0.0, -0.05, np.nan])
def test_non_positive_volatility_is_rejected(sigma):
    # Regression: a negative volatility once sent the premium-adjusted maximiser into an endless loop.
    conv = DeltaConvention(spot=True, premium_adjusted=True)
    with pytest.raises(ValueError):
        strike_from_delta(0.25, F, sigma, TAU, CALL, conv, DB)
    with pytest.raises(ValueError):
        pa_call_delta_maximiser(F, sigma, TAU)
