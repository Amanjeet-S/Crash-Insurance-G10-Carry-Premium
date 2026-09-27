import functools
import math

import numpy as np
import pytest
from scipy.stats import norm

from qef.fx.conventions import G10
from qef.fx.gk import CALL, PUT, strike_from_delta_smile
from qef.fx.identification import SharpMomentSet, boundary
from qef.fx.moments import (
    CONTRACTS,
    contract_intervals,
    contract_weight,
    middle_contracts,
    otm_price,
    tail_bounds,
    weight_sign_changes,
)
from qef.validation.models import (
    Heston,
    Market,
    Merton,
    fourier_price,
    hazard_rates,
    model_smile,
    price_bound,
    quantlib_otm_prices,
    quantlib_prices,
    reference_moments,
    tail_exponent,
)

# An EURUSD-type pair (USD quote currency) and a USD-base pair, one month; FX-typical parameters.
MARKETS = {"EUR": Market(1.10, 0.040, 0.020), "JPY": Market(150.0, 0.005, 0.045)}
MODELS = {
    "merton_down": Merton(sigma=0.08, lam=2.0, mu=-0.04, delta=0.06),
    "merton_sym": Merton(sigma=0.08, lam=1.0, mu=0.0, delta=0.05),
    "heston_neg": Heston(v0=0.01, kappa=2.0, theta=0.01, xi=0.5, rho=-0.3),
    "heston_pos": Heston(v0=0.01, kappa=2.0, theta=0.01, xi=0.4, rho=0.3),
}


@functools.lru_cache(maxsize=None)
def setup(model_name, ccy):
    """Model, market, smile and the smile's 10Δ strikes in the pair's convention."""
    model, market = MODELS[model_name], MARKETS[ccy]
    smile = model_smile(model, market)
    conv = G10[ccy].delta
    K_min = strike_from_delta_smile(-0.10, market.forward, market.tau, PUT, conv, smile, market.df_base)
    K_max = strike_from_delta_smile(0.10, market.forward, market.tau, CALL, conv, smile, market.df_base)
    return model, market, smile, K_min, K_max


@pytest.mark.parametrize("name", ["merton_down", "merton_sym"])
def test_quantlib_bates_with_constant_variance_is_merton(name):
    model, market = MODELS[name], MARKETS["EUR"]
    F = market.forward
    K = F * np.exp(np.linspace(-0.8, 0.8, 81))
    diff = quantlib_otm_prices(model, market, K) - model.series_price(K, F, market.tau)
    assert np.max(np.abs(diff)) < 1e-12 * F


@pytest.mark.parametrize("name", ["heston_neg", "heston_pos"])
def test_quantlib_heston_agrees_with_fourier_prices(name):
    model, market = MODELS[name], MARKETS["EUR"]
    F = market.forward
    K = F * np.exp(np.linspace(-0.8, 0.8, 41))
    diff = quantlib_otm_prices(model, market, K) - fourier_price(model, K, F, market.tau)
    assert np.max(np.abs(diff)) < 1e-12 * F


@pytest.mark.parametrize("name", list(MODELS))
def test_put_call_parity_and_the_forward(name):
    model, market = MODELS[name], MARKETS["JPY"]
    F = market.forward
    K = F * np.exp(np.linspace(-0.3, 0.3, 13))
    c, p = quantlib_prices(model, market, K, CALL), quantlib_prices(model, market, K, PUT)
    np.testing.assert_allclose(c - p, F - K, atol=1e-11 * F)
    assert float(np.real(np.exp(model.cgf(1.0, market.tau)))) == pytest.approx(1.0, abs=1e-12)  # E^q S_T = F


@pytest.mark.parametrize("name", list(MODELS))
def test_distribution_matches_the_strike_derivative_of_prices(name):
    # G(K) = P'(K) and 1 − G(K) = −C'(K) (R2), from the model's own distribution function.
    model, market = MODELS[name], MARKETS["EUR"]
    F, tau = market.forward, market.tau
    price = model.series_price if isinstance(model, Merton) else (lambda K, F, tau: fourier_price(model, K, F, tau))
    for x in (-0.12, -0.05, 0.05, 0.12):
        K, h = F * math.exp(x), 2e-4 * F * math.exp(x)
        d = lambda h: (price(np.array([K + h]), F, tau)[0] - price(np.array([K - h]), F, tau)[0]) / (2 * h)
        slope = (4 * d(h / 2) - d(h)) / 3  # Richardson: error O(h⁴)
        lf, lF, lS = model.log_distribution(np.array([x]), tau)
        target = math.exp(lF[0]) if x < 0 else -math.exp(lS[0])
        assert slope == pytest.approx(target, rel=1e-8)
        assert math.exp(lF[0]) + math.exp(lS[0]) == pytest.approx(1.0, abs=1e-12)


def test_heston_hazard_rates_approach_the_critical_moments():
    # The tails of x are exponential with the critical moments as rates (Keller-Ressel, 2011).
    model, tau = MODELS["heston_neg"], MARKETS["EUR"].tau
    lo, hi = model.strip(tau)
    assert lo < 0 < 1 < hi
    r_far, _ = hazard_rates(model, np.array([-6.0]), tau)
    _, h_far = hazard_rates(model, np.array([6.0]), tau)
    assert r_far[0] == pytest.approx(-lo, rel=0.1) and h_far[0] == pytest.approx(hi, rel=0.1)


def omitted_tail_bound(model, market, smile, usd_base):
    """Bounds on the spanning contributions beyond the smile's node range, from Lee (2004a) through R4."""
    F, tau = market.forward, market.tau
    out = np.zeros(len(CONTRACTS))
    for i, k in enumerate(CONTRACTS):
        w = lambda K, k=k: contract_weight(k, K, F, usd_base)
        breaks = weight_sign_changes(k, F, usd_base)
        for side, K_e in (("lower", smile.K_lo), ("upper", smile.K_hi)):
            B, e = price_bound(model, tau, F, K_e, side)
            out[i] += sum(tail_bounds(w, side, K_e, B, e, breaks)[:2])
    return out


@pytest.mark.parametrize("ccy", ["EUR", "JPY"])
@pytest.mark.parametrize("name", list(MODELS))
def test_untruncated_spanning_reproduces_model_moments(name, ccy):
    model, market, smile, *_ = setup(name, ccy)
    usd_base = G10[ccy].usd_base
    true = reference_moments(model, market.tau, usd_base)
    spanned = middle_contracts(market.forward, smile.K_lo, smile.K_hi, market.tau, smile, usd_base)["values"]
    slack = omitted_tail_bound(model, market, smile, usd_base)
    assert np.all(np.abs(spanned - true) <= slack + 1e-7 * np.abs(true))


@pytest.mark.parametrize("ccy", ["EUR", "JPY"])
@pytest.mark.parametrize("name", list(MODELS))
def test_truncated_intervals_contain_the_truth_at_the_models_own_exponents(name, ccy):
    # R4 and R10 at the largest exponents the model's tails satisfy: containment is a theorem,
    # so a failure would be a bug in the moment code.
    model, market, smile, K_min, K_max = setup(name, ccy)
    F, tau, usd_base = market.forward, market.tau, G10[ccy].usd_base
    x_min, x_max = math.log(K_min / F), math.log(K_max / F)
    gamma = 0.999 * tail_exponent(model, tau, x_min, "lower")["inf"]
    eta = 1.0 + 0.999 * (tail_exponent(model, tau, x_max, "upper")["inf"] - 1.0)
    true = reference_moments(model, tau, usd_base)
    mid = middle_contracts(F, K_min, K_max, tau, smile, usd_base)["values"]
    P_min, C_max = float(otm_price(K_min, F, tau, smile)), float(otm_price(K_max, F, tau, smile))
    lo, hi, *_ = contract_intervals(mid, F, K_min, P_min, K_max, C_max, usd_base, gamma, eta)
    tol = 1e-7 * np.abs(true)
    assert np.all(lo - tol <= true) and np.all(true <= hi + tol)
    bd = boundary(F, tau, smile, K_min, K_max)
    assert gamma <= bd.gamma_bar * (1 + 1e-6) and eta <= bd.eta_bar * (1 + 1e-6)  # R10(a)
    S = SharpMomentSet(F, tau, smile, usd_base, bd, min(gamma, bd.gamma_bar), min(eta, bd.eta_bar), n_c=400, mid=mid)
    v, k3 = S.variance_range(), S.third_central_range()
    var, third = true[1] - true[0] ** 2, true[2] - 3 * true[0] * true[1] + 2 * true[0] ** 3
    assert v[0] - 1e-6 * var <= var <= v[1] + 1e-6 * var
    assert k3[0] - 1e-6 * abs(third) <= third <= k3[1] + 1e-6 * abs(third)
