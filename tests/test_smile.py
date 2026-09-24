import numpy as np
import pytest

from qef.fx.arbitrage import check_smile
from qef.fx.gk import DeltaConvention
from qef.fx.sabr import sabr_vol
from qef.fx.smile import SabrSmile, SmileQuotes, calibrate_sabr, quotes_from_smile

PIPS = DeltaConvention(spot=True, premium_adjusted=False)
PA = DeltaConvention(spot=True, premium_adjusted=True)


def test_sabr_rho_zero_is_symmetric_in_log_moneyness():
    F, tau = 1.1, 0.25
    x = np.array([0.02, 0.05, 0.1])
    up = sabr_vol(F * np.exp(x), F, tau, 0.1, 0.0, 1.2)
    down = sabr_vol(F * np.exp(-x), F, tau, 0.1, 0.0, 1.2)
    np.testing.assert_allclose(up, down, rtol=1e-12)


def test_sabr_flat_when_vol_of_vol_zero_and_beta_one():
    K = np.linspace(0.8, 1.4, 7)
    np.testing.assert_allclose(sabr_vol(K, 1.1, 0.5, 0.12, 0.3, 0.0), 0.12, rtol=1e-14)


def test_sabr_continuous_at_the_money():
    F, tau = 1.1, 0.25
    eps = 1e-9
    a = sabr_vol(F * (1 + eps), F, tau, 0.1, -0.3, 1.5)
    b = sabr_vol(F, F, tau, 0.1, -0.3, 1.5)
    assert a == pytest.approx(b, rel=1e-7)


@pytest.mark.parametrize("conv", [PIPS, PA])
@pytest.mark.parametrize("reading", ["market", "smile"])
@pytest.mark.parametrize("params", [(0.08, -0.25, 1.2), (0.12, 0.35, 2.5), (0.10, 0.0, 0.8)])
def test_synthetic_recovery(conv, reading, params):
    F, tau, db = 1.1, 1 / 12, 0.998
    truth = SabrSmile(F, tau, *params)
    atm, rr, bf = quotes_from_smile(truth, tau, conv, 0.25, db, reading)
    res = calibrate_sabr(SmileQuotes(F, tau, atm, rr, bf, 0.25, db, conv), bf_reading=reading)
    assert res.success, res.message
    np.testing.assert_allclose([res.smile.alpha, res.smile.rho, res.smile.nu], params, rtol=1e-6, atol=1e-8)


def test_market_and_smile_readings_differ_for_skewed_smile():
    F, tau = 1.1, 1 / 12
    truth = SabrSmile(F, tau, 0.10, -0.5, 2.0)
    _, _, bf_ms = quotes_from_smile(truth, tau, PA, 0.25, 1.0, "market")
    _, _, bf_ss = quotes_from_smile(truth, tau, PA, 0.25, 1.0, "smile")
    assert abs(bf_ms - bf_ss) > 1e-5


def test_flat_smile_is_arbitrage_free():
    rep = check_smile(lambda K: np.full_like(np.asarray(K, float), 0.1), 1.1, 0.25)
    assert rep.convex_ok and rep.slope_bounds_ok and rep.density_ok


def test_arbitrageable_smile_detected():
    # A sharp volatility spike produces negative density.
    spike = lambda K: 0.1 + 0.5 * np.exp(-((np.log(np.asarray(K) / 1.1)) / 0.01) ** 2)
    rep = check_smile(spike, 1.1, 0.25)
    assert not rep.density_ok


def test_ten_delta_market_quotes_for_low_volatility_smile():
    # Regression: the butterfly bracket must keep sigma_ATM + BF positive when ATM volatility is low.
    # Synthetic values; the code before the fix does not terminate on them.
    truth = SabrSmile(1.30, 1 / 12, 0.04, 0.10, 2.0)
    atm, rr, bf = quotes_from_smile(truth, 1 / 12, PA, 0.10, 0.999, "market")
    assert 0 < atm < 0.1 and bf > 0
