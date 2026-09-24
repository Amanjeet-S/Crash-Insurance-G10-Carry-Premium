import numpy as np
import pytest

from qef.fx.crash import forward_discount, leg_skew_cost, protective_option, rank_legs
from qef.fx.gk import CALL, PUT
from qef.fx.sabr import sabr_vol


def test_protective_option_mapping():
    # Long EUR loses when EURUSD falls (put); long JPY loses when USDJPY rises (USD call).
    assert protective_option("EUR", True) == PUT and protective_option("EUR", False) == CALL
    assert protective_option("JPY", True) == CALL and protective_option("JPY", False) == PUT


def test_forward_discount_sign_is_interest_differential():
    # A currency at a forward discount to USD (higher interest) has positive fd in both quotations.
    assert forward_discount(0.65, 0.649, usd_base=False) > 0  # AUDUSD forward below spot
    assert forward_discount(150.0, 150.5, usd_base=True) > 0  # USDJPY forward above spot: JPY at a forward discount


def test_flat_smile_has_zero_skew_cost():
    flat = lambda K: np.full_like(np.asarray(K, float), 0.09)
    for ccy, long_leg in (("EUR", True), ("JPY", True), ("AUD", False), ("CHF", False)):
        _, v_s, v_f = leg_skew_cost(flat, 1.1 if ccy in ("EUR", "AUD") else 150.0, 1 / 12, 0.998, ccy, 0.09, long_leg)
        assert v_s == pytest.approx(v_f, rel=1e-12)


def test_downside_skew_makes_crash_protection_expensive():
    F, tau = 0.65, 1 / 12
    skewed = lambda K: sabr_vol(K, F, tau, 0.10, -0.5, 2.0)  # puts on AUDUSD rich
    _, v_s, v_f = leg_skew_cost(skewed, F, tau, 0.998, "AUD", float(skewed(F)), True)
    assert v_s > v_f > 0


def test_rank_legs():
    longs, shorts = rank_legs({"A": 0.3, "B": -0.1, "C": 0.2, "D": 0.0}, 1)
    assert longs == ["A"] and shorts == ["B"]


def test_hedged_leg_loss_is_floored_at_the_strike():
    # EURUSD long, put at K: hedged value at delivery is at least (K - F)/F per USD notional.
    from qef.fx.crash import forward_return, option_payoff
    F, K = 1.10, 1.06
    for S_end in (0.8, 1.0, 1.06, 1.2):
        hedged = forward_return("EUR", True, S_end, F) + option_payoff("EUR", True, K, F, S_end)
        assert hedged >= (K - F) / F - 1e-12
    # USDJPY long JPY (USD call at K): hedged value is at least (F - K)/S_end.
    F, K = 150.0, 155.0
    for S_end in (140.0, 155.0, 170.0, 200.0):
        hedged = forward_return("JPY", True, S_end, F) + option_payoff("JPY", True, K, F, S_end)
        assert hedged >= (F - K) / S_end - 1e-12
