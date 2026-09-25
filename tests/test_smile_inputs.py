import numpy as np
import pandas as pd
import pytest

from qef.data.panel import CURRENCIES, instrument_catalogue, parse_ric
from qef.data.smile_inputs import build_month_end_inputs, input_rics, locate_raw, option_dates, rate_ric


def _d(s):
    return pd.Timestamp(s)


def test_spot_rolls_over_us_holiday():
    # Friday 30 June 2023: T+2 skips 4 July.
    spot, delivery, expiry = option_dates(_d("2023-06-30"), "EUR")
    assert spot == _d("2023-07-05") and delivery == _d("2023-08-07") and expiry == _d("2023-08-03")


def test_weekend_delivery_rolls_forward():
    spot, delivery, expiry = option_dates(_d("2024-01-31"), "EUR")
    assert spot == _d("2024-02-02") and delivery == _d("2024-03-04") and expiry == _d("2024-02-29")


def test_modified_following_rolls_back_at_month_end():
    # Spot 31 July 2024; +1M is Saturday 31 August; the following business day
    # (3 September, after Labor Day) is in the next month, so delivery rolls back to 30 August.
    spot, delivery, expiry = option_dates(_d("2024-07-29"), "EUR")
    assert spot == _d("2024-07-31") and delivery == _d("2024-08-30") and expiry == _d("2024-08-28")


def test_cad_spot_lag_is_one_day():
    spot, _, _ = option_dates(_d("2019-06-28"), "CAD")
    assert spot == _d("2019-07-01")


# ---------------------------------------------------------------------------
# Three-month tenor


def test_three_month_holiday_delivery_and_expiry():
    # Trade Friday 31 March 2023; spot is Tuesday 4 April (T+2 over the weekend).
    # 4 April + 3 months is Tuesday 4 July, Independence Day, so delivery is
    # Wednesday 5 July. Two business days back skips 4 July and the weekend:
    # Monday 3 July, then Friday 30 June.
    spot, delivery, expiry = option_dates(_d("2023-03-31"), "EUR", months=3)
    assert spot == _d("2023-04-04") and delivery == _d("2023-07-05") and expiry == _d("2023-06-30")


def test_three_month_modified_following_rolls_back_at_month_end():
    # Trade Wednesday 29 May 2024; spot Friday 31 May. +3M is Saturday 31 August;
    # 1 September is a Sunday and 2 September Labor Day, so the following business
    # day (3 September) is in the next month and delivery rolls back to Friday 30 August.
    spot, delivery, expiry = option_dates(_d("2024-05-29"), "EUR", months=3)
    assert spot == _d("2024-05-31") and delivery == _d("2024-08-30") and expiry == _d("2024-08-28")


def test_three_month_cad_one_day_spot_lag():
    # Same spot and delivery as above from a trade one day later (T+1); expiry is
    # one business day before delivery.
    spot, delivery, expiry = option_dates(_d("2024-05-30"), "CAD", months=3)
    assert spot == _d("2024-05-31") and delivery == _d("2024-08-30") and expiry == _d("2024-08-29")
    # Quarter-end trade, Friday 28 June 2019: CAD spot Monday 1 July (Canada Day is
    # not a New York holiday), delivery Tuesday 1 October, expiry Monday 30 September.
    # EUR spot is Tuesday 2 July and delivery Wednesday 2 October; both expire on 30 September.
    assert option_dates(_d("2019-06-28"), "CAD", months=3) == (_d("2019-07-01"), _d("2019-10-01"), _d("2019-09-30"))
    assert option_dates(_d("2019-06-28"), "EUR", months=3) == (_d("2019-07-02"), _d("2019-10-02"), _d("2019-09-30"))


def test_three_month_end_of_month_clipping():
    # Spot Thursday 30 November 2023; 30 February does not exist, so +3M is the
    # last day of February 2024, Thursday 29 February, a business day.
    spot, delivery, expiry = option_dates(_d("2023-11-28"), "EUR", months=3)
    assert spot == _d("2023-11-30") and delivery == _d("2024-02-29") and expiry == _d("2024-02-27")


def test_delivery_must_fall_inside_holiday_table():
    # The holiday table ends on 31 December 2040.
    option_dates(_d("2040-10-15"), "EUR", months=1)
    with pytest.raises(ValueError):
        option_dates(_d("2040-10-15"), "EUR", months=3)
    # One month. Trade Wednesday 28 November 2040: spot Friday 30 November; +1M is
    # Sunday 30 December, so delivery is Monday 31 December, the last day of the table.
    assert option_dates(_d("2040-11-28"), "EUR")[1] == _d("2040-12-31")
    # Trade Thursday 29 November: spot Monday 3 December, delivery 3 January 2041.
    # Two business days back would give 1 January 2041, a holiday the table does not hold.
    with pytest.raises(ValueError):
        option_dates(_d("2040-11-29"), "EUR")
    # Three months. Trade Wednesday 26 September 2040: spot Friday 28 September,
    # delivery Friday 28 December. A day later, spot is Monday 1 October and +3M is 1 January 2041.
    assert option_dates(_d("2040-09-26"), "EUR", months=3)[1] == _d("2040-12-28")
    with pytest.raises(ValueError):
        option_dates(_d("2040-09-27"), "EUR", months=3)
    with pytest.raises(ValueError):
        option_dates(_d("1989-12-29"), "EUR")


# ---------------------------------------------------------------------------
# Rate RICs

# The quote currency is USD for EURUSD, GBPUSD, AUDUSD and NZDUSD, and the
# non-USD currency for USDJPY, USDCHF, USDCAD, USDNOK and USDSEK. USD uses OIS;
# the others use deposits.
_QUOTE_CCY = {"EUR": "USD", "GBP": "USD", "AUD": "USD", "NZD": "USD",
              "JPY": "JPY", "CHF": "CHF", "CAD": "CAD", "NOK": "NOK", "SEK": "SEK"}


@pytest.mark.parametrize("tenor", ["1M", "3M"])
def test_rate_ric_has_the_option_tenor(tenor):
    for ccy, q in _QUOTE_CCY.items():
        expected = f"USD{tenor}OIS=" if q == "USD" else f"{q}{tenor}D="
        assert input_rics(ccy, tenor=tenor)["rate"] == ("rates", expected)
    assert rate_ric("USD", "3M") == ("USD3MOIS=", "rates")
    assert rate_ric("JPY", "3M") == ("JPY3MD=", "rates")
    assert rate_ric("USD", "1M") == ("USD1MOIS=", "rates")


def test_rate_ric_rejects_unknown_tenor():
    with pytest.raises(ValueError):
        rate_ric("USD", "6M")


@pytest.mark.parametrize("tenor", ["1M", "3M"])
@pytest.mark.parametrize("contributor", ["", "FN"])
def test_input_rics_are_in_the_catalogue(tenor, contributor):
    catalogue = {i.ric: i for i in instrument_catalogue()}
    for ccy in CURRENCIES:
        for block, ric in input_rics(ccy, contributor, tenor).values():
            assert ric in catalogue
            assert parse_ric(ric).block == block
            assert parse_ric(ric).tenor in (None, tenor)
    rics = input_rics("EUR", contributor, tenor)
    assert rics["fwd"] == ("forward", f"EUR{tenor}=")
    assert rics["rr10"] == ("vol_rr10", f"EUR{tenor}R10=" + contributor)


# ---------------------------------------------------------------------------
# Raw roots searched in order


def _write(root, block, ric, rows, columns=("BID", "ASK")):
    path = root / block / f"{ric}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["Date", *columns]).to_csv(path, index=False)
    return path


def test_locate_raw_uses_the_first_root_holding_the_file(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    only_b = _write(b, "forward", "EUR3M=", [("2023-03-31", 1.0, 2.0)])
    in_a = _write(a, "rates", "USD3MOIS=", [("2023-03-31", 1.0, 2.0)])
    _write(b, "rates", "USD3MOIS=", [("2023-03-31", 3.0, 4.0)])
    assert locate_raw([a, b], "forward", "EUR3M=") == only_b
    assert locate_raw([a, b], "rates", "USD3MOIS=") == in_a
    with pytest.raises(FileNotFoundError):
        locate_raw([a, b], "rates", "GBP3MD=")


@pytest.mark.parametrize("days_mat", [93.0, None])
def test_three_month_inputs_from_two_roots(tmp_path, days_mat):
    """Synthetic EURUSD month-end, 31 March 2023, with the three-month series in a second root."""
    main, extra = tmp_path / "main", tmp_path / "extra"
    t = "2023-03-31"
    _write(main, "spot", "EUR=", [(t, 1.0860, 1.0862)])
    # One-month series in the main root must not be used at three months.
    _write(main, "forward", "EUR1M=", [(t, 10.0, 12.0)])
    _write(main, "rates", "USD1MOIS=", [(t, 9.0, 9.0)])
    for block, code, (bid, ask) in (("vol_atm", "O", (8.0, 8.4)), ("vol_rr25", "RR", (-0.6, -0.4)),
                                    ("vol_bf25", "BF", (0.2, 0.3))):
        _write(main, block, f"EUR3M{code}=", [(t, bid, ask)])
    _write(main, "vol_rr10", "EUR3MR10=", [("2023-03-24", -1.2, -0.8)])  # five business days earlier
    _write(main, "vol_bf10", "EUR3MB10=", [("2023-03-01", 0.7, 0.9)])  # too old
    if days_mat is None:
        _write(extra, "forward", "EUR3M=", [(t, 50.0, 52.0)])
    else:
        _write(extra, "forward", "EUR3M=", [(t, 50.0, 52.0, days_mat)], columns=("BID", "ASK", "DAYS_MAT"))
    _write(extra, "rates", "USD3MOIS=", [(t, 4.80, 4.82)])

    out = build_month_end_inputs(main, [_d(t)], ["EUR"], tenor="3M", extra_roots=[extra], quote_status=True)
    r = out.iloc[0]
    S = 1.0861
    assert r.S == pytest.approx(S, rel=1e-15)
    assert r.F == pytest.approx(S + 51.0e-4, rel=1e-15)
    assert r.r_quote == pytest.approx(0.0481, rel=1e-14)
    # Expiry 30 June 2023 (test above): 91 days after the trade date.
    assert r.tau == pytest.approx(91 / 365, rel=1e-15)
    # Spot 4 April to delivery 5 July: 92 days.
    days = 92.0 if days_mat is None else days_mat
    assert r.days == days and r.days_source == ("calendar" if days_mat is None else "DAYS_MAT")
    assert r.df_quote == pytest.approx(1 / (1 + 0.0481 * days / 360), rel=1e-14)
    assert r.df_base == pytest.approx(r.F * r.df_quote / S, rel=1e-14)
    assert (r.atm, r.rr25, r.bf25, r.rr10) == pytest.approx((0.082, -0.005, 0.0025, -0.01), rel=1e-12)
    assert np.isnan(r.bf10)
    assert (r.status_rr10, r.status_bf10, r.status_rate) == ("substituted", "missing", "observed")
    assert r.n_substituted == 1
    assert r.complete_calib and not r.has_10d and not r.complete

    plain = build_month_end_inputs(main, [_d(t)], ["EUR"], tenor="3M", extra_roots=[extra])
    assert not any(c.startswith("status_") for c in plain.columns)
    with pytest.raises(FileNotFoundError):
        build_month_end_inputs(main, [_d(t)], ["EUR"], tenor="3M")
