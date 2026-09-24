import pandas as pd

from qef.data.smile_inputs import option_dates


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
