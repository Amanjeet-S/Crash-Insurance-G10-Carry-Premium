import numpy as np
import pandas as pd
import pytest

from qef.data.panel import (
    forward_outright,
    instrument_catalogue,
    invert_bid_ask,
    ny_month_ends,
    outright_bid_ask,
    parse_ric,
    sample_month_ends,
    stale_mask,
    two_sided,
    unchanged_runs,
    usd_per_unit,
)


def test_parse_ric_classifies_each_block():
    i = parse_ric("NOK1MRR=FN")
    assert (i.block, i.currency, i.tenor, i.quote, i.contributor) == ("vol_rr25", "NOK", "1M", "rr25", "fenics")
    assert parse_ric("EUR3MB10=").quote == "bf10"
    assert parse_ric("EUR1MO=TIFO").contributor == "tifo"
    assert (parse_ric("EUR1M=").block, parse_ric("EUR1M=").quote) == ("forward", "points")
    assert parse_ric("EUR1MD=").quote == "deposit"
    assert parse_ric("EUR1MOIS=").quote == "ois"
    assert (parse_ric("EUR3M=").block, parse_ric("EUR3M=").tenor) == ("forward", "3M")
    assert (parse_ric("USD3MOIS=").quote, parse_ric("JPY3MD=").tenor) == ("ois", "3M")
    assert parse_ric("USDSROIS3M=").quote == "sofr_ois"
    assert parse_ric("JPY=").quote == "spot"
    assert parse_ric("VXc2").block == "vix"
    with pytest.raises(ValueError):
        parse_ric("EUR6MO=")


def test_catalogue_is_unique_and_complete():
    rics = [i.ric for i in instrument_catalogue()]
    assert len(rics) == len(set(rics))
    # 9 spot + 9 forward + 9 x 2 tenors x 5 quotes x 2 contributors + 3 TIFO + 3 USD + 18 rates + 2 VIX
    # + three-month: 9 forwards + 3 USD rates + 9 deposits
    assert len(rics) == 9 + 9 + 180 + 3 + 3 + 18 + 2 + 9 + 3 + 9


def test_ny_month_ends_holidays_and_weekends():
    me = ny_month_ends("2021-05-01", "2022-01-31")
    assert pd.Timestamp("2021-05-28") in me  # 31 May 2021 is Memorial Day
    assert pd.Timestamp("2021-07-30") in me  # 31 July 2021 is a Saturday
    assert pd.Timestamp("2021-12-31") in me  # Saturday New Year's Day is not moved to Friday
    assert len(me) == 9
    assert ny_month_ends("2026-10-01", "2026-10-29").empty  # month not complete


def _quotes(dates, bid, ask):
    return pd.DataFrame({"BID": bid, "ASK": ask}, index=pd.DatetimeIndex(dates))


def test_two_sided_requires_both_sides():
    q = _quotes(["2024-01-02", "2024-01-03", "2024-01-04"], [1.0, np.nan, 1.1], [1.2, 1.3, np.nan])
    assert two_sided(q).tolist() == [True, False, False]


def test_month_end_substitution_rule():
    # Month-ends 31 January, 29 February and 29 March 2024 (Good Friday is not a
    # Federal Reserve holiday).
    q = _quotes(
        ["2024-01-31", "2024-02-22", "2024-02-29", "2024-03-21", "2024-03-29"],
        [1.0, 2.0, 2.5, 3.0, 4.0],
        [1.2, 2.2, np.nan, 3.2, np.nan],
    )
    me = ny_month_ends("2024-01-01", "2024-03-31")
    out = sample_month_ends(q, me, window=5)
    assert out["status"].tolist() == ["observed", "substituted", "missing"]
    assert out.loc["2024-02-29", "source_date"] == pd.Timestamp("2024-02-22")
    assert out.loc["2024-02-29", "lag"] == 5
    assert out.loc["2024-01-31", "mid"] == pytest.approx(1.1)
    # 21 March lies six business days before 29 March, outside the window
    assert np.isnan(out.loc["2024-03-29", "bid"])


def test_stale_runs_flag_five_or_more():
    dates = pd.bdate_range("2024-01-01", periods=12)
    bid = [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 4, 4]
    q = _quotes(dates, bid, [b + 0.5 for b in bid])
    runs = unchanged_runs(q)
    assert runs["length"].tolist() == [5, 4, 1, 2]
    mask = stale_mask(q, min_run=5)
    assert mask.sum() == 5 and mask.iloc[:5].all()


def test_forward_outright_and_usd_inversion():
    # USDJPY: points are in units of 0.01 yen
    assert forward_outright(150.0, -40.0, 1e-2) == pytest.approx(149.6)
    fb, fa = outright_bid_ask(149.98, 150.02, -40.2, -39.8, 1e-2)
    assert fa > fb
    xb, xa = usd_per_unit(fb, fa, usd_base=True)
    assert xb == pytest.approx(1 / fa) and xa == pytest.approx(1 / fb) and xb < xa
    assert usd_per_unit(1.1, 1.2, usd_base=False) == (1.1, 1.2)
    b, a = invert_bid_ask(np.array([2.0]), np.array([4.0]))
    assert (b[0], a[0]) == (0.25, 0.5)
