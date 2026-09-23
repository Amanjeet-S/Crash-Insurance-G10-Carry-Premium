"""Instrument names and sampling rules for the LSEG FX panel.

The raw layer holds one table per RIC, indexed by observation date, with BID
and ASK among its columns. The helpers here name the instruments, select New
York month-ends, apply the two-sided and five-business-day substitution rules
of the research design (section 4), measure runs of unchanged quotes, and
build forward outrights in USD per unit of the non-USD currency. They hold no
LSEG values.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.tseries.holiday import (
    AbstractHolidayCalendar,
    Holiday,
    USColumbusDay,
    USLaborDay,
    USMartinLutherKingJr,
    USMemorialDay,
    USPresidentsDay,
    USThanksgivingDay,
    sunday_to_monday,
)

CURRENCIES = ("EUR", "GBP", "AUD", "NZD", "JPY", "CHF", "CAD", "NOK", "SEK")
TENORS = ("1M", "3M")
# RIC code -> quote type. R10 and B10 are the 10-delta risk reversal and butterfly.
VOL_QUOTES = {"O": "atm", "RR": "rr25", "BF": "bf25", "R10": "rr10", "B10": "bf10"}
CONTRIBUTORS = {"": "composite", "FN": "fenics", "TIFO": "tifo"}


@dataclass(frozen=True)
class Instrument:
    ric: str
    block: str  # storage block under raw/
    currency: str | None = None
    tenor: str | None = None
    quote: str | None = None  # spot, points, atm, rr25, bf25, rr10, bf10, deposit, ois, sofr_ois, future
    contributor: str = "composite"


_VOL = re.compile(r"^([A-Z]{3})(1M|3M)(O|RR|BF|R10|B10)=(FN|TIFO)?$")
_SPOT = re.compile(r"^([A-Z]{3})=$")
_FWD = re.compile(r"^([A-Z]{3})1M=$")
_DEP = re.compile(r"^([A-Z]{3})1MD=$")
_OIS = re.compile(r"^([A-Z]{3})1MOIS=$")
_VIX = re.compile(r"^VXc(\d)$")


def parse_ric(ric: str) -> Instrument:
    """Classify a RIC of the project panel by block, currency, tenor, quote and contributor."""
    if m := _VOL.match(ric):
        ccy, tenor, code, contrib = m.groups()
        quote = VOL_QUOTES[code]
        return Instrument(ric, f"vol_{quote}", ccy, tenor, quote, CONTRIBUTORS[contrib or ""])
    if m := _SPOT.match(ric):
        return Instrument(ric, "spot", m[1], None, "spot")
    if m := _FWD.match(ric):
        return Instrument(ric, "forward", m[1], "1M", "points")
    if m := _DEP.match(ric):
        return Instrument(ric, "rates", m[1], "1M", "deposit")
    if m := _OIS.match(ric):
        return Instrument(ric, "rates", m[1], "1M", "ois")
    if ric == "USDSROIS1M=":
        return Instrument(ric, "rates", "USD", "1M", "sofr_ois")
    if m := _VIX.match(ric):
        return Instrument(ric, "vix", None, f"c{m[1]}", "future", "exchange")
    raise ValueError(f"unrecognised RIC {ric!r}")


def instrument_catalogue() -> list[Instrument]:
    """All RICs requested for project 01 (data plan, LSEG Workspace table).

    USD1MD= is included so that the forward-point sign check compares deposit
    rates with deposit rates.
    """
    rics = [f"{c}=" for c in CURRENCIES]
    rics += [f"{c}1M=" for c in CURRENCIES]
    for suffix in ("", "FN"):
        rics += [f"{c}{t}{q}={suffix}" for c in CURRENCIES for t in TENORS for q in VOL_QUOTES]
    rics += ["EUR1MO=TIFO", "EUR1MRR=TIFO", "EUR1MBF=TIFO"]
    rics += ["USD1MOIS=", "USDSROIS1M=", "USD1MD="]
    rics += [f"{c}1M{k}=" for c in CURRENCIES for k in ("D", "OIS")]
    rics += ["VXc1", "VXc2"]
    return [parse_ric(r) for r in rics]


def raw_path(root: Path, inst: Instrument) -> Path:
    return Path(root) / "raw" / inst.block / f"{inst.ric}.csv"


def read_raw(path: Path) -> pd.DataFrame:
    """Read a raw table written by the acquisition script (date index, full precision)."""
    return pd.read_csv(path, index_col=0, parse_dates=True)


# ---------------------------------------------------------------------------
# New York business calendar


class NewYorkCalendar(AbstractHolidayCalendar):
    """Federal Reserve holidays.

    A holiday on Sunday is observed on Monday. A holiday on Saturday is not
    moved, so the preceding Friday (for example 31 December 2021) remains a
    business day, as in the FX market.
    """

    rules = [
        Holiday("New Year's Day", month=1, day=1, observance=sunday_to_monday),
        USMartinLutherKingJr,
        USPresidentsDay,
        USMemorialDay,
        Holiday("Juneteenth", month=6, day=19, start_date="2022-01-01", observance=sunday_to_monday),
        Holiday("Independence Day", month=7, day=4, observance=sunday_to_monday),
        USLaborDay,
        USColumbusDay,
        Holiday("Veterans Day", month=11, day=11, observance=sunday_to_monday),
        USThanksgivingDay,
        Holiday("Christmas Day", month=12, day=25, observance=sunday_to_monday),
    ]


def ny_business_days(start, end) -> pd.DatetimeIndex:
    holidays = NewYorkCalendar().holidays(start, end)
    return pd.bdate_range(start, end, freq="C", holidays=holidays)


def ny_month_ends(start, end) -> pd.DatetimeIndex:
    """Last New York business day of each calendar month, restricted to [start, end]."""
    start, end = pd.Timestamp(start), pd.Timestamp(end)
    days = ny_business_days(start.to_period("M").start_time, end.to_period("M").end_time.normalize())
    last = days.to_series().groupby(days.to_period("M")).max()
    out = pd.DatetimeIndex(last.values)
    return out[(out >= start) & (out <= end)]


def on_business_days(frame: pd.DataFrame) -> pd.DataFrame:
    """Rows dated on New York business days, with duplicate dates reduced to the last row."""
    frame = frame[~frame.index.duplicated(keep="last")].sort_index()
    if frame.empty:
        return frame
    days = ny_business_days(frame.index.min(), frame.index.max())
    return frame.loc[frame.index.isin(days)]


# ---------------------------------------------------------------------------
# Two-sided quotes and the month-end rule


def two_sided(frame: pd.DataFrame, bid: str = "BID", ask: str = "ASK") -> pd.Series:
    """True where both bid and ask are present and finite.

    Crossed and zero-width quotes count as two-sided here; the audit flags
    them separately.
    """
    if bid not in frame or ask not in frame:
        return pd.Series(False, index=frame.index)
    b = pd.to_numeric(frame[bid], errors="coerce").astype(float)
    a = pd.to_numeric(frame[ask], errors="coerce").astype(float)
    return pd.Series(np.isfinite(b.to_numpy()) & np.isfinite(a.to_numpy()), index=frame.index)


def mid(frame: pd.DataFrame, bid: str = "BID", ask: str = "ASK") -> pd.Series:
    """(bid + ask)/2 where the quote is two-sided, missing otherwise."""
    ok = two_sided(frame, bid, ask)
    if not ok.any():
        return pd.Series(np.nan, index=frame.index)
    return ((frame[bid].astype(float) + frame[ask].astype(float)) / 2).where(ok)


def sample_month_ends(
    frame: pd.DataFrame, month_ends, window: int = 5, bid: str = "BID", ask: str = "ASK"
) -> pd.DataFrame:
    """Apply the month-end rule of the research design (section 4).

    For each month-end m the quote dated m is used if it is two-sided
    ("observed"). Otherwise the last two-sided quote dated within the
    ``window`` New York business days before m is used ("substituted").
    Otherwise the value is "missing". Duplicate dates keep the last row.
    Returns one row per month-end: source date, lag in business days,
    status, and bid, ask and mid of the selected quote.
    """
    month_ends = pd.DatetimeIndex(month_ends)
    cols = ["source_date", "lag", "status", "bid", "ask"]
    ok = frame.loc[two_sided(frame, bid, ask), [bid, ask]].astype(float)
    ok = ok[~ok.index.duplicated(keep="last")].sort_index()
    if ok.empty or month_ends.empty:
        out = pd.DataFrame({"source_date": pd.NaT, "lag": np.nan, "status": "missing",
                            "bid": np.nan, "ask": np.nan}, index=month_ends)[cols]
        out["mid"] = np.nan
        return out
    dates = ok.index
    bdays = ny_business_days(min(dates.min(), month_ends.min()) - pd.Timedelta(days=14), month_ends.max())
    pos = dates.searchsorted(month_ends, side="right") - 1
    rows = []
    for m, i in zip(month_ends, pos):
        if i < 0:
            rows.append((pd.NaT, np.nan, "missing", np.nan, np.nan))
            continue
        d = dates[i]
        # business days in (d, m]: 0 when the quote is dated m itself
        lag = int(bdays.searchsorted(m, side="right") - bdays.searchsorted(d, side="right"))
        if d == m:
            status = "observed"
        elif lag <= window:
            status = "substituted"
        else:
            rows.append((pd.NaT, np.nan, "missing", np.nan, np.nan))
            continue
        b, a = ok.iloc[i]
        rows.append((d, lag, status, b, a))
    out = pd.DataFrame(rows, index=month_ends, columns=cols)
    out["mid"] = (out["bid"] + out["ask"]) / 2
    return out


# ---------------------------------------------------------------------------
# Staleness


def _run_ids(frame: pd.DataFrame, columns) -> pd.Series:
    x = frame.loc[:, list(columns)].dropna()
    x = x[~x.index.duplicated(keep="last")].sort_index()
    changed = x.ne(x.shift()).any(axis=1)
    return changed.cumsum()


def unchanged_runs(frame: pd.DataFrame, columns=("BID", "ASK")) -> pd.DataFrame:
    """Runs of consecutive identical quotes.

    Rows with a missing column are dropped first, so a run continues across a
    gap in the data. Length is counted in observations; pass a frame
    restricted to business days (``on_business_days``) to count business days.
    """
    ids = _run_ids(frame, columns)
    if ids.empty:
        return pd.DataFrame(columns=["start", "end", "length"])
    dates = ids.index.to_series()
    g = dates.groupby(ids.to_numpy())
    return pd.DataFrame({"start": g.min().to_numpy(), "end": g.max().to_numpy(), "length": g.size().to_numpy()})


def stale_mask(frame: pd.DataFrame, min_run: int = 5, columns=("BID", "ASK")) -> pd.Series:
    """True for observations inside a run of at least ``min_run`` identical quotes."""
    ids = _run_ids(frame, columns)
    if ids.empty:
        return pd.Series(False, index=frame.index)
    length = ids.map(ids.value_counts())
    return (length >= min_run).reindex(frame.index, fill_value=False)


# ---------------------------------------------------------------------------
# Forwards and USD-base inversion


def pip_factor_from_scaling(scaling: float) -> float:
    """Pip factor from LSEG's ScalingFactor (points per unit of the quoted price)."""
    return 1.0 / float(scaling)


def forward_outright(spot, points, pip_factor):
    """Outright forward in the pair's own quotation: spot + points × pip factor."""
    return spot + points * pip_factor


def outright_bid_ask(spot_bid, spot_ask, points_bid, points_ask, pip_factor):
    """Two-sided outright: bid side from the bids, ask side from the asks."""
    return forward_outright(spot_bid, points_bid, pip_factor), forward_outright(spot_ask, points_ask, pip_factor)


def invert_bid_ask(bid, ask):
    """Invert a two-sided quote. The bid of 1/X is 1/ask and the ask is 1/bid."""
    return 1.0 / ask, 1.0 / bid


def usd_per_unit(bid, ask, usd_base: bool):
    """Bid and ask in USD per unit of the non-USD currency (research design, section 3)."""
    return invert_bid_ask(bid, ask) if usd_base else (bid, ask)
