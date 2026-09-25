"""Month-end smile inputs from the raw LSEG panel.

For each currency and New York month-end the builder returns spot, the
outright forward, the quote-currency rate, discount factors and the five
volatility quotes at one tenor, each selected under the five-business-day
rule of the research design (section 4). The tenor is one month (primary) or
three months (robustness variant 2 of the design, section 8).

Timing. Dates are built on the New York business calendar: spot is the
trade date plus the spot lag (one business day for USDCAD, two otherwise),
delivery is spot plus the tenor in calendar months under the
modified-following rule, and expiry is delivery minus the spot lag. Volatility
time is (expiry − trade)/365. Discounting from spot to delivery uses DAYS_MAT
of the forward where the provider reports it and the computed spot-to-delivery
days otherwise (``days_source``), on the money-market basis of the quote
currency. Currency-specific holidays are not applied, so dates can differ from
market dates by a business day around non-US holidays.

Rates. The quote-currency rate has the option's tenor: the USD OIS rate
(USD1MOIS=, USD3MOIS=) when the quote currency is USD, and the deposit rate
(<CCY>1MD=, <CCY>3MD=) otherwise, as recorded in the research log of
24 September 2026 (the deviation on discount rates, and the choices fixed for
the three-month variant). The provider metadata of the 23 September
retrieval describe USD1MOIS= as the Fed Funds OIS; the 24 September
retrieval, which holds USD3MOIS=, carries no metadata.

Completeness. ``complete_calib`` requires spot, forward, rate, ATM, 25Δ risk
reversal and 25Δ butterfly (the calibration set; used for calibration and for
the extended sample). ``has_10d`` records the 10Δ quotes, and ``complete``
requires both (the primary-sample rule).
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from ..fx.conventions import G10
from .panel import NewYorkCalendar, read_raw, sample_month_ends

VOL_BLOCKS = {"atm": ("vol_atm", "O"), "rr25": ("vol_rr25", "RR"), "bf25": ("vol_bf25", "BF"),
              "rr10": ("vol_rr10", "R10"), "bf10": ("vol_bf10", "B10")}
MM_BASIS = {"USD": 360, "JPY": 360, "CHF": 360, "CAD": 365, "NOK": 360, "SEK": 360}
SPOT_LAG = {"CAD": 1}
TENOR_MONTHS = {"1M": 1, "3M": 3}
CALIB_FIELDS = ("spot", "fwd", "rate", "atm", "rr25", "bf25")
QUOTE_FIELDS = CALIB_FIELDS + ("rr10", "bf10")
_CAL_START, _CAL_END = pd.Timestamp("1990-01-01"), pd.Timestamp("2040-12-31")
# Holidays are generated once; building the calendar inside the offset on every call is slow.
_BDAY = pd.offsets.CustomBusinessDay(holidays=NewYorkCalendar().holidays(_CAL_START, _CAL_END))


def option_dates(trade: pd.Timestamp, currency: str, months: int = 1):
    """Spot, delivery and expiry dates on the New York business calendar.

    The holiday table runs from 1990 to 31 December 2040. A trade date before
    it, or a delivery date after its end, raises: past the end, holidays are
    unknown, so a rolled delivery or the expiry could fall on one.
    """
    if trade < _CAL_START:
        raise ValueError("trade date outside the holiday table")
    bday = _BDAY
    lag = SPOT_LAG.get(currency, 2)
    spot = trade + lag * bday
    target = spot + pd.DateOffset(months=months)
    delivery = target if bday.is_on_offset(target) else target + bday
    if delivery.month != target.month:  # modified following
        delivery = target - bday
    if delivery > _CAL_END:
        raise ValueError("delivery date outside the holiday table")
    expiry = delivery - lag * bday
    return spot, delivery, expiry


def tenor_months(tenor: str) -> int:
    """Calendar months of a tenor label ("1M" or "3M")."""
    if tenor not in TENOR_MONTHS:
        raise ValueError(f"tenor must be one of {sorted(TENOR_MONTHS)}")
    return TENOR_MONTHS[tenor]


def rate_ric(quote_ccy: str, tenor: str = "1M") -> tuple[str, str]:
    """(RIC, block) of the quote-currency rate at the tenor: USD OIS for USD, the deposit rate otherwise."""
    tenor_months(tenor)
    return (f"USD{tenor}OIS=", "rates") if quote_ccy == "USD" else (f"{quote_ccy}{tenor}D=", "rates")


def input_rics(currency: str, contributor: str = "", tenor: str = "1M") -> dict:
    """Field name -> (block, RIC) of every raw series behind one currency's inputs."""
    quote_ccy = currency if G10[currency].usd_base else "USD"
    ric, block = rate_ric(quote_ccy, tenor)
    suffix = f"={contributor}" if contributor else "="
    out = {"spot": ("spot", f"{currency}="), "fwd": ("forward", f"{currency}{tenor}="), "rate": (block, ric)}
    out.update({name: (blk, f"{currency}{tenor}{code}{suffix}") for name, (blk, code) in VOL_BLOCKS.items()})
    return out


def locate_raw(roots: Sequence, block: str, ric: str) -> Path:
    """Path of a RIC's raw table in the first of ``roots`` that holds it."""
    roots = [Path(r) for r in roots]
    for r in roots:
        path = r / block / f"{ric}.csv"
        if path.exists():
            return path
    raise FileNotFoundError(f"{block}/{ric}.csv in none of {[str(r) for r in roots]}")


def _selected(path: Path, month_ends) -> pd.DataFrame:
    return sample_month_ends(read_raw(path), month_ends)


def build_month_end_inputs(raw_root, month_ends, currencies=None, contributor: str = "", tenor: str = "1M",
                           extra_roots: Sequence = (), quote_status: bool = False) -> pd.DataFrame:
    """One row per (currency, month-end). Volatilities are returned as decimals.

    Provenance: each series is read whole from the first of ``raw_root``,
    ``*extra_roots`` that holds its file, so a series is never spliced across
    retrievals. In the project, spot, all volatility quotes and the one-month
    forwards and rates come from the retrieval of 23 September 2026; the
    three-month forward points and rates are absent there and come from the
    retrieval of 24 September 2026 (its manifest.json lists the files).
    ``locate_raw`` reports the file used for any series.

    ``quote_status`` adds one ``status_<field>`` column per quote (observed,
    substituted or missing); the default leaves the columns unchanged.
    """
    roots = [Path(raw_root), *(Path(r) for r in extra_roots)]
    months = tenor_months(tenor)
    currencies = currencies or list(G10)
    frames = []
    for ccy in currencies:
        conv = G10[ccy]
        quote_ccy = ccy if conv.usd_base else "USD"
        paths = {name: locate_raw(roots, block, ric) for name, (block, ric) in input_rics(ccy, contributor, tenor).items()}
        spot = _selected(paths["spot"], month_ends)
        fwd = _selected(paths["fwd"], month_ends)
        rate = _selected(paths["rate"], month_ends)
        fwd_raw = read_raw(paths["fwd"])
        fwd_raw = fwd_raw[~fwd_raw.index.duplicated(keep="last")]
        days = fwd_raw["DAYS_MAT"].reindex(fwd["source_date"]).to_numpy() if "DAYS_MAT" in fwd_raw else np.full(len(fwd), np.nan)
        out = pd.DataFrame(index=pd.DatetimeIndex(month_ends, name="date"))
        out["currency"] = ccy
        out["S"] = spot["mid"].to_numpy()
        out["F"] = out["S"] + fwd["mid"].to_numpy() * conv.pip_factor
        dates = [option_dates(t, ccy, months) for t in out.index]
        computed = np.array([(dl - sp).days for sp, dl, _ in dates], dtype=float)
        reported = np.isfinite(days) & (days > 0)
        out["days"] = np.where(reported, days, computed)
        out["days_source"] = np.where(reported, "DAYS_MAT", "calendar")
        out["tau"] = np.array([(ex - t).days for t, (_, _, ex) in zip(out.index, dates)], dtype=float) / 365.0
        out["r_quote"] = rate["mid"].to_numpy() / 100.0
        out["df_quote"] = 1.0 / (1.0 + out["r_quote"] * out["days"] / MM_BASIS[quote_ccy])
        out["df_base"] = out["F"] * out["df_quote"] / out["S"]
        status = {"spot": spot["status"], "fwd": fwd["status"], "rate": rate["status"]}
        for name in VOL_BLOCKS:
            sel = _selected(paths[name], month_ends)
            out[name] = sel["mid"].to_numpy() / 100.0
            status[name] = sel["status"]
        st = pd.DataFrame(status, index=out.index)
        out["n_substituted"] = (st == "substituted").sum(axis=1).to_numpy()
        out["complete_calib"] = (st[list(CALIB_FIELDS)] != "missing").all(axis=1).to_numpy()
        out["has_10d"] = (st[["rr10", "bf10"]] != "missing").all(axis=1).to_numpy()
        out["complete"] = out["complete_calib"] & out["has_10d"]
        if quote_status:
            for name in QUOTE_FIELDS:
                out[f"status_{name}"] = st[name].to_numpy()
        frames.append(out.reset_index())
    return pd.concat(frames, ignore_index=True)
