"""Month-end smile inputs from the raw LSEG panel.

For each currency and New York month-end the builder returns spot, the
one-month outright forward, the quote-currency one-month rate, discount
factors and the five one-month volatility quotes, each selected under the
five-business-day rule of the research design (section 4).

Timing. Dates are built on the New York business calendar: spot is the
trade date plus the spot lag (one business day for USDCAD, two otherwise),
delivery is spot plus one calendar month under the modified-following rule,
and expiry is delivery minus the spot lag. Volatility time is
(expiry − trade)/365. Discounting from spot to delivery uses DAYS_MAT of the
forward where the provider reports it and the computed spot-to-delivery days
otherwise (``days_source``), on the money-market basis of the quote currency.
Currency-specific holidays are not applied, so dates can differ from market
dates by a business day around non-US holidays.

Completeness. ``complete_calib`` requires spot, forward, rate, ATM, 25Δ risk
reversal and 25Δ butterfly (the calibration set; used for calibration and for
the extended sample). ``has_10d`` records the 10Δ quotes, and ``complete``
requires both (the primary-sample rule).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..fx.conventions import G10
from .panel import NewYorkCalendar, read_raw, sample_month_ends

VOL_BLOCKS = {"atm": ("vol_atm", "O"), "rr25": ("vol_rr25", "RR"), "bf25": ("vol_bf25", "BF"),
              "rr10": ("vol_rr10", "R10"), "bf10": ("vol_bf10", "B10")}
MM_BASIS = {"USD": 360, "JPY": 360, "CHF": 360, "CAD": 365, "NOK": 360, "SEK": 360}
SPOT_LAG = {"CAD": 1}
CALIB_FIELDS = ("spot", "fwd", "rate", "atm", "rr25", "bf25")
# Holidays are generated once; building the calendar inside the offset on every call is slow.
_BDAY = pd.offsets.CustomBusinessDay(holidays=NewYorkCalendar().holidays("1990-01-01", "2040-12-31"))


def option_dates(trade: pd.Timestamp, currency: str, months: int = 1):
    """Spot, delivery and expiry dates on the New York business calendar (trade dates 1990 to 2040)."""
    if not pd.Timestamp("1990-01-01") <= trade <= pd.Timestamp("2040-11-30"):
        raise ValueError("trade date outside the holiday table")
    bday = _BDAY
    lag = SPOT_LAG.get(currency, 2)
    spot = trade + lag * bday
    target = spot + pd.DateOffset(months=months)
    delivery = target if bday.is_on_offset(target) else target + bday
    if delivery.month != target.month:  # modified following
        delivery = target - bday
    expiry = delivery - lag * bday
    return spot, delivery, expiry


def _rate_ric(quote_ccy: str) -> tuple[str, str]:
    return ("USD1MOIS=", "rates") if quote_ccy == "USD" else (f"{quote_ccy}1MD=", "rates")


def _selected(root: Path, block: str, ric: str, month_ends) -> pd.DataFrame:
    path = root / block / f"{ric}.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    return sample_month_ends(read_raw(path), month_ends)


def build_month_end_inputs(raw_root, month_ends, currencies=None, contributor: str = "", tenor: str = "1M") -> pd.DataFrame:
    """One row per (currency, month-end). Volatilities are returned as decimals."""
    root = Path(raw_root)
    currencies = currencies or list(G10)
    frames = []
    for ccy in currencies:
        conv = G10[ccy]
        quote_ccy = ccy if conv.usd_base else "USD"
        spot = _selected(root, "spot", f"{ccy}=", month_ends)
        fwd = _selected(root, "forward", f"{ccy}{tenor}=", month_ends)
        rate_ric, rate_block = _rate_ric(quote_ccy)
        rate = _selected(root, rate_block, rate_ric, month_ends)
        fwd_raw = read_raw(root / "forward" / f"{ccy}{tenor}=.csv")
        fwd_raw = fwd_raw[~fwd_raw.index.duplicated(keep="last")]
        days = fwd_raw["DAYS_MAT"].reindex(fwd["source_date"]).to_numpy() if "DAYS_MAT" in fwd_raw else np.full(len(fwd), np.nan)
        out = pd.DataFrame(index=pd.DatetimeIndex(month_ends, name="date"))
        out["currency"] = ccy
        out["S"] = spot["mid"].to_numpy()
        out["F"] = out["S"] + fwd["mid"].to_numpy() * conv.pip_factor
        dates = [option_dates(t, ccy) for t in out.index]
        computed = np.array([(dl - sp).days for sp, dl, _ in dates], dtype=float)
        reported = np.isfinite(days) & (days > 0)
        out["days"] = np.where(reported, days, computed)
        out["days_source"] = np.where(reported, "DAYS_MAT", "calendar")
        out["tau"] = np.array([(ex - t).days for t, (_, _, ex) in zip(out.index, dates)], dtype=float) / 365.0
        out["r_quote"] = rate["mid"].to_numpy() / 100.0
        out["df_quote"] = 1.0 / (1.0 + out["r_quote"] * out["days"] / MM_BASIS[quote_ccy])
        out["df_base"] = out["F"] * out["df_quote"] / out["S"]
        status = {"spot": spot["status"], "fwd": fwd["status"], "rate": rate["status"]}
        suffix = f"={contributor}" if contributor else "="
        for name, (block, code) in VOL_BLOCKS.items():
            sel = _selected(root, block, f"{ccy}{tenor}{code}{suffix}", month_ends)
            out[name] = sel["mid"].to_numpy() / 100.0
            status[name] = sel["status"]
        st = pd.DataFrame(status, index=out.index)
        out["n_substituted"] = (st == "substituted").sum(axis=1).to_numpy()
        out["complete_calib"] = (st[list(CALIB_FIELDS)] != "missing").all(axis=1).to_numpy()
        out["has_10d"] = (st[["rr10", "bf10"]] != "missing").all(axis=1).to_numpy()
        out["complete"] = out["complete_calib"] & out["has_10d"]
        frames.append(out.reset_index())
    return pd.concat(frames, ignore_index=True)
