"""Cboe VIX futures (VX) settlement files: loading and price look-up."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load_vx(folder: Path) -> tuple[dict, list]:
    """Settlement series per contract, keyed by expiry.

    Files named VX_<date>.csv carry the expiry in the name; archive files
    (CFE_*) contain only expired contracts and are keyed by their last trading
    date. The archive quotes early contracts at ten times the index level; a
    day-on-day fall by a factor near ten within a contract is taken as that
    rescaling and earlier settlements are divided by ten. Detected rescaling
    dates are returned for the record.
    """
    contracts, rescaled = {}, []
    for f in sorted(folder.glob("*.csv")):
        df = pd.read_csv(f)
        df["date"] = pd.to_datetime(df["Trade Date"], format="mixed")
        s = df.set_index("date")["Settle"].astype(float)
        s = s[s > 0].sort_index()
        s = s[~s.index.duplicated(keep="last")]
        if s.empty:
            continue
        ratio = (s / s.shift(1)).to_numpy()
        jumps = np.where((ratio > 0.07) & (ratio < 0.14))[0]
        if jumps.size:
            k = jumps[0]
            s.iloc[:k] = s.iloc[:k] / 10.0
            rescaled.append((f.name, s.index[k].date().isoformat()))
        key = pd.Timestamp(f.stem.split("_")[1]) if f.stem.startswith("VX_") else s.index[-1]
        contracts[key] = s
    # Contracts that expired before the detected rescaling date are entirely at the old scale.
    if rescaled:
        cut = pd.Timestamp(min(d for _, d in rescaled))
        for key, s in contracts.items():
            if s.index[-1] < cut:
                contracts[key] = s / 10.0
    return dict(sorted(contracts.items())), rescaled


def price_on(s: pd.Series, date) -> float:
    prior = s.loc[:date]
    return float(prior.iloc[-1]) if not prior.empty and (date - prior.index[-1]).days <= 7 else np.nan
