"""Quotation and delta conventions for the G10 currency set.

Market conventions follow Reiswich and Wystup (2012), which reports the
conventions table of Clark (2011). NZDUSD, USDNOK and USDSEK are not in that
table; they follow its general rule (premium-adjusted delta when the premium is
paid in the base currency). The audit confirmed each instrument's underlying
pair and forward-point scaling; provider documentation does not state the delta
convention (research log, 24 September 2026).
"""

from __future__ import annotations

from dataclasses import dataclass

from .gk import DeltaConvention


@dataclass(frozen=True)
class PairConvention:
    currency: str  # the non-USD currency
    pair: str  # market quotation, base then quote
    usd_base: bool  # True if quoted as foreign units per USD
    pip_factor: float  # forward outright = spot + points * pip_factor
    delta: DeltaConvention

    def to_usd_per_unit(self, quote):
        """Convert a market quote to USD per unit of the non-USD currency."""
        return 1.0 / quote if self.usd_base else quote


_PIPS = DeltaConvention(spot=True, premium_adjusted=False)
_PA = DeltaConvention(spot=True, premium_adjusted=True)

G10 = {
    c.currency: c
    for c in (
        PairConvention("EUR", "EURUSD", False, 1e-4, _PIPS),
        PairConvention("GBP", "GBPUSD", False, 1e-4, _PIPS),
        PairConvention("AUD", "AUDUSD", False, 1e-4, _PIPS),
        PairConvention("NZD", "NZDUSD", False, 1e-4, _PIPS),
        PairConvention("JPY", "USDJPY", True, 1e-2, _PA),
        PairConvention("CHF", "USDCHF", True, 1e-4, _PA),
        PairConvention("CAD", "USDCAD", True, 1e-4, _PA),
        PairConvention("NOK", "USDNOK", True, 1e-4, _PA),
        PairConvention("SEK", "USDSEK", True, 1e-4, _PA),
    )
}
