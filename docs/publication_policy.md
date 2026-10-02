# Publication policy

## 1. Research findings and source observations

The LSEG data were obtained under a university student licence. The licence holder advised on 24 September 2026 that they may be used for individual study and research and may not be redistributed, and confirmed on 1 October 2026 that manipulated or transformed data may be published while the raw data may not. I apply the test of LSEG's redistribution guidance for derived data: a released output must be unrecognisable, non-reversible and impossible to trace back to the original content without exceptional effort.

| Material | Repository treatment |
| --- | --- |
| Raw provider exports, quote histories and saved responses | Retain privately |
| Cleaned subsets, rounded or rescaled quotes, per-currency inputs, calibrated smile parameters and any other reversible transformation | Retain privately |
| Plain averages of raw quotes in quote units | Retain privately |
| Portfolio-level series that average several currencies whose identities are not published, with fewer values per month than the quotes they depend on | Eligible after the reversibility review in section 2 |
| Estimates, test statistics, intervals, coverage facts and written findings | Eligible after checking that they do not disclose or permit reconstruction of observations |
| Original methods and source code | Eligible after checking attribution and embedded data |

## 2. Review of a release

I review the files together, because a series may be harmless alone while another file supplies what is needed to recover the observations. For the portfolio series, `scripts/check_public_reversibility.py` assumes an adversary who knows everything about a month except the smile quotes, counts the published values against the unknown quotes, and constructs, month by month, quote values that differ materially from the true ones and reproduce every published value; the [public series README](../data/public/fx_carry_portfolio_series/README.md) reports the result. The review withdrew a plain average of risk-reversal quotes; the option payoff and premium of the long at-the-money sample, which in some months depend on one currency alone; and three secondary skew-price series, with which a worst-case attacker could have narrowed single risk-reversal quotes. It also checks that the quote values on the provider's quote grid that reproduce a month's published numbers remain numerous, which for the coarser Fenics grid required rounding the extended sample to four significant figures.

A release states what is available, what remains private and what access is required for full reproduction.

## 3. Licences and attribution

Original software is released under [MIT](../LICENSE). Original research writing and the portfolio-level series are released under [CC BY 4.0](../LICENSES/CC-BY-4.0.txt), within the [licensing scope](../LICENSING.md). These licences cover the original released material and confer no rights to the underlying LSEG observations.
