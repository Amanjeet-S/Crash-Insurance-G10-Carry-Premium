# Data plan

This plan specifies what is acquired, from where, how it is stored and how it is audited; the sample rules are in the [research design](research_design.md). Quote values, per-currency month-level series and calibrated parameters are LSEG-derived and stay in `data/private/`; only portfolio-level month series that pass the review of the [publication policy](docs/publication_policy.md) are published, in `data/public/fx_carry_portfolio_series/`. Coverage facts from the audit (sample windows, first available dates, and counts and shares of stale, substituted or missing quotes) describe the dataset without revealing any value, and the reports publish those they need.

## Sources

### LSEG Workspace

Access is through the LSEG Data Library for Python in a desktop session, using the separate environment pinned in `requirements-lseg.lock`. The licence is an academic Workspace licence, and its entitlements were probed on 23 September 2026; the probe records are held privately.

| Block | Instrument pattern (RIC) | Contributors | Fields |
| --- | --- | --- | --- |
| Spot | `EUR=`, `GBP=`, `AUD=`, `NZD=`, `JPY=`, `CHF=`, `CAD=`, `NOK=`, `SEK=` | Composite | Daily bid, ask, mid |
| One-month forward points | `<CCY>1M=` | Composite | Bid, ask |
| At-the-money volatility | `<CCY><T>O=`, T ∈ {1M, 3M} | Composite `=`, Fenics `=FN` | Bid, ask, mid |
| 25Δ risk reversal | `<CCY><T>RR=` | Composite, Fenics | Bid, ask, mid |
| 25Δ butterfly | `<CCY><T>BF=` | Composite, Fenics | Bid, ask, mid |
| 10Δ risk reversal | `<CCY><T>R10=` | Composite, Fenics | Bid, ask, mid |
| 10Δ butterfly | `<CCY><T>B10=` | Composite, Fenics | Bid, ask, mid |
| USD one-month rates | `USD1MOIS=` (Fed Funds OIS), `USDSROIS1M=` (SOFR OIS), `USD1MD=` (deposit) | Composite | Bid, ask |
| Other one-month rates | `<CCY>1MD=` (deposit), `<CCY>1MOIS=` where available (no such RIC exists for EUR or SEK) | Composite | Bid, ask |
| Three-month forward points and rates (robustness variant 2; retrieved 24 September 2026) | `<CCY>3M=`; `USD3MOIS=`, `USDSROIS3M=`, `USD3MD=`; `<CCY>3MD=` | Composite | Bid, ask |
| VIX futures | `VXc1`, `VXc2` | Exchange | Settlement or close |
| One-month 5Δ risk reversal and butterfly (out-of-sample test of the tail hypotheses; retrieved 27 September 2026) | `<CCY>1MR5=FN`, `<CCY>1MB5=FN` | Fenics only | Bid, ask, mid |

Named-broker contributors (`=TIFO`, `=BGCP`, `=TPI`) also exist but are short, sparse or mid-only, so only TIFO is used, as a check over its available period. One-month 5Δ quotes exist only from Fenics (two-sided, from June 2022) and from Tullett Prebon (mid-only, from September 2022); the Fenics quotes are used. WM/Reuters fixings, the S&P 500 index and the Cboe VIX index are not licensed.

The provider does not document the time of day of daily history. Matching daily values against intraday bars in summer and winter 2026 places the composite daily value at about 21:18 UTC throughout the year and the Fenics value at about 17:16 London time; the research log records the method, and month-end alignment uses these times.

### Free public sources

| Source | Use | Licence to record |
| --- | --- | --- |
| Cboe VX historical data: one daily settlement file per contract (`scripts/acquire_cboe_vx.py`) | VIX roll-down factor in E5; checked against LSEG `VXc2` | Redistribution not established: stored in `data/private/cboe/` |
| Verdelhan currency portfolios (web.mit.edu/adrienv, `CurrencyPortfolios.xls`, retrieved 27 September 2026; `scripts/acquire_verdelhan.py`) | Sign and magnitude check of HML_FX | No terms stated: stored in `data/private/verdelhan/` |
| Federal Reserve H.10 and H.15 releases, if needed | Rate and spot cross-checks | Public domain (US government) |

Each download keeps its original bytes and is stored privately under `data/private/cboe/<retrieval-date>/` or `data/private/verdelhan/<retrieval-date>/`, with a `manifest.json` that records the URL, retrieval time and SHA-256 hash of each file (for the Verdelhan file also its size, its HTTP Last-Modified header and the absence of stated terms). `data/public/` holds only this project's own portfolio-level series, in `data/public/fx_carry_portfolio_series/` with a `manifest.json` of SHA-256 hashes and no date subfolder.

## Storage and provenance

```
data/private/lseg/<retrieval-date>/
    raw/<block>/<RIC>.csv            returned tables, unmodified (full precision)
    requests.jsonl                   one line per request: RIC, fields, interval, start, end, time, library version
    metadata.csv                     search metadata per RIC: title, currencies, underlying, scaling
    manifest.json                    file list with SHA-256 hashes, row counts, first and last dates
data/private/audit/<retrieval-date>/
    coverage.csv, two_sidedness.csv, staleness.csv, conventions.csv, splice.csv, plausibility.csv, samples.csv
    audit_report.md                  full audit with counts and dates (restricted)
data/private/results/<retrieval-date>/
    month-end smile inputs, calibrated smiles, estimates, summaries and logs
```

Raw tables are never edited, and no cleaned copy of them is stored. Each script reads the raw tables of an identified retrieval and applies the transformations in code (`src/qef/data/panel.py` and `src/qef/data/smile_inputs.py`), which records them: conversion of percentage quotes to decimals, pip factors and outright forwards, inversion of USD-base quotes, and New York month-end selection under the two-sided and five-business-day rules. Derived outputs, among them the month-end smile inputs written by the calibration, go to `data/private/results/<retrieval-date>/`, and the audit tables to `data/private/audit/<retrieval-date>/`. A missing value is never converted to zero; substitutions follow the rules of the research design and are counted. The retrieval time and the observation date are recorded separately, because a current retrieval does not reconstruct what was knowable historically. That distinction matters only for quote revisions, which `scripts/check_quote_revisions.py` checks by comparing repeated retrievals.

## Audit checks (Stage 1)

1. Coverage: the first and last dates with a two-sided quote, for every instrument and contributor.
2. Two-sidedness: the share of month-ends with both bid and ask, and the substitutions required.
3. Staleness: runs of unchanged quotes per instrument, with butterfly runs of five or more business days flagged.
4. Underlying pair: the document title and underlying of every volatility RIC, confirmed through LSEG search metadata, in particular whether the NOK and SEK instruments are against USD or against EUR.
5. Quote semantics: the units (volatility points, or forward points and their pip factor), the sign convention of risk reversals, and whether 10Δ instruments share the 25Δ conventions.
6. Delta and premium conventions, documented per pair with the source of each convention: provider documentation where available, otherwise Reiswich and Wystup (2012), who report the conventions table of Clark (2011).
7. Butterfly convention: provider documentation first; as supporting evidence only, a diagnostic that calibrates the 25Δ smile under each reading on the primary sample's month-ends and compares the absolute errors in predicting the 10Δ risk reversal and butterfly.
8. Splice: Fenics against composite on their overlap, with the mean and dispersion of the differences per instrument.
9. Plausibility: no crossed quotes, non-negative volatilities, and forward points consistent in sign with the rate differential; an observation that fails a check is flagged, not removed.

## Publication

The public audit record states the checks performed and the rules applied, and the reports publish pooled estimates and the coverage facts described above. Estimates and figures that depend on LSEG inputs are generated into `data/private/` and reviewed before any pooled summary is published. The acquisition script lets a reader with their own entitlement reproduce the private layer. On 1 October 2026 the licence holder confirmed that manipulated or transformed data may be published and that the raw data may not. I read a transformation as publishable only if the raw quotes cannot be recovered from it, so calibrated smile parameters, which reproduce the quotes they were fitted to, and per-currency inputs stay private.
