# Data plan

This plan specifies what is acquired, from where, how it is stored and how it is audited. Sample rules are in the [research design](research_design.md). Coverage dates and counts that result from the audit are LSEG-derived. They are kept in `data/private/` until the licence terms for publication are confirmed.

## Sources

### LSEG Workspace

Access is through the LSEG Data Library for Python (desktop session), using the separate environment pinned in `requirements-lseg.lock`. The licence is an academic Workspace licence. Its entitlements were probed on 23 September 2026; the probe records are held privately.

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
| VIX futures | `VXc1`, `VXc2` | Exchange | Settlement or close |

**Also available:**

- **Additional contributors.** Named-broker contributors (`=TIFO`, `=BGCP`, `=TPI`) exist but are short, sparse or mid-only. Only TIFO is used, as a check over its available period.
- **Not licensed.** WM/Reuters fixings, the S&P 500 index and the Cboe VIX index.

### Free public sources

| Source | Use | Licence to record |
| --- | --- | --- |
| Cboe VX historical data: one daily settlement file per contract (`scripts/acquire_cboe_vx.py`) | VIX roll-down factor in E5; checked against LSEG `VXc2` | Redistribution not established: stored in `data/private/cboe/` |
| Verdelhan currency portfolios (web.mit.edu/adrienv) | Sign and magnitude check of HML_FX | Author terms |
| Federal Reserve H.10 and H.15 releases, if needed | Rate and spot cross-checks | Public domain (US government) |

Each public snapshot is stored under `data/public/<source>/<date>/` with its original bytes, URL, retrieval time, SHA-256 hash and licence note.

**Snapshot time.** The provider does not document the time of day of daily history. Matching daily values against intraday bars (summer and winter 2026) places the composite daily value at about 21:18 UTC throughout the year and the Fenics value at about 17:16 London time. The research log records the method; month-end alignment uses these times.

## Storage and provenance

```
data/private/lseg/<retrieval-date>/
    raw/<block>/<RIC>.csv            returned tables, unmodified (full precision)
    requests.jsonl                   one line per request: RIC, fields, interval, start, end, time, library version
    manifest.json                    file list with SHA-256 hashes, row counts, first and last dates
data/private/audit/<retrieval-date>/
    coverage.csv, two_sidedness.csv, staleness.csv, conventions.csv
    audit_report.md                  full audit with counts and dates (restricted)
```

**Rules.**

- **Raw tables are never edited.** Cleaning is done by code into `data/private/clean/`. Every transformation is recorded: unit conversion, inversion of USD-base quotes, pip factors and month-end selection.
- **Missing values stay missing.** A missing value is never converted to zero. Substitutions follow the rules in the research design and are counted.
- **Recording times.** The retrieval time and the observation date are recorded separately. A current retrieval does not reconstruct what was knowable historically. That matters only for quote revisions, which the audit checks by comparing repeated pulls.

## Audit checks (Stage 1)

1. **Coverage.** First and last dates with a two-sided quote, for every instrument and contributor.
2. **Two-sidedness.** The share of month-ends with both bid and ask, and the substitutions required.
3. **Staleness.** Runs of unchanged quotes, per instrument. Butterfly runs of five or more business days are flagged.
4. **Underlying pair.** The document title and underlying of every volatility RIC, confirmed through LSEG search metadata. Of particular concern is whether the NOK and SEK instruments are against USD or EUR.
5. **Quote semantics.**
   - Units: volatility points, or forward points and their pip factor.
   - The sign convention of risk reversals.
   - Whether 10Δ instruments share the 25Δ conventions.
6. **Delta and premium conventions.** Documented per pair, with the source of each convention (provider documentation where available, otherwise Reiswich and Wystup, 2012, which reports the conventions table of Clark, 2011).
7. **Butterfly convention.**
   - Provider documentation first.
   - The supporting diagnostic calibrates the 25Δ smile under each reading and compares the absolute errors in predicting the 10Δ risk reversal and butterfly. It uses the primary sample's month-ends.
   - The comparison is recorded as supporting evidence only.
8. **Splice.** Fenics against composite on their overlap: mean and dispersion of differences per instrument.
9. **Plausibility.**
   - No crossed quotes.
   - Non-negative volatilities.
   - Forward points consistent in sign with the rate differential.
   - Where a check fails, the observation is flagged, not removed.

## Publication handling

- **Public summary.** The public audit summary states the checks performed and the rules applied. It does not report LSEG-derived values, counts or dates until they are cleared.
- **Figures and statistics.** Figures and estimates that depend on LSEG inputs are generated into `data/private/` or `results/private/` until they are reviewed.
- **Reproduction.** The acquisition script lets a reader with their own entitlement reproduce the private layer.
