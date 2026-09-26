# Data audit record

This record covers the Stage 1 audit of the LSEG panel retrieved on 23 September 2026. It states which checks were run, the rules applied and the scripts that reproduce them. It reports no quote values; the per-instrument dates and counts stay in the private audit (`data/private/audit/2026-09-23/`), and the other reports publish only the pooled sample windows and counts they need.

## Reproduction

```bash
.venv-lseg/bin/python scripts/acquire_lseg_fx.py --retrieval-date 2026-09-23
.venv/bin/python scripts/audit_fx_panel.py --retrieval-date 2026-09-23
```

Acquisition needs an LSEG Workspace entitlement. The parsing and sampling rules are in `src/qef/data/panel.py` and are tested on synthetic data in `tests/test_panel.py`.

## Acquisition

The instruments are those listed in the [data plan](../data_plan.md), together with `USD1MD=`, the one-month USD deposit rate, so that the forward-point sign check compares deposit rates with deposit rates. Each series was requested at a daily interval with all default fields from 1 January 1995 to 22 September 2026, using lseg-data 2.1.1. If the first returned date lay more than seven days after the requested start, the earlier range was requested again until a request returned no rows.

Returned tables are stored unmodified at full precision. Each request is logged with its RIC, fields, interval, range, UTC time, library version, row count and any error, and a manifest records the SHA-256 hash, row count and first and last date of every file. The provider's search metadata are recorded for every RIC: document title, first and second currency, underlying RIC, contributor and forward-point scaling factor.

## Definitions

A quote is two-sided when bid and ask are both present; crossed and zero-width quotes count as two-sided and are flagged separately. A New York month-end is the last business day of the calendar month under the Federal Reserve holiday calendar, with a Saturday holiday not moved to the Friday. Under the month-end rule of the research design (section 4), a quote is observed if it is two-sided on the month-end, substituted if the last two-sided quote lies within the five preceding business days, and missing otherwise. A stale run is a sequence of consecutive business days with identical bid and ask, and a run of five or more days is flagged. The mid is (bid + ask)/2, computed from the two sides rather than taken from the vendor field.

## Checks

| Check | Rule | Output |
| --- | --- | --- |
| 1. Coverage | First and last date with a two-sided quote, per RIC and contributor; weekend, holiday and duplicate rows counted | `coverage.csv` |
| 2. Two-sidedness | Month-ends observed, substituted and missing under the month-end rule, over each series' own window and over the primary and pre-primary windows | `two_sidedness.csv` |
| 3. Staleness | Stale runs per RIC; share of stale days and of stale month-end quotes. The flag rule of the design applies to 25Δ and 10Δ butterflies | `staleness.csv` |
| 4. Underlying pair | First and second currency and underlying RIC from search metadata, for every volatility RIC. Supporting test: RMSE of one-month ATM volatility against realised volatility over the next 21 business days, for the USD pair and for the EUR cross | `conventions.csv` |
| 5. Quote semantics | Units from quote magnitudes. Point scaling from the provider scaling factor, tested by (F/S − 1)·360/days against the one-month deposit differential on days when that differential is at least one percentage point. Risk-reversal sign: the sign of the 25Δ risk reversal against the sign of the correlation between daily changes in log spot and ATM volatility, over the full sample and by calendar year. 10Δ against 25Δ: sign agreement of risk reversals and the ratios of 10Δ to 25Δ quotes | `conventions.csv` |
| 8. Splice | Fenics minus composite mids on common business days and month-ends: mean, standard deviation, share of equal mids, spreads. TIFO is compared with both for the EUR one-month ATM, 25Δ risk reversal and butterfly | `splice.csv` |
| 9. Plausibility | Crossed quotes, zero-width quotes, non-positive ATM volatility and spot, negative butterfly mids, and forward points whose sign differs from the one-month deposit differential when that differential is at least 0.10 percentage points. Flagged, not removed | `plausibility.csv` |

Checks 6 and 7, on the delta and premium conventions and on the butterfly convention, rely on provider documentation and on the calibrated smiles; they are recorded in the research log of 24 September 2026 and in [smile_calibration.md](smile_calibration.md).

## Samples

The sample rules of the research design (section 4) are applied as follows, with month-end availability recorded in `samples.csv`. The primary start is the first month-end at which all 45 composite one-month smile series (ATM, 25Δ and 10Δ risk reversals and butterflies for nine currencies) and the nine forward-point series are available under the month-end rule; the date without substitution and the binding series are also recorded. For the extended sample, the audit records for each currency the first month-end with Fenics one-month ATM, 25Δ and 10Δ risk reversals and butterflies and composite forward points, and a second version that requires only the ATM, 25Δ risk reversal and 25Δ butterfly used in calibration; the extended sample needs at least four currencies, two long and two short. For the long ATM sample, the audit records for each currency the first month-end with ATM volatility (composite or Fenics) and forward points.

## Outcomes that set a method

The provider metadata give every volatility instrument as the currency against USD, including NOK and SEK, and the realised-volatility test agrees for every currency, so no triangulation through EURUSD is applied. The provider scaling factors equal the pip factors in `src/qef/fx/conventions.py`, and the implied rate differentials confirm them. The evidence supports reading risk reversals as call minus put on the quoted pair, with 10Δ quotes following the 25Δ convention. Provider documentation does not define the delta, premium or butterfly conventions, so market conventions apply (check 6), and the market-strangle reading of the butterfly is primary under the pre-registered order of evidence (check 7).

## Quote revisions

On 26 September 2026 I retrieved every instrument again with the same requested range, 1 January 1995 to 22 September 2026, and compared each file with the earlier retrieval that holds it (23 September, or 24 September for the three-month forwards and rates):

```bash
.venv-lseg/bin/python scripts/acquire_lseg_fx.py --retrieval-date 2026-09-26
.venv/bin/python scripts/check_quote_revisions.py --new 2026-09-26 --old 2026-09-23 2026-09-24
```

All 243 instruments with data returned byte-identical files, so no bid, ask or other field was revised, added or removed on any date, and no month-end input changed. The two instruments without data are the one-month OIS series for EUR and SEK, which do not exist. The search metadata of the 226 instruments searched in both retrievals agree in every descriptive field, and the new search covers the three-month instruments, which the 24 September retrieval did not search: the provider describes the three-month USD OIS series used in the three-month variant as the Fed Funds OIS rate. The check can only detect revisions made between my retrievals; a revision made before 23 September would be present in both. The comparison rules are tested on synthetic data in `tests/test_quote_revisions.py`.
