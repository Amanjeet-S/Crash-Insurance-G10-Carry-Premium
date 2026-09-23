# Data audit record

Stage 1 audit of the LSEG panel, retrieval of 23 September 2026. This record states which checks were run, the rules applied and the scripts that reproduce them. It reports no LSEG-derived values, dates or counts. Those are in the private audit (`data/private/audit/2026-09-23/`) until the licence terms for publication are confirmed.

## Reproduction

```bash
.venv-lseg/bin/python scripts/acquire_lseg_fx.py --retrieval-date 2026-09-23
.venv/bin/python scripts/audit_fx_panel.py --retrieval-date 2026-09-23
```

Acquisition needs an LSEG Workspace entitlement. Parsing and sampling rules are in `src/qef/data/panel.py` and are tested on synthetic data in `tests/test_panel.py`.

## Acquisition

- **Instruments.** Those listed in the [data plan](../data_plan.md), plus `USD1MD=`, the one-month USD deposit rate, so that the forward-point sign check compares deposit rates with deposit rates.
- **Request.** Daily interval, all default fields, 1 January 1995 to 22 September 2026, lseg-data 2.1.1.
- **Truncation.** If the first returned date lies more than seven days after the requested start, the earlier range is requested again, until a request returns no rows.
- **Storage.** Returned tables are stored unmodified at full precision. Each request is logged with its RIC, fields, interval, range, UTC time, library version, row count and any error. A manifest records the SHA-256 hash, row count and first and last date of every file.
- **Metadata.** Provider search metadata are recorded for every RIC: document title, first and second currency, underlying RIC, contributor and forward-point scaling factor.

## Definitions

- **Two-sided quote:** bid and ask both present. Crossed and zero-width quotes count as two-sided and are flagged separately.
- **New York month-end:** the last business day of the calendar month under the Federal Reserve holiday calendar. A Saturday holiday is not moved to the Friday.
- **Month-end rule** (research design, section 4). A quote is observed if it is two-sided on the month-end. It is substituted if the last two-sided quote lies within the five preceding business days. Otherwise it is missing.
- **Stale run:** consecutive business days with identical bid and ask. A run of five or more days is flagged.
- **Mid:** (bid + ask)/2, computed from the two sides rather than taken from the vendor field.

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

## Samples

The sample rules of the research design (section 4) are applied as follows. Month-end availability is recorded in `samples.csv`.

- **Primary start.** The first month-end at which all 45 composite one-month smile series (ATM, 25Δ and 10Δ risk reversals and butterflies for nine currencies) and the nine forward-point series are available under the month-end rule. The date without substitution and the binding series are also recorded.
- **Extended sample.** For each currency, the first month-end with Fenics one-month ATM, 25Δ and 10Δ risk reversals and butterflies, and composite forward points. A second version requires only the ATM, 25Δ risk reversal and 25Δ butterfly used in calibration. The extended sample needs at least four currencies (two long, two short).
- **Long ATM sample.** For each currency, the first month-end with ATM volatility (composite or Fenics) and forward points.

## Outcomes that set a method

- **Underlying pair.** Provider metadata give every volatility instrument as the currency against USD, including NOK and SEK. The realised-volatility test agrees for every currency. No triangulation through EURUSD is applied.
- **Forward points.** The provider scaling factors equal the pip factors in `src/qef/fx/conventions.py`, and the implied rate differentials confirm them.
- **Risk reversals.** The evidence supports call minus put on the quoted pair, with 10Δ quotes following the 25Δ convention.

## Not yet done

- Check 6: delta and premium conventions from provider documentation. The configuration remains provisional.
- Check 7: the butterfly convention, from provider documentation and the 10Δ diagnostic.
- Quote revisions, which need a second retrieval.
