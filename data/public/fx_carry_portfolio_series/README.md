# Portfolio-level month series: crash insurance and the G10 carry premium

Amanjeet Singh. This folder holds the month-level portfolio series behind the portfolio-level results of my paper
"Crash insurance and the G10 carry premium" (`paper/main.tex`).
With them, anyone can rerun the paper's portfolio-level tests without an LSEG licence, apart from E4
and the 25-delta comparisons, whose series are withheld:
`scripts/reproduce_from_public.py` recomputes those results from these files alone and compares each
with the number printed in the paper.

Every value is a portfolio-level quantity: an average over the legs of a carry portfolio (six legs in
the primary sample, four in the extended sample, four or six in the long at-the-money sample), an
average over the nine currencies (`sigma_fx`), or a ratio or sum of such averages. No value refers to a
single currency, and no file says which currencies were held.

## Licence and provenance

The series are derived from quotes that I obtained from LSEG Workspace under a student licence provided
by my university. They are transformed data. Each value is computed from the quotes of four to nine
currencies at once and averaged over the legs of a portfolio whose composition is not published. Most
values are nonlinear functions of calibrated option smiles, forward rates and realised spot rates (option
premia and payoffs, hedged returns, hedge-cost terms and moment bounds); `FD` is an average of log forward
discounts with signs set by the unpublished legs, and `sigma_fx` an average, over the nine currencies
and the business days of a return window, of absolute daily log spot changes. My university, the
licence holder, confirmed on 1 October 2026 that manipulated or transformed data may be published and
that raw data may not (research log, 1 October 2026). I publish these series under that confirmation.
They contain no raw quote, no per-currency value, no leg identity and no calibrated parameter. The raw
LSEG data may not be redistributed, and anyone who wants to rebuild these series from the quotes needs
their own LSEG Workspace licence.

I release the series under [Creative Commons Attribution 4.0 International](../../../LICENSES/CC-BY-4.0.txt)
as my original research results; the grant covers my contribution to them and confers no right to the
underlying LSEG observations (`LICENSING.md`). For attribution, cite Amanjeet Singh, "Crash insurance
and the G10 carry premium", and this repository.

I wrote the files with `scripts/export_public_series.py` from my private results of the retrieval of
23 September 2026. `manifest.json` gives the SHA-256 hash, the size, the number
of rows, the columns and the first and last date of each CSV file.

## Files

| File | Sample | Quotes | Legs per side | Rows | Formation month-ends | Rows with status ok |
| --- | --- | --- | --- | --- | --- | --- |
| `primary.csv` | primary | composite | 3 | 160 | 2013-05-31 to 2026-08-31 | 159 |
| `extended.csv` | extended | Fenics | 2 | 76 | 2007-01-31 to 2013-04-30 | 72 |
| `long_atm.csv` | long at-the-money | composite (spot, forward points, ATM volatility) | 2 or 3 | 379 | 1995-01-31 to 2026-07-31 | 379 |

All three files are comma-separated with a header row. Numbers are rounded to six significant figures in
`primary.csv` and `long_atm.csv` and to four in `extended.csv` (section on traceability below);
a missing value is an empty field. Rows are sorted by date.

## Samples, regimes and dates

- `date`: the formation month-end t, the last New York business day of the calendar month (Federal
  Reserve holiday calendar), written YYYY-MM-DD. Quotes are sampled at t under the five-business-day
  substitution rule of the research design. A return dated t runs from t to the expiry of the one-month
  option traded at t, about one month later, and is evaluated at the end-of-day spot on that expiry date
  (`estimate_stage4.month_rows` and `estimate_stage4.spot_on`; `qef.data.smile_inputs.option_dates`).
- Primary sample: composite quotes, three long and three short legs, 2013-05-31 to 2026-08-31: 160 month-ends,
  104 in the zero-rate regime and 56 in the hiking regime. The return window of the last
  month-end ends after the last spot observation of the retrieval, so 159 month-ends have realised returns.
- Extended sample: Fenics quotes, two long and two short legs, 2007-01-31 to 2013-04-30: 76 month-ends, of which
  72 have the four currencies with the calibration set that the sample needs. Fenics quotes before 2010
  update infrequently, butterflies in particular, so this sample is secondary evidence.
- Long at-the-money sample: composite spot, one-month forward points and one-month ATM volatility only,
  1995-01-31 to 2026-07-31: 379 month-ends with realised returns (36 with two legs per
  side, 343 with three). The portfolio holds three and three legs when at least six
  currencies are available and two and two when four or five are (`scripts/estimate_long_atm.py`, rules 1
  to 3). The month-end of August 2026, whose return is not yet realised, is not included.
- `regime`: `zero_rate` for formation month-ends before 1 January 2022 and `hiking` from that date
  (`estimate_e1.REGIME_BREAK`). The labels follow the research design and denote calendar periods, not
  rate paths. Every month-end of the extended sample falls before the break, so its label is `zero_rate`
  throughout; the paper's regime comparisons use the primary sample only.

## Units

Premia, returns, hedge-cost terms, `C_skew` and `FD` are decimal fractions per USD of forward notional,
not basis points: 0.0017 is 17 bp. Premia are carried to the forward's delivery date, so they are in USD
per USD of forward notional at delivery (paper, Section 4.2). `phi` is dimensionless. The variance bounds `var_lo` and `var_hi` are in squared
log-return units (the variance of the one-month log return, not annualised); the oriented skewness bounds
`oskew_lo` and `oskew_hi` are dimensionless. `sigma_fx` is a mean absolute daily log change.

## Columns of primary.csv

Each leg buys a one-month protective option that pays when the leg loses: a put on the currency for a
long leg and a call on the currency for a short leg (for the dollar-base pairs, the corresponding call or
put on USD), with notional equal to the forward notional. Strikes are found from deltas in each pair's
own convention on the calibrated SABR smile (beta = 1), and premia are Garman-Kohlhagen prices in forward
form at the stated volatility (`qef.fx.crash`, `qef.fx.gk`, `qef.fx.sabr`).

| Column | Definition | Code |
| --- | --- | --- |
| `date`, `regime` | As above. | |
| `status` | `ok`: every column is present. `missing_spot`: the spot at the option expiry is not observed in the retrieval (the month-end of August 2026, whose return window ends after the retrieval), so `U`, `H10`, `H25`, `Hatm`, `c1`, `c2`, `c3` and `sigma_fx` are empty; the E1 columns, `FD` and the moment predictors are present. Every primary month-end has E1 status ok. | `estimate_stage4.month_rows` |
| `C_skew_market_10d` | The ex-ante skew price of protection, (1/3) times the sum over the six legs of V_smile minus V_flat, where V_smile is the premium of the leg's protective option at its 10-delta strike on the calibrated smile (market-strangle reading of the 25-delta butterfly, the primary reading of E1), priced at the smile volatility, and V_flat is the premium at the same strike priced at the ATM volatility. | `estimate_e1.e1_series`, `qef.fx.crash.leg_skew_cost` |
| `phi_market_10d` | phi = C_skew / FD, the skew price per unit of forward discount: the E1 estimand. | `estimate_e1.e1_series` |
| `FD` | Forward-discount spread, (1/3) times (the sum of fd over the long legs minus the sum over the short legs), with fd = log(X/F), X the spot and F the one-month outright forward, both as USD per unit of the currency. The legs are the three highest and the three lowest fd at t. | `qef.fx.crash.forward_discount`, `rank_legs` |
| `U` | Unhedged carry return HML^U: (1/3) times the sum over the legs of the signed forward return, rx = X_T/F - 1 for a long leg and minus that for a short leg, with X_T the spot on the option expiry date. | `estimate_stage4.month_rows`, `qef.fx.crash.forward_return` |
| `H10`, `H25`, `Hatm` | Hedged carry return HML^H: U plus (1/3) times the sum over the legs of the option payoff less its premium. `H10` and `H25` use the smile 10-delta and 25-delta strikes and the smile premium; `Hatm` uses the delta-neutral-straddle strike and the premium at the ATM volatility. | `estimate_stage4.month_rows` |
| `c1` | E3 term (i), payoff: (1/3) times the sum over the legs of the 10-delta payoff less the Garman-Kohlhagen premium at the realised volatility, the annualised standard deviation of daily log spot changes over the 21 New York business days ending at t. | `estimate_stage4.month_rows`, `realised_vol` |
| `c2` | E3 term (ii), volatility level: minus (1/3) times the sum over the legs of V_flat less that realised-volatility premium, at the 10-delta strike. | `estimate_stage4.month_rows` |
| `c3` | E3 term (iii), skew: minus (1/3) times the sum over the legs of V_smile less V_flat at the 10-delta strike. | `estimate_stage4.month_rows` |
| `var_lo`, `var_hi` | Lower and upper endpoints of the risk-neutral variance of y = ln(X_T/F) under the USD forward measure, from contracts spanned by the calibrated smile between its 10-delta strikes and bounded beyond them by result R4 under tail setting (i) (boundary elasticities), averaged endpoint by endpoint over the six legs. | `estimate_moments.portfolio_predictors`, `qef.fx.moments.implied_moment_intervals` |
| `oskew_lo`, `oskew_hi` | Endpoints of the oriented skewness under the same setting: the mean over the legs of the skewness interval of a long leg and of minus the skewness interval of a short leg (a short leg's interval is [-skew_hi, -skew_lo]), averaged endpoint by endpoint. | `estimate_moments.portfolio_predictors` |
| `sigma_fx` | Global FX-volatility level over the return window: the mean over the window's business days (after t, up to the one-month expiry of EURUSD) of the cross-sectional mean absolute daily log spot change of the nine currencies against USD (those with a spot on the day). It uses daily spot only. The E5 factor is the residual of an AR(1) fitted by least squares to this series over the months with status `ok` of `extended.csv` and `primary.csv` together, in date order. | `estimate_e5.fx_abs_returns`, `estimate_e5.main` |

## Columns of extended.csv

The same definitions as `primary.csv`, with two long and two short legs (weights 1/2 in place of 1/3),
Fenics smiles and the market-strangle reading at 10 delta only.

| Column | Definition | Code |
| --- | --- | --- |
| `date`, `regime` | As above; `regime` is `zero_rate` in every row. | |
| `status` | `ok`: every column is present. `too_few_currencies`: fewer than four currencies have the calibration set at the month-end (4 month-ends), so every other column is empty. | `estimate_e1.e1_series`, `estimate_stage4.month_rows` |
| `C_skew`, `phi`, `FD` | As `C_skew_market_10d`, `phi_market_10d` and `FD` in `primary.csv`. | `estimate_e1.e1_series` |
| `U`, `H10`, `H25`, `Hatm`, `c1`, `c2`, `c3` | As in `primary.csv`. | `estimate_stage4.month_rows` |
| `sigma_fx` | As in `primary.csv`, over the return window of the extended month-end. | `estimate_e5.fx_abs_returns`, `estimate_e5.main` |

## Columns of long_atm.csv

| Column | Definition | Code |
| --- | --- | --- |
| `date`, `status` | As above; `status` is `ok` in every row. | `estimate_long_atm.month_rows` |
| `n_legs` | Legs per side, 2 or 3 (weights 1/n). | `estimate_long_atm.n_legs` |
| `U` | Unhedged carry return HML^U, as in `primary.csv`. | `estimate_long_atm.month_rows` |
| `Hatm` | ATM-hedged return: U plus (1/n) times the sum over the legs of the payoff less the premium of the protective option at the delta-neutral-straddle strike, the premium being the Garman-Kohlhagen premium at the quoted ATM volatility, carried to delivery. | `estimate_long_atm.month_rows` |

## Exact identities

Some columns are exact functions of others and add no information; I publish them for convenience.
Before rounding, `phi` = `C_skew` / `FD` in both files; `c3` = -`C_skew_market_10d` in
`primary.csv` and `c3` = -`C_skew` in `extended.csv` (to 1e-16); and `H10` = `U` + `c1` + `c2` + `c3` in
both files (to 3e-16). After rounding these hold to the rounding error of each file's significant figures.

## The diffusive null of the 10-delta hedge

The reference theta0 of result R6 for the 10-delta hedge is, for each month, the average over the six
legs of the absolute forward delta of the hedging options (`estimate_stage4.month_rows`). The paper uses
only its mean over the 159 return months of the primary sample, so only the mean is
published:

    theta0_mean_10d = 0.0999517

## What is not published, and why

- The skew price and phi of the market-strangle reading at 25 delta and of the smile-strangle reading at
  10 and 25 delta (`C_skew_market_25d`, `phi_market_25d`, `C_skew_smile_10d`, `phi_smile_10d`,
  `C_skew_smile_25d`, `phi_smile_25d`: E4 and the 25-delta comparison of E1). Each is another average over the same legs of a
  function of the same quotes, and together with the other columns they would let an adversary who knew
  the legs and the exact forwards narrow single risk-reversal quotes; without them that adversary cannot
  (section on traceability below).
- The oriented 10-delta risk reversal `rr_or`, the secondary E2 predictor of the paper: it is a plain
  signed average of six raw 10-delta risk-reversal quotes, in quote units, and so the least transformed
  of the portfolio-level series.
- The payoff and the premium of the at-the-money hedge in `long_atm.csv`, which `Hatm` combines: in
  42 of the 379 months exactly one leg's at-the-money option ends in the money, so
  the published payoff would depend on that leg alone. `Hatm` - `U`, the payoff less the premium, always
  mixes all legs, because every leg pays a premium.
- Leg identities (which currencies are long and short at each month-end), currency counts other than
  `n_legs`, and any per-currency quantity: with them a reader could recover or closely approximate
  individual quotes.
- Calibrated smile parameters and the smile inputs: a calibrated smile reproduces the quotes it was fitted
  to, so its parameters are not a substantive transformation.
- theta0 for each month: the paper uses only its mean, given above.
- Counts of substituted quotes and other audit detail: they enter no test.
- The VIX roll-down factor `r_vix`: it comes from Cboe VX futures settlements, whose redistribution terms
  are not established. A reader downloads them with `scripts/acquire_cboe_vx.py`, and
  `scripts/reproduce_from_public.py` recomputes the factor.
- The series of the robustness grid and of the three-month tenor: the base, 25-delta and ATM rows of the
  robustness grid can be recomputed from `primary.csv`; the other variants need the quotes.
- Verdelhan's currency portfolios: they are downloaded with `scripts/acquire_verdelhan.py`.

## Reproduction

From the root of a clone, after the installation in the repository's `README.md`:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/reproduce_from_public.py \
        --out reproduction_report.md \
        [--cboe-dir data/private/cboe/<retrieval date>/raw] \
        [--verdelhan-file data/private/verdelhan/<retrieval date>/CurrencyPortfolios.xls]

By default the report goes to the git-ignored `data/private/results/public_reproduction.md`. The
script checks the hashes in `manifest.json`, recomputes the paper's portfolio-level results with
the project's own estimation functions, compares each number with the paper at its printed precision
(and, where my private results exist, with my private summaries) and exits with a non-zero status if any
number that these files determine does not match. E5 needs the Cboe files, and the comparison with
published portfolios needs Verdelhan's file and the `xlrd` package (not in `requirements.lock`); each
is used only if present. `tests/test_public_series.py`
checks the files and a few key numbers in the test suite.

The E5 regressions of both samples are reproduced: their FX-volatility innovation is the residual of an
AR(1) fitted to `sigma_fx` over the months with status `ok` of `extended.csv` and `primary.csv` together,
and the VIX roll-down is recomputed from the Cboe files.

Some results cannot be recomputed from these files: everything that needs a single currency's smile or
quotes (the identification of option-implied moments, the per-currency moment intervals and tail exponents,
the out-of-sample test with 5-delta quotes, the leg-by-leg comparison of contributors, the SOFR and
named-broker checks, the calibration diagnostics and the theta0 of the at-the-money hedge), E4 and the
25-delta rows of E1 (their skew prices are withheld), the robustness variants other than the base,
25-delta and ATM rows (of the 25-delta row only the returns and theta_UB), and the three-month tenor. Every result that uses the
oriented risk reversal is also among them, because `rr_or` is not published: its rows of Table 2 in sample
and out of sample, its persistence and its correlations with `phi` and with the return innovations
(Section 5.2), and its change between the regimes (Section 5.1). So are the mean payoff and mean premium
of the at-the-money hedge in 2008 (Section 5.3), because the long sample's payoff and premium are not
published; the long sample's `U` and `Hatm`, and everything the paper builds from them, are reproduced.
The model validation of the moment code is synthetic and already public
(`scripts/validate_moments_models.py`).
`scripts/reproduce_from_public.py` lists each case and its reason.

To rebuild these files from the quotes, a reader with their own LSEG licence runs `scripts/reproduce.sh`
and then `scripts/export_public_series.py`.

## Why the series cannot be traced back to the quotes

LSEG's redistribution guidance treats an output as derived data if it is unrecognisable, non-reversible
and cannot be traced back to the original content without exceptional effort. I checked these files
against that test with `scripts/check_public_reversibility.py`, which reads my private data and reports
pooled statistics only. Its adversary is stronger than any reader: it knows everything about a month
except the option quotes, including the unpublished leg identities and the exact LSEG spot, forwards and
rates.

- Unrecognisable. No column is a quote, the price of a single option or a value of a vendor series.
  Every value averages at least four currencies (nine for `sigma_fx`) whose identities are not published,
  and most are nonlinear functions of calibrated smiles, forwards and realised spot rates.
- Non-reversible. In every month the published values that depend on the option quotes are fewer than
  the quotes they depend on: 9 for 18 in the primary sample, 5 for 12 in the extended sample and 1 for 4
  or 6 in the long sample. A continuous map from more unknowns to fewer values cannot be one-to-one, so
  infinitely many quote vectors reproduce every published value; the sensitivity matrix of the published
  values has full rank in every month, leaving 9, 7 and 3 to 5 directions unconstrained.
- Constructed alternatives. For three months of each sample, spread over its span, I constructed quote
  vectors that reproduce every published value of the month at the file's precision while differing
  from the true quotes materially: by at least 0.5 volatility points in every month, by 2 in two of the
  three primary months and by 4 in every extended and long month. Every alternative calibrates, passes
  the project's arbitrage checks and stays within the range of quotes seen in the data.
- Not traceable on the quote grid. Real quotes lie on a grid (0.0005 volatility points for the composite
  quotes, 0.00625 to 0.025 for Fenics). Within one volatility point of the true quotes, the quote vectors
  on that grid that reproduce every rounded value of a month number at least about 10^3 in the least
  protected month of the long sample, 10^5.9 in the extended sample and 10^12.9 in the primary sample
  (medians 10^8.8, 10^10.1 and 10^18.4). For the extended sample, whose Fenics quotes sit on the coarsest
  grid, this required rounding `extended.csv` to four significant figures instead of six.
- Best attainable estimate. Given the leg identities, the exact spot and forwards, a prior for each quote
  from public information (the realised volatility for at-the-money volatility, zero for the risk
  reversal) and even the exact sensitivities of the published values at the true quotes, the best linear
  estimate leaves a median uncertainty of 1.3 volatility points for at-the-money volatility and 0.7 for
  the 25-delta risk reversal in the primary sample (2.4 and 1.6 in the extended sample; 2.4 for
  at-the-money volatility in the long sample). No risk reversal is pinned down to within 0.10
  volatility points, and about 1% of at-the-money and risk-reversal quotes to within 0.25. The 25-delta
  butterfly varies so little that public information alone predicts it to about 0.09 volatility points;
  the series reduce that only to 0.08.

Before release I withdrew the columns that failed this review: `rr_or`, the long sample's payoff and
premium, and the skew prices and phi of the market reading at 25 delta and of the smile reading. With
those three skew-price series, the same attacker would have pinned about a fifth of primary risk
reversals to within 0.25 volatility points.
