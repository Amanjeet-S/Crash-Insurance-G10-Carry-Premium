# Stage 4: portfolio returns, E2, E3 and E5

Stage 4 of the [research design](../research_design.md). Aggregate results only; month-level series stay in `data/private/results/`.

## Reproduction

```bash
.venv/bin/python scripts/estimate_stage4.py
.venv/bin/python scripts/acquire_cboe_vx.py --retrieval-date 2026-09-24
.venv/bin/python scripts/estimate_e5.py
.venv/bin/python scripts/estimate_moments.py
```

## Returns

- **Legs.** Legs are those of the E1 series: market reading, 3L/3S in the primary sample and 2L/2S in the extended sample. The portfolio sorts on the forward discount, in the manner of Lustig, Roussanov and Verdelhan (2011), who sort currencies into six portfolios on f − s at each month-end and take the highest minus the lowest. With nine currencies the analogue is the top three minus the bottom three.
- **Excess returns.** Arithmetic, per USD of forward notional. They are evaluated at the spot quote on the option expiry date, which is the spot rate for value on the forward's delivery date (`src/qef/fx/crash.py`).
- **Hedged returns.** Each leg adds its protective option's payoff less its premium carried to delivery. Hedges are at 10Δ (primary), 25Δ and ATM.
- **Return windows.** A return window that ends after the last available spot observation has not been realised, and is excluded. This removes the final month-end of the E1 sample, leaving 159 primary month-ends from May 2013 to July 2026, with returns realised from June 2013 to August 2026.

## E2: predictability of the unhedged carry return

The regression is HML^U_{t+1} = a + b x_t, with x_t = φ_t (primary) or the oriented 10Δ risk reversal (quoted; protective-option volatility minus opposite-option volatility, averaged over legs).

- **In sample.** Newey–West errors (automatic bandwidth), and a residual bootstrap under the null of no predictability (9,999 draws).
  - The bootstrap draws the pairs (û_k, v̂_{k+1}) independently. û are residuals of the unrestricted regression; v̂ are AR(1) residuals of the predictor.
  - The bootstrap bias is checked against the first-order approximation of Stambaugh (1999, eq. 18), −(σ_uv/σ_v²)(1 + 3ρ)/T.
- **Out of sample.** Expanding window. The first forecast, made at end-December 2016, is for January 2017; training starts with the extended sample. The test is the MSPE-adjusted statistic of Clark and West (2007) against the historical mean. It is one-sided, with the least-squares standard error recommended for one-step forecasts; Newey–West is reported alongside.

| Predictor | b | HAC t | Bootstrap bias (analytic) | Bias-corrected b | One-sided p (b > 0) | OOS R² | Clark–West t (p) | Clark–West t, NW (p) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| φ | 0.0062 | 5.95 | 0.00055 (0.00049) | 0.0057 | 0.021 | 3.7% | 1.639 (0.051) | 1.50 (0.067) |
| Oriented 10Δ RR | 0.397 | 4.13 | 0.0183 (0.0148) | 0.379 | 0.004 | −5.5% | −1.11 (0.87) | −1.18 (0.88) |

There are 116 out-of-sample forecasts. The correlation between the two predictors is 0.56. The predictors are persistent, with AR(1) coefficients of 0.86 for φ and 0.52 for the risk reversal. Their innovations are negatively correlated with returns (−0.39 and −0.53).

**Interpretation under the pre-registered rules.**

- **In sample.** A higher skew price per unit of carry predicts higher carry returns in the following month, after the small-sample bias correction (one-sided p = 0.021). The sign is the one a crash-compensation reading implies.
- **Out of sample.** The Clark–West statistic of 1.639 falls just short of the 5% one-sided critical value of 1.645, so predictability is not established.
- **The risk reversal alone** predicts in sample but not out of sample.
- **The regime-shift question.** The claim that crash-insurance pricing explains the regime shift requires the E1 difference, a bias-corrected b > 0 and an out-of-sample rejection. Only the second holds, so the evidence is reported as partial support.

## E2, secondary predictors: option-implied variance and skewness

The design treats option-implied variance and skewness as secondary predictors, reported as intervals because quotes stop at the 10Δ strikes (`src/qef/fx/moments.py`; theory results R1, R3, R4).

**Construction.**
- **Moments.** The moments are those of log(X_T/F_X) under the USD one-month forward measure. For USD-base pairs, quote-currency prices are converted with the change of numeraire of R1, E^USD h(S_T) = E^q[h(S_T) S_T/F].
- **Spanning.** The contracts are spanned by out-of-the-money options (R3), adapted from Bakshi, Kapadia and Madan (2003) to contracts centred at the forward.
- **Middle part.** Between the 10Δ strikes the contracts come from the calibrated smile, by composite Simpson's rule with the forward as a node. The step is halved until every contract changes by less than 1e-10 relative, which holds in all 1,440 currency-months.
- **Tails.** Beyond the 10Δ strikes each contract is bounded by R4, with sign-changing weights split into positive and negative parts. Variance and skewness intervals are the exact ranges over the box of contract intervals.
- **Tail exponents, two settings fixed before estimation.**
  - Setting (i): the local elasticities of the smile's density at the 10Δ strikes (median γ = 61, η = 67). This assumes the tails beyond the quotes are no heavier than at the quotes.
  - Setting (ii): γ = η = 2.
- **Portfolio predictors.** Variance, and skewness oriented to each leg (sign reversed for short legs), both averaged over the six E1 legs.

**Identification.**
- **Setting (ii).** No currency-month or portfolio month has a skewness interval that excludes zero, and the variance interval is about fifty times wider than under (i). Truncated quotes therefore do not identify skewness without a tail assumption of the strength of (i).
- **Setting (i).** The median currency-month skewness interval is [−0.70, 0.14], and 32% of intervals exclude zero.
- **A caveat on setting (i).** The calibrated smile's own extrapolation contradicts assumption (i) in about three quarters of currency-months: its elasticity falls below the boundary value within two ATM standard deviations beyond the boundary. Setting (i) intervals are therefore conditional on an assumption the smile itself does not support.

**Predictive regressions, setting (i).** HML^U_{t+1} on the predictor at the lower endpoint, midpoint and upper endpoint of its interval, over 159 months.

| Predictor | b (lower, mid, upper) | Newey–West t | One-sided bootstrap p |
| --- | --- | --- | --- |
| Variance | 4.98, 4.81, 4.64 | 0.87, 0.90, 0.93 | 0.185, 0.171, 0.163 (b > 0) |
| Oriented skewness | −0.024, −0.028, −0.028 | −2.74, −2.60, −1.87 | 0.013, 0.015, 0.045 (b < 0) |

- **Variance.** It does not predict carry returns at any endpoint.
- **Oriented skewness.** More negative oriented skewness, meaning more left-tail risk in the carry position, predicts higher returns: the crash-compensation sign, the same economic direction as the E2 test on φ. The one-sided bootstrap test rejects at 5% at all three endpoints; the Newey–West t does not at the upper endpoint.
- **Reading.** The result holds across the interval only under the bootstrap, and only under setting (i), which the smile's own wings contradict. It is therefore reported as conditional supporting evidence, not as a finding.

## E3: decomposition of the hedge cost

Terms are in basis points of notional per month (primary sample, 159 months, 10Δ hedge), with Newey–West standard errors.

- **Term (i)** is the realised payoff minus the Garman–Kohlhagen value at forecast realised volatility σ̂P. σ̂P is the realised volatility over the previous 21 business days.
- **Term (ii)** is the volatility-level term: minus the flat-ATM premium less that value.
- **Term (iii)** is the skew term: minus the smile premium less the flat-ATM premium.
- The three terms add up to HML^H − HML^U.

| Quantity | Mean | s.e. |
| --- | --- | --- |
| HML^U (unhedged) | 17.7 | 12.1 |
| HML^H, 10Δ | 8.3 | 11.9 |
| HML^H, 25Δ | 9.1 | 11.3 |
| HML^H, ATM | 9.7 | 8.9 |
| (i) payoff less σ̂P value | 0.5 | 3.8 |
| (ii) volatility level | 1.4 | 1.4 |
| (iii) skew | −11.2 | 0.65 |

- **Hedge-cost ratio.** θ_UB = 1 − mean(HML^H)/mean(HML^U) is 0.53 at 10Δ, 0.49 at 25Δ and 0.45 at ATM. The diffusive null of result R6 is θ₀ ≈ 0.10 at 10Δ: the average forward delta of the hedges.
- **Confidence set.** The test-inversion 95% set for θ_UB is unbounded, because the mean unhedged carry return is not significantly different from zero in this sample (t = 1.47). As the design anticipated, the ratio is weakly identified, so no point value is claimed.
- **What the decomposition shows.** The hedged-minus-unhedged difference over 2013–2026 consists almost entirely of the skew term, which equals the ex-ante skew price of E1 carried to delivery. Term (iii) is the ex-ante skew premium, known at each month-end, so its small standard error reflects how little that premium varies over time, not the precision of realised crash compensation. The realised terms, the payoff term and the volatility-level term, are small and not distinguishable from zero.

**Extended sample.** Fenics quotes, 2007–2013, 72 months. The unhedged mean is 8.4 bp (s.e. 62.7). Term (i) is 16.7 bp (s.e. 10.2), reflecting option payoffs in 2008, and the skew term is −22.3 bp (s.e. 3.0).

**Supplementary (post hoc): realised returns by regime.** The design measures the change in carry through the forward-discount spread (E1). The realised means, with Newey–West standard errors and the stationary-bootstrap 95% interval for the difference (9,999 draws), are:

| Portfolio | Zero-rate (104 months) | Hiking (55 months) | Difference (s.e.) | Bootstrap 95% |
| --- | --- | --- | --- | --- |
| HML^U | 7.4 (12.4) | 37.0 (17.8) | 29.6 (23.8) | [−19.5, 80.0] |
| HML^H, 10Δ | −2.7 (15.7) | 29.3 (16.0) | 32.0 (22.4) | [−16.8, 80.3] |

Realised carry returns are higher in the hiking regime, but the difference is not significant; monthly returns are too noisy over 55 months to measure the change in the realised premium.

## E5: systemic-risk exposure

Both factors are measured over each carry return window, from the month-end to the one-month option expiry.

- **VIX roll-down.** Minus the percentage change of the second-month VIX future's settlement price, holding the same contract.
  - This is a month-end analogue of the roll-down strategy of Caballero and Doyle (2012), who short the VIX future whose expiry matches the one-month forward's maturity and hold it to expiry.
  - Settlements are Cboe daily files, one per contract. The pre-2007 archive quotes contracts at ten times the index level; the change of scale is detected in the files (26 March 2007) and removed.
  - At every month-end, Cboe second-month settlements equal the LSEG second-month continuation (223 of 223).
- **Global FX-volatility innovation.** The FX-volatility level follows Menkhoff, Sarno, Schmeling and Schrimpf (2012, eq. 4): the average over the window's days of the cross-sectional mean absolute daily log spot change of the nine currencies. The innovation is the residual of an AR(1) fitted to that series (coefficient 0.76).

Regressions are on both factors, with Newey–West t-statistics in parentheses. The primary sample has 159 months. α is in basis points per month.

| Portfolio | α | β, VIX roll-down | β, FX-volatility innovation | R² |
| --- | --- | --- | --- | --- |
| HML^U | −7.3 (−0.61) | 0.049 (7.41) | −3.11 (−2.16) | 0.33 |
| HML^H, 10Δ | −10.6 (−0.79) | 0.039 (4.69) | −1.76 (−1.21) | 0.22 |

- **Exposure.** Unhedged carry returns load strongly on the VIX roll-down (correlation 0.55), consistent with Caballero and Doyle (2012). The 10Δ hedge lowers the loading by about a fifth but leaves most of it.
- **Adjusted hedge-cost ratio.** After adjusting for both factors, neither portfolio earns a significant α. The exposure-adjusted ratio 1 − α_H/α_U is −0.45, and its test-inversion set is unbounded.
- **Extended sample (65 months).** The loadings are larger: 0.124 unhedged (t = 5.39) and 0.112 hedged (t = 4.53). The αs are not significant.
