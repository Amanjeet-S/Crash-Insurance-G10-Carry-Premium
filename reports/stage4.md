# Stage 4: portfolio returns, E2, E3 and E5

This record covers Stage 4 of the [research design](../research_design.md). It reports aggregate results only; month-level series stay in `data/private/results/`.

## Reproduction

```bash
.venv/bin/python scripts/estimate_stage4.py
.venv/bin/python scripts/acquire_cboe_vx.py --retrieval-date 2026-09-24
.venv/bin/python scripts/estimate_e5.py
.venv/bin/python scripts/estimate_moments.py
# Computed after the main results (sections below marked accordingly)
.venv/bin/python scripts/estimate_long_atm.py
.venv/bin/python scripts/estimate_design_checks.py
```

## Returns

The legs are those of the E1 series under the market reading, three long and three short in the primary sample and two and two in the extended sample. The portfolio sorts on the forward discount in the manner of Lustig, Roussanov and Verdelhan (2011), who sort currencies into six portfolios on f − s at each month-end and take the highest minus the lowest; with nine currencies the analogue is the top three minus the bottom three.

Excess returns are arithmetic, per USD of forward notional, and are evaluated at the spot quote on the option expiry date, which is the spot rate for value on the forward's delivery date (`src/qef/fx/crash.py`). A hedged return adds each leg's protective-option payoff less its premium carried to delivery, with hedges at 10Δ (primary), 25Δ and ATM. A return window that ends after the last available spot observation has not been realised and is excluded. This removes the final month-end of the E1 sample, leaving 159 primary month-ends from May 2013 to July 2026, with returns realised from June 2013 to August 2026.

## E2: predictability of the unhedged carry return

The regression is HML^U_{t+1} = a + b x_t, with x_t = φ_t (primary) or the oriented 10Δ risk reversal, the quoted volatility of the protective option minus that of the opposite option, averaged over legs.

In sample the slope has Newey–West errors with the automatic bandwidth and is tested by a residual bootstrap under the null of no predictability with 9,999 draws. The bootstrap draws the pairs (û_k, v̂_{k+1}) independently, where û are residuals of the unrestricted regression and v̂ are AR(1) residuals of the predictor, and its bias is checked against the first-order approximation of Stambaugh (1999, eq. 18), −(σ_uv/σ_v²)(1 + 3ρ)/T. Out of sample the forecasts use an expanding window whose training starts with the extended sample; the first forecast, made at end-December 2016, is for January 2017. The test is the one-sided MSPE-adjusted statistic of Clark and West (2007) against the historical mean, with the least-squares standard error they recommend for one-step forecasts, and the Newey–West version is reported alongside.

| Predictor | b | HAC t | Bootstrap bias (analytic) | Bias-corrected b | One-sided p (b > 0) | OOS R² | Clark–West t (p) | Clark–West t, NW (p) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| φ | 0.0062 | 5.95 | 0.00055 (0.00049) | 0.0057 | 0.021 | 3.7% | 1.639 (0.051) | 1.50 (0.067) |
| Oriented 10Δ RR | 0.397 | 4.13 | 0.0183 (0.0148) | 0.379 | 0.004 | −5.5% | −1.11 (0.87) | −1.18 (0.88) |

There are 116 out-of-sample forecasts, and the correlation between the two predictors is 0.56. Both predictors are persistent, with AR(1) coefficients of 0.86 for φ and 0.52 for the risk reversal, and their innovations are negatively correlated with returns (−0.39 and −0.53).

In sample, a higher skew price per unit of carry predicts higher carry returns in the following month after the small-sample bias correction (one-sided p = 0.021), which is the sign a crash-compensation reading implies. Out of sample the Clark–West statistic of 1.639 does not reach the 5% one-sided critical value of 1.645, so predictability is not established. The risk reversal on its own predicts in sample but not out of sample. The claim that crash-insurance pricing explains the regime shift requires the E1 difference, a bias-corrected b > 0 and an out-of-sample rejection; only the second holds, so the evidence is reported as partial support.

## Secondary predictors: option-implied variance and skewness

The design treats option-implied variance and skewness as secondary predictors, reported as intervals because the quotes stop at the 10Δ strikes (`src/qef/fx/moments.py`; theory results R1, R3 and R4). The moments are those of log(X_T/F_X) under the USD one-month forward measure. For USD-base pairs, quote-currency prices are converted with the change of numeraire of R1, E^USD h(S_T) = E^q[h(S_T) S_T/F]. The contracts are spanned by out-of-the-money options (R3), adapted from Bakshi, Kapadia and Madan (2003) to contracts centred at the forward.

Between the 10Δ strikes the contracts come from the calibrated smile by composite Simpson's rule with the forward as a node, and the step is halved until every contract changes by less than 1e-10 relative, which holds in all 1,440 currency-months. The design's closed-form checks of this code were computed after the main results (`scripts/estimate_design_checks.py`): in all 1,440 currency-months the adaptive tail quadrature agrees with the closed-form R4 bounds to 1.6e-15 relative in setting (i) and 6.7e-16 in setting (ii), no tail integral stopped at the subdivision limit before its relative tolerance of 1e-12, and the Simpson rule needed a median of 512 and at most 4,096 subintervals per side, with a largest relative change at the last halving of 9.97e-11, below the tolerance of 1e-10. Beyond the 10Δ strikes each contract is bounded by R4, with sign-changing weights split into positive and negative parts, and the variance and skewness intervals are the exact ranges over the box of contract intervals. Two settings for the tail exponents were fixed before estimation. Setting (i) uses the local elasticities of the smile's density at the 10Δ strikes (median γ = 61, η = 67), which assumes that the tails beyond the quotes are no heavier than at the quotes; setting (ii) uses γ = η = 2. The portfolio predictors are the variance and the skewness oriented to each leg, with the sign reversed for short legs, both averaged over the six E1 legs.

Under setting (ii) no currency-month or portfolio month has a skewness interval that excludes zero, and the variance interval is about fifty times wider than under (i), so truncated quotes do not identify skewness without a tail assumption of the strength of (i). Under setting (i) the median currency-month skewness interval is [−0.70, 0.14], and 32% of intervals exclude zero. The calibrated smile's own extrapolation, however, contradicts assumption (i) in about three quarters of currency-months, where its elasticity falls below the boundary value within two ATM standard deviations beyond the boundary. Setting (i) intervals are therefore conditional on an assumption that the smile itself does not support.

The predictive regressions use setting (i), with HML^U_{t+1} regressed on the predictor at the lower endpoint, midpoint and upper endpoint of its interval over 159 months.

| Predictor | b (lower, mid, upper) | Newey–West t | One-sided bootstrap p |
| --- | --- | --- | --- |
| Variance | 4.98, 4.81, 4.64 | 0.87, 0.90, 0.93 | 0.185, 0.171, 0.163 (b > 0) |
| Oriented skewness | −0.024, −0.028, −0.028 | −2.74, −2.59, −1.87 | 0.013, 0.015, 0.045 (b < 0) |

Variance does not predict carry returns at any endpoint. More negative oriented skewness, meaning more left-tail risk in the carry position, predicts higher returns, which is the crash-compensation sign and the same economic direction as the E2 test on φ. The one-sided bootstrap test rejects at 5% at all three endpoints, but the Newey–West t does not at the upper endpoint. Because the result holds across the interval only under the bootstrap, and only under a setting that the smile's own wings contradict, I report it as conditional supporting evidence rather than as a finding.

## Post hoc: how strong a tail assumption the sign of skewness needs

This section was added on 27 September 2026, after the results above were recorded, and is post hoc. It applies theory result R10 (`theory/notes.tex`), which gives, under the hypotheses of R4, the exponents compatible with the smile's boundary price and slope and the sharp identified set of the raw moments.

```bash
.venv/bin/python scripts/estimate_identification.py
```

For each of the 1,440 currency-months the largest compatible exponents are γ̄ = K_min G0/P0 − 1 and η̄ = 1 + K_max Ḡ0/C0, with quartiles 47.7, 59.7 and 75.4 for γ̄ and 54.5, 67.8 and 84.0 for η̄. Setting (i) requires its exponents to be at most these. It fails the requirement in the lower tail in 54.4% of currency-months, in the upper tail in 49.3% and in one or both in 75.3%; the ratio γ_i/γ̄ has quartiles 0.97, 1.01 and 1.05, so setting (i) sits at the strongest assumption the boundary prices admit. Every currency-month in which the check refutes setting (i) also fails the diagnostic above, which extrapolates the smile two ATM standard deviations beyond the boundary, while the diagnostic fails in a further 309 lower and 357 upper tails that pass the check; the check uses only the boundary price, slope and curvature, and it is necessary, not sufficient.

The identified sets are computed exactly for a grid of step points (1,000 per tail); refining the grid fourfold in 40 currency-months moved no bound on the third central moment by more than 1.2×10⁻³ of its interval width and left every breakdown value unchanged. The table gives the share of currency-months in which the sign of the skewness is identified, with θ the fraction of the maximal exponents (γ = θγ̄, η = 1 + θ(η̄ − 1)).

| Tail exponents | Admissible | Sign identified (sharp set) | Of which negative | Sign identified (R4 box) |
| --- | --- | --- | --- | --- |
| Setting (ii), γ = η = 2 | 100% | 0.0% | 0.0% | 0.0% |
| Setting (i), boundary elasticities | 24.7% | 93.3% | 86.8% | 32.1% of all currency-months |
| θ = 0.25 | 100% | 0.7% | 0.6% | |
| θ = 0.5 | 100% | 26.7% | 21.9% | |
| θ = 0.75 | 100% | 70.3% | 57.5% | |
| θ = 0.9 | 100% | 87.6% | 70.3% | |

The shares of negative signs are shares of all admissible currency-months. The breakdown value θ*, the smallest fraction at which the sign is identified, has quartiles 0.49, 0.63 and 0.80; it is at most 0.5 in 26.2% of currency-months and at most 0.75 in 69.9%. At θ = 1, where both tails are pure power laws and the moments are point-identified, the skewness has quartiles −0.535, −0.337 and −0.082 and is negative in 78.5% of currency-months. By currency, the median θ* is lowest for AUD (0.54) and NZD (0.58) and highest for CHF (0.76) and EUR (0.71).

For comparison, in Merton and Heston models with one-month FX-typical parameters the models' own tail exponents lie between 0.58 and 0.97 of the maximal ones, and the sign of the skewness is identified at those exponents whenever the skewness is not close to zero ([model validation](model_validation.md)). Truncated one-month quotes therefore identify the sign of risk-neutral skewness only under tail assumptions close to the strongest the boundary prices admit: in the median currency-month the tail exponents beyond the 10Δ strikes must be at least 63% of the largest exponents that the boundary price and slope admit. Standard models meet that requirement in some cases and not in others, and setting (i), on which the secondary E2 regressions above rest, is excluded by the boundary prices themselves in three quarters of currency-months.

## Post hoc: an out-of-sample test of the tail hypotheses with 5Δ quotes

The hypotheses behind the results above concern the tails beyond the 10Δ strikes, which the quotes used so far do not reach. The Fenics contributor quotes one-month 5Δ risk reversals and butterflies for all nine currencies from June 2022, and these quotes had entered no calibration. I recorded the test in the research log and committed it before retrieving their history. The boundary is the Fenics market-reading smile at its 10Δ strikes; the 5Δ put and call volatilities are σ_ATM + BF5 ∓ RR5/2 (the smile-strangle reading), and the tests use the closed-form identified interval of a price beyond the boundary that follows from R10.

```bash
.venv-lseg/bin/python scripts/acquire_lseg_fx.py --retrieval-date 2026-09-27 --skip-search --only <the eighteen RICs <CCY>1MR5=FN and <CCY>1MB5=FN>
.venv/bin/python scripts/estimate_wing_test.py
```

All 450 currency-months from July 2022 to August 2026 have both 5Δ quotes and a converged Fenics smile. The 5Δ risk reversal has the sign of the Fenics 10Δ risk reversal in 99.8% of them, and every 5Δ strike lies beyond the smile's 10Δ strike. The 5Δ prices lie below the smile's own extrapolation, with quote-to-smile ratios of 0.94, 0.96 and 1.00 (quartiles) for puts and 0.93, 0.96 and 0.98 for calls, so the SABR smile extrapolated from 25Δ quotes overprices 5Δ protection by about 4% in the median currency-month.

| Out-of-sample result (450 currency-months) | Lower tail | Upper tail |
| --- | --- | --- |
| 5Δ price consistent with the weakest hypothesis of the class | 99.1% | 99.6% |
| θ_5, largest consistent fraction of the maximal exponent (quartiles) | 0.82, 0.88, 0.95 | 0.81, 0.86, 0.93 |
| θ_5 = 1 (tail exactly the power law matching the 10Δ boundary) | 0.0% | 0.0% |

The class is consistent with the 5Δ prices in both tails in 99.1% of currency-months; the four exceptions are all USDJPY months. With both 5Δ prices as constraints and the exponents at the largest values they allow, the sign of the skewness is identified in 90.6% of the 446 consistent currency-months, negative in 71.5% of them. Along the common path θ, with the 5Δ constraints, the sign is identified in 394 currency-months, and the breakdown point there has quartiles 0.36, 0.45 and 0.57. As a supplementary comparison not in the plan, the breakdown point without the 5Δ prices on the same smiles and months has quartiles 0.54, 0.64 and 0.79, close to the composite smiles' 0.54, 0.65 and 0.78 over the same months; the 5Δ prices lower it in every currency-month where the sign is identified with them, by a median of 0.16. The boundary check refutes setting (i) on these smiles in 57.3% of currency-months.

By the rule fixed in the plan, in 90.6% of currency-months the 5Δ quotes do not rule identification out: the sign of the skewness is identified under the strongest hypothesis of the class they allow. In the remaining 9.4% it is identified under no hypothesis of the class consistent with them. Consistency with one extra price is necessary for a hypothesis, not sufficient, so the favourable share is an upper limit. The tail strength the sign needs, θ* of about 0.64 in the median currency-month without the 5Δ prices and 0.45 with them, lies below the strength the 5Δ prices allow, θ_5 of about 0.87. By currency, USDJPY has the heaviest tails relative to the power-law benchmark (median θ_5 of 0.62 in the lower and 0.72 in the upper tail) and all four inconsistent months. The sample is short, lies entirely in the hiking regime, and the 5Δ quotes are dealer indications from one contributor.

## E3: decomposition of the hedge cost

The terms below are in basis points of notional per month, for the primary sample of 159 months and the 10Δ hedge, with Newey–West standard errors. Term (i) is the realised payoff minus the Garman–Kohlhagen value at the forecast realised volatility σ̂P, the realised volatility over the previous 21 business days. Term (ii), the volatility-level term, is minus the flat-ATM premium less that value, and term (iii), the skew term, is minus the smile premium less the flat-ATM premium. The three terms add up to HML^H − HML^U.

| Quantity | Mean | s.e. |
| --- | --- | --- |
| HML^U (unhedged) | 17.7 | 12.1 |
| HML^H, 10Δ | 8.3 | 11.9 |
| HML^H, 25Δ | 9.1 | 11.3 |
| HML^H, ATM | 9.7 | 8.9 |
| (i) payoff less σ̂P value | 0.5 | 3.8 |
| (ii) volatility level | 1.4 | 1.4 |
| (iii) skew | −11.2 | 0.65 |

The hedge-cost ratio θ_UB = 1 − mean(HML^H)/mean(HML^U) is 0.53 at 10Δ, 0.49 at 25Δ and 0.45 at ATM, against the diffusive null of result R6, θ₀ ≈ 0.10 at 10Δ, the average forward delta of the hedges. The test-inversion 95% set for θ_UB, Fieller's (1954) construction with a Newey–West variance, is unbounded, because the mean unhedged carry return is not significantly different from zero in this sample (t = 1.47). As the design anticipated, the ratio is weakly identified and no point value is claimed.

The total 10Δ hedge cost, HML^H − HML^U, is −9.3 bp per month with a Newey–West standard error of 4.0 (t = −2.35). This test, and the shape of the unbounded sets described next, were added after a first run of `scripts/estimate_design_checks.py` and are post hoc. HML^H − HML^U is the net profit of the protective options, payoff less premium carried to delivery, and it is strongly right-skewed (skewness 3.6), which gives the t-statistic a heavier left tail than the normal. A studentised stationary bootstrap (bootstrap-t; Hall, 1992, Chapter 3), with the expected block length of Politis and White (2004) as corrected by Patton, Politis and White (2009), 9,999 draws, seed 20260924, and the Newey–West standard error recomputed at the automatic lag in each draw, puts the 2.5% quantile of t at −2.78. Its 95% interval for the mean is [−16.4, 1.7] bp, which contains zero, so a two-sided test at 5% does not reject (the one-sided bootstrap p-value is 0.04), and I do not claim that the realised hedge cost differs from zero. For the 25Δ and ATM hedges the means are −8.6 and −8.0 bp, with t = −1.07 and −0.89 and bootstrap-t intervals [−24.0, 9.6] and [−26.8, 10.4] bp. The 10Δ result depends on how many large payoffs the sample contains: without the month with the largest net profit, t = −3.55, while over the extended sample, which contains 2008 but holds two legs per side priced from Fenics smiles, the same t is −0.25.

At 10Δ the test-inversion set is the line less the interval from −1.231 to 0.093. It retains θ₀ = 0.100 only narrowly: t = −1.93 at the automatic Newey–West lag of 6, and the test rejects θ₀ at every fixed lag from 0 to 5. The 25Δ and ATM sets contain the whole scanned grid [−10, 10]. The delta-method intervals that the design names as secondary are [−0.21, 1.27] at 10Δ, [−0.30, 1.27] at 25Δ and [−0.30, 1.20] at ATM (delta-method standard errors 0.38, 0.40 and 0.38). They are bounded because the delta method holds the standard error of the mean of H − (1 − v)U at its value at θ_UB, while the test-inversion set uses its value at each v; as |v| grows that standard error grows like |v| times the standard error of mean(HML^U), so every v far from θ_UB is retained when mean(HML^U) is insignificant. Procedures that always return bounded sets have zero worst-case coverage for such ratios (Gleser and Hwang, 1987; Dufour, 1997).

The hedged-minus-unhedged difference over 2013–2026 consists almost entirely of the skew term, which equals the ex-ante skew price of E1 carried to delivery. Term (iii) is the ex-ante skew premium, known at each month-end, so its small standard error reflects how little that premium varies over time, not the precision of realised crash compensation. The realised terms, the payoff term and the volatility-level term, are small and not distinguishable from zero.

In the extended sample of Fenics quotes, 72 months over 2007–2013, the unhedged mean is 8.4 bp (s.e. 62.7). Term (i) is 16.7 bp (s.e. 10.2), reflecting option payoffs in 2008, and the skew term is −22.3 bp (s.e. 3.0).

## The long at-the-money sample (computed after the main results)

The design defines a long ATM sample: ATM-hedged carry as in Burnside et al. (2011), from the start of the ATM and forward series. `scripts/estimate_long_atm.py` states its rules. It uses composite spot, one-month forward points and one-month ATM volatility only, with no Fenics splice. A currency is available at a month-end if all three quotes are observed or substituted under the five-business-day rule; the portfolio holds three long and three short legs when at least six currencies are available, two and two when four or five are, and the month is excluded otherwise. The hedge is at the delta-neutral-straddle strike in each pair's convention, which does not depend on the base-currency discount factor, and the premium is carried to delivery as in Stage 4, so no interest rate enters. The first month-end is January 1995, the first of the provider's one-month ATM series: a coverage retrieval on 1 October 2026 for history from 1970 returned no composite ATM quote before 6 January 1995, although spot quotes go back to 1971 and forward points to 1982; SEK enters in 1998, EUR in 1999, NZD in 1999 and NOK in 2004. There are 379 month-ends with realised returns, to July 2026: 36 with two legs per side (to December 1997) and 343 with three. Twenty ATM quotes and 12 leg quotes were substituted, and 2 currency-months were excluded after a currency's first available month-end.

I fixed these rules after the Stage 4 results had been recorded, so they were not fixed blind: the 159 primary-window months reproduce the Stage 4 HML^U and ATM-hedged series exactly (largest difference 2e-12 bp, the same legs in every month), and the Fenics extended-sample ATM series was known. No return before May 2013 had been computed under these rules. A draft of the script had been run once, and I overwrote its outputs without reading them; relative to that draft, two sub-periods and a staleness diagnostic were added.

Means in bp per month, Newey–West standard errors in parentheses. The 95% set for θ_UB is by test inversion; it is unbounded on both sides when |t| of mean(HML^U) is below 1.96, because the automatic bandwidth does not depend on the scale of the series.

| Window (formation month-ends) | Months | HML^U | HML^H, ATM | H − U | θ_UB | 95% set |
| --- | --- | --- | --- | --- | --- | --- |
| Full, 1995-01 to 2026-07 | 379 | 37.2 (13.2) | 16.3 (7.7) | −20.9 (7.9) | 0.56 | [0.27, 0.94] |
| 1995-01 to 2007-12 | 156 | 61.4 (19.4) | 27.4 (13.4) | −34.0 (10.6) | 0.55 | [0.33, 0.94] |
| 2008 | 12 | −222.3 (156.4) | −98.5 (39.1) | 123.8 (136.9) | 0.56 | unbounded |
| 2009-01 to 2026-07 | 211 | 34.1 (13.5) | 14.6 (8.0) | −19.5 (11.6) | 0.57 | [−0.23, 1.09] |
| 2009-01 to 2021-12 | 156 | 33.1 (16.9) | 10.3 (10.2) | −22.8 (12.9) | 0.69 | unbounded |
| 2022-01 to 2026-07 | 55 | 37.0 (17.8) | 26.8 (11.7) | −10.2 (20.8) | 0.28 | [−6.02, 0.73] ∪ [1.25, 1.38] |
| Full without 2008 | 367 | 45.7 (11.5) | 20.0 (7.8) | −25.7 (7.5) | 0.56 | [0.33, 0.85] |
| Before the primary start | 220 | 51.4 (20.7) | 21.0 (11.6) | −30.3 (14.3) | 0.59 | [0.224, 0.239] ∪ [0.243, 1.066] |
| Primary window, 2013-05 to 2026-07 | 159 | 17.7 (12.1) | 9.7 (8.9) | −8.0 (8.9) | 0.45 | unbounded |

The diffusive null for an ATM hedge, θ₀, the average absolute forward delta of the hedging options, is between 0.494 and 0.498 in every window, and every set contains it. For dollar-base legs the hedge ratio in dollar returns is the premium-adjusted forward delta, so the pips form of the design is a leading-order approximation; with the premium-adjusted form θ₀ rises by 0.002 to 0.006 (to 0.4999 over the full sample). Over 1995–2026 the unhedged premium is significant (t = 2.83) and the set is bounded, but it contains θ₀, so the ATM hedge gives up about what a diffusion predicts; an ATM hedge has no skew term and the set is wide, so this sample is uninformative about priced crash risk. In 2008 mean(HML^U) is negative and the hedge raised the mean return: the option payoff averaged 502 bp against a premium of 378 bp, so the positive θ_UB of that year is not a cost share. Over the extended sample the Fenics ATM-hedged series (two and two) correlates at 0.94 with this one.

At the primary sample's ratio of the mean unhedged return to its long-run standard deviation, the unhedged t-statistic reaches 1.96, and so the θ_UB set becomes bounded, after about 159 × (1.96/1.47)² ≈ 280 months. The one-month risk-reversal quotes that an out-of-the-money hedge needs begin in January 2007 from Fenics and between July 2010 and May 2013 from the composite in the licensed panel, against January 1995 for the ATM volatility.

The design measures the change in carry through the forward-discount spread of E1. A supplementary comparison of realised returns by regime, added after the results above and therefore post hoc, gives the following means, with Newey–West standard errors and the stationary-bootstrap 95% interval for the difference (9,999 draws).

| Portfolio | Zero-rate (104 months) | Hiking (55 months) | Difference (s.e.) | Bootstrap 95% |
| --- | --- | --- | --- | --- |
| HML^U | 7.4 (12.4) | 37.0 (17.8) | 29.6 (23.8) | [−19.5, 80.0] |
| HML^H, 10Δ | −2.7 (15.7) | 29.3 (16.0) | 32.0 (22.4) | [−16.8, 80.3] |

Realised carry returns are higher in the hiking regime, but the difference is not significant; monthly returns are too noisy over 55 months to measure the change in the realised premium.

## E5: systemic-risk exposure

Both factors are measured over each carry return window, from the month-end to the one-month option expiry. The VIX roll-down factor is minus the percentage change of the second-month VIX future's settlement price, holding the same contract. It is a month-end analogue of the roll-down strategy of Caballero and Doyle (2012), who short the VIX future whose expiry matches the one-month forward's maturity and hold it to expiry. Settlements come from Cboe daily files, one per contract. The pre-2007 archive quotes contracts at ten times the index level; the change of scale is detected in the files (26 March 2007) and removed, after which the Cboe second-month settlements equal the LSEG second-month continuation at every month-end (223 of 223).

The global FX-volatility level follows Menkhoff, Sarno, Schmeling and Schrimpf (2012, eq. 4): the average over the window's days of the cross-sectional mean absolute daily log spot change of the nine currencies. The innovation is the residual of an AR(1) fitted to that series (coefficient 0.76).

The regressions are on both factors, with Newey–West t-statistics in parentheses. The primary sample has 159 months, and α is in basis points per month.

| Portfolio | α | β, VIX roll-down | β, FX-volatility innovation | R² |
| --- | --- | --- | --- | --- |
| HML^U | −7.3 (−0.61) | 0.049 (7.41) | −3.11 (−2.16) | 0.33 |
| HML^H, 10Δ | −10.6 (−0.79) | 0.039 (4.69) | −1.76 (−1.21) | 0.22 |

Unhedged carry returns load strongly on the VIX roll-down (correlation 0.55), consistent with Caballero and Doyle (2012), and the 10Δ hedge lowers the loading by about a fifth but leaves most of it. After adjusting for both factors neither portfolio earns a significant α; the exposure-adjusted ratio 1 − α_H/α_U is −0.45, and its test-inversion set is unbounded. In the extended sample of 65 months the loadings are larger, 0.124 unhedged (t = 5.39) and 0.112 hedged (t = 4.53), and the αs are not significant.
