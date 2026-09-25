# Robustness grid

Section 8 of the [research design](../research_design.md). Aggregate results only; month-level series stay in `data/private/results/`.

## Reproduction

```bash
.venv/bin/python scripts/calibrate_vv.py
.venv/bin/python scripts/robustness.py
.venv/bin/python scripts/calibrate_smiles.py --tenor 3M --extra-raw-root data/private/lseg/2026-09-24/raw
.venv/bin/python scripts/estimate_3m.py
```

The base variant reproduces the E1 regime difference and the Stage 4 means of HML^U, HML^H and the skew term to machine precision, and the script stops if it does not. Bootstrap intervals and the E2 bootstrap use 1,999 draws, so they differ from the 9,999-draw base figures in [e1.md](e1.md) and [stage4.md](stage4.md) in the second decimal.

## Definitions

The design lists the variants. The rules below were fixed in `scripts/robustness.py` before any variant was estimated.

- **Stale butterflies.** The audit rule: a currency is dropped at a month-end if its selected 25Δ butterfly quote lies in a run of five or more identical two-sided business-day quotes. 109 currency-months are flagged.
- **Exclusions.** A month-end is removed if its quote date or return window falls in the excluded period. For March 2020, this removes the February and March 2020 month-ends. For CHF, CHF is removed from the ranking at the November 2014 to March 2015 month-ends.
- **Implementable returns.** Forwards are traded at the quoted bid or ask and closed at the spot bid or ask on the expiry date. Options are bought at the smile volatility plus k times the ATM half-spread, converted to premium by vega (design, section 5). The skew term keeps its mid-price definition.
- **Dollar carry.** Lustig, Roussanov and Verdelhan (2014, Section 3.2): long all available currencies with equal weights when their average forward discount is positive, short all otherwise. Their signal is the developed-country average; here it is the nine-currency average.
- **Ten currencies.** USD enters the ranking with fd = 0. A USD leg has zero return and no option.
- **Log returns.** Each leg's forward return is replaced by ±(s_{t+1} − f_t). Option payoffs and premia are unchanged.
- **Previous-day spot.** Option payoffs use the spot of the New York business day before expiry. Forward returns keep the expiry-date spot.
- **Vanna–volga smiles.** The vanna–volga price smile of Castagna and Mercurio (2007, eqs. (6)–(7)), with Greeks at the ATM volatility, pillars at the ATM delta-neutral-straddle strike and the 25Δ strikes at the pillar volatilities in the pair's convention, and the smile strangle solved so that the smile prices the market strangle (Reiswich, 2010, Section 3.3.3; Bossens et al., 2010, Section 3.3). The first- and second-order approximations are checks only, because they are expansions of the price and the first-order form overvalues the wings (Castagna and Mercurio, 2007, p. 10). Protective options are struck at the 10Δ strikes of the vanna–volga smile.
- **Out of sample.** Every variant, the base included, uses an expanding window that starts at the primary start, with the first forecast for January 2017 (116 forecasts). Stage 4 trains on the extended sample from 2007, so its base statistic (1.639) is not comparable with the column below.

The three-month tenor is reported in its own section below.

## Results

Primary sample. φ is the E1 ratio (market reading, 10Δ unless stated), and Δφ is its hiking-minus-zero-rate difference, with the Newey–West standard error in parentheses and the stationary-bootstrap 95% interval. b is the bias-corrected E2 slope with its one-sided bootstrap p-value. CW is the Clark–West t-statistic. Returns and the skew term are in basis points per month, with Newey–West t-statistics in parentheses.

| Variant | Months (φ, returns) | Δφ (s.e.) | Bootstrap 95% | b (p) | CW | HML^U | HML^H | Skew term | θ_UB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Base | 160, 159 | −0.355 (0.204) | [−0.82, 0.01] | 0.0058 (0.017) | 0.96 | 17.7 (1.47) | 8.3 (0.70) | −11.2 | 0.53 |
| 25Δ hedge | 160, 159 | −0.272 (0.166) | [−0.65, 0.02] | 0.0074 (0.013) | 0.75 | 17.7 (1.47) | 9.1 (0.81) | −8.4 | 0.49 |
| ATM hedge | 159 | not defined | | | | 17.7 (1.47) | 9.7 (1.09) | 0 | 0.45 |
| Smile-strangle reading | 160, 159 | −0.345 (0.196) | [−0.80, 0.01] | 0.0059 (0.018) | 0.99 | 17.7 (1.47) | 8.7 (0.73) | −10.8 | 0.51 |
| Fenics quotes | 160, 159 | −0.347 (0.203) | [−0.81, 0.02] | 0.0058 (0.019) | 0.98 | 17.7 (1.47) | 8.5 (0.71) | −11.2 | 0.52 |
| Stale butterflies missing | 156, 155 | −0.440 (0.265) | [−1.04, 0.02] | 0.0037 (0.023) | 0.19 | 17.7 (1.53) | 8.8 (0.76) | −11.2 | 0.50 |
| Excluding March 2020 | 158, 157 | −0.329 (0.193) | [−0.76, 0.02] | 0.0053 (0.035) | 1.31 | 19.9 (1.73) | 8.8 (0.77) | −11.0 | 0.56 |
| Excluding CHF, Dec 2014 to Mar 2015 | 160, 159 | −0.353 (0.204) | [−0.82, 0.01] | 0.0056 (0.016) | 0.84 | 20.4 (1.87) | 10.4 (0.89) | −11.2 | 0.49 |
| Implementable, k = 1 | 160, 159 | as base | as base | 0.0059 (0.016) | 0.92 | 5.4 (0.45) | −6.7 (−0.55) | −11.2 | not meaningful |
| Implementable, k = 1.5 | 160, 159 | as base | as base | 0.0059 (0.016) | 0.92 | 5.4 (0.45) | −8.0 (−0.66) | −11.2 | not meaningful |
| Implementable, k = 2 | 160, 159 | as base | as base | 0.0059 (0.016) | 0.92 | 5.4 (0.45) | −9.4 (−0.78) | −11.2 | not meaningful |
| Dollar carry | 160, 159 | −1.060 (0.559) | [−2.12, −0.37] | 0.0002 (0.373) | −1.51 | 0.3 (0.02) | −2.6 (−0.17) | −2.8 | not meaningful |
| Ten currencies | 160, 159 | −0.300 (0.149) | [−0.63, −0.02] | 0.0053 (0.052) | 1.17 | 15.3 (1.37) | 8.2 (0.74) | −9.7 | 0.47 |
| Log returns | 160, 159 | as base | as base | 0.0056 (0.020) | 0.86 | 16.7 (1.38) | 7.4 (0.62) | −11.2 | 0.56 |
| Previous-day spot for payoffs | 160, 159 | as base | as base | as base | as base | 17.7 (1.47) | 5.2 (0.43) | −11.2 | 0.70 |
| Vanna–volga smiles | 160, 159 | −0.337 (0.190) | [−0.77, 0.00] | 0.0060 (0.021) | 1.08 | 17.7 (1.47) | 8.7 (0.73) | −10.8 | 0.51 |

θ_UB is marked not meaningful where the mean unhedged return is within half a standard error of zero.

## Interpretation under the pre-registered rules

- **E1.** φ is lower in the hiking regime in every variant. The difference is significant at 5% under both inferences only in the ten-currency ranking (t = −2.01, interval [−0.63, −0.02]). The dollar-carry interval excludes zero, but that φ divides by the absolute average forward discount, which is below 5 bp in 40 zero-rate months, so its difference comes mainly from the denominator. In the base and ten-currency portfolios the skew price per month is nearly the same in both regimes (11.2 and 11.2 bp; 10.0 and 9.3 bp), and the fall in φ comes from the rise in the forward-discount spread, as in the post hoc split of [e1.md](e1.md). The pre-registered E1 conclusion, no significant regime difference, therefore depends on the portfolio definition, and the evidence points to a shift in carry rather than in the price of crash insurance.
- **E2.** The bias-corrected slope is positive with one-sided p ≤ 0.035 in every HML variant, and p = 0.052 for ten currencies. It is zero for dollar carry, where φ is dominated by the denominator. No variant rejects out of sample at 5% when training starts in 2013. The Stage 4 statistic of 1.639 used the longer 2007 training window, so the out-of-sample evidence is sensitive to the training start.
- **E3.** The skew term is −8.4 to −11.2 bp per month across the HML variants. It is the ex-ante premium, so its small standard error reflects its stability over time, not the precision of realised crash compensation. θ_UB lies between 0.45 and 0.70 at mid prices, and the hedged mean is insignificant everywhere.
- **Execution costs.** The unhedged mean is not significant even at mid prices (t = 1.47). Quoted composite spreads on forwards and spot reduce it from 17.7 to 5.4 bp, and option spreads make the hedged mean negative for every k.
- **Payoff timing.** With the previous day's spot, the hedged mean falls by 3.1 bp. The change comes from 43 months with an option in the money on one of the two days, and no single month contributes more than 0.5 bp to the mean. The 10:00 New York cut lies between the two closes in time, and neither close is the cut, so the payoff term is reported under both.
- **Smile construction.** Vanna–volga smiles give φ within 5 per cent of SABR in each regime and the same E1, E2 and E3 conclusions. The vanna–volga price smile is calibrated in all 1,440 currency-months under both readings; it has 10Δ strikes in all of them and passes the static-arbitrage check over ±4 ATM standard deviations in 99.7% (market reading). No leg in the variant uses a smile that fails the check. On the held-out 10Δ quotes it is closer than SABR: mean absolute errors of 0.10 and 0.13 volatility points for the 10Δ risk reversal and butterfly, against 0.13 and 0.16 for SABR (market reading).
- **Butterfly reading, contributor, stale quotes, exclusions, log returns.** None changes a conclusion. The smile-strangle reading and Fenics quotes move φ by less than 4 per cent.

## Three-month tenor

Quarterly, non-overlapping positions: 3L/3S on the three-month forward discount, protected at the 10Δ strike of the three-month market-reading SABR smile, and valued at the spot on the three-month expiry. The primary rebalancing dates are the calendar-quarter month-ends from June 2013; the two other quarterly phases are supplementary.

- **Data.** Three-month forward points and rates were retrieved on 24 September 2026 (research log). Rates are USD OIS for USD and deposits otherwise, as at one month. The three-month quotes are complete for all nine currencies from May 2013; 11 of 11,520 quote fields are substituted and 2 missing.
- **Calibration.** 1,436 of 1,440 smiles converge under each reading, and every converged smile passes both static-arbitrage checks.
- **Timing.** An expiry can fall up to four days after the next rebalancing date, because delivery follows the spot convention. This is also the case at one month.
- **Units.** Returns, C_skew, FD and the skew term are in basis points per quarter. The E3 split reports the skew term only.
- **Inference.** As at one month: Newey–West t-statistics, the stationary-bootstrap interval and the Stambaugh bootstrap, each with 1,999 draws.

| Phase (quarter-end months) | Quarters (φ; zero-rate, hiking) | Δφ (s.e.) | Bootstrap 95% | b (p) | HML^U | HML^H | Skew term | θ_UB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Mar/Jun/Sep/Dec (primary) | 53; 35, 18 | −0.301 (0.166) | [−0.63, −0.04] | 0.0228 (0.045) | 47.5 (1.29) | 36.6 (1.12) | −24.9 | 0.23 |
| Jan/Apr/Jul/Oct | 53; 34, 19 | −0.275 (0.167) | [−0.63, −0.03] | 0.0142 (0.126) | 57.8 (1.64) | 42.8 (1.25) | −23.6 | 0.26 |
| Feb/May/Aug/Nov | 54; 35, 19 | −0.246 (0.130) | [−0.52, −0.04] | 0.0119 (0.169) | 43.0 (1.18) | 9.7 (0.26) | −22.7 | 0.77 |

- **E1.** φ falls from about 0.6 to 0.3 in every phase. The bootstrap interval excludes zero in all three phases; the Newey–West t-statistic is between −1.65 and −1.89. The difference is therefore not significant under both inferences, as at one month. In the post hoc split, the skew price per quarter changes little (25.5 to 23.6 bp in the primary phase), while the forward-discount spread rises from 55 to 82 bp.
- **E2.** The bias-corrected slope is positive in every phase and significant at 5% only in the primary phase (p = 0.045).
- **E3.** The skew term is −23 to −25 bp per quarter (t between −14 and −24); as at one month, its standard error reflects the stability of an ex-ante premium. θ_UB ranges from 0.23 to 0.77 across phases, and its test-inversion set is unbounded in each, so the three-month hedge-cost ratio is not identified.
- **Reading.** The three-month evidence repeats the one-month pattern: a stable skew price, a wider carry spread in the hiking regime, and a regime difference in φ that one inference supports and the other does not.
