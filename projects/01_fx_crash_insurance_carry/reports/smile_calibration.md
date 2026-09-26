# Smile calibration and butterfly diagnostic

This record covers the Stage 2 panel calibration and Stage 1 check 7, run on the primary sample, and gives the method and the outcomes that set a method. Calibrated parameters, currency-month errors and all quote values are LSEG-derived and stay in `data/private/results/`; only pooled counts and summary statistics are published here.

## Reproduction

```bash
.venv/bin/python scripts/calibrate_smiles.py --start 2010-07-30 --end 2026-08-31
```

## Inputs

Month-end inputs are built by `src/qef/data/smile_inputs.py` under the five-business-day rule of the research design (section 4). Spot is the mid quote, and the one-month outright forward is built from the mid forward points and the pair's pip factor. The one-month quote-currency rate (Fed Funds OIS for USD, deposit rates otherwise) gives the discount factor D_q over the forward's spot-to-delivery period on the currency's money-market basis, and the base-currency factor follows as D_b = F D_q / S. Volatility time is (expiry − trade)/365, with the expiry date equal to the delivery date less the spot lag on the New York calendar; an earlier version used the spot-to-delivery period and was corrected on 24 September 2026 (research log). The quotes are the mid ATM, 25Δ and 10Δ risk-reversal and butterfly quotes.

## Method

For each currency, month-end and butterfly reading (market strangle or smile strangle), a SABR smile with β = 1 (Hagan, Kumar, Lesniewski and Woodward, 2002) is calibrated to the ATM volatility, the 25Δ risk reversal and the 25Δ butterfly. The system of theory result R8 is solved by Levenberg–Marquardt in (ln α, atanh ρ, ln ν), starting each month from the previous month's solution with three fixed starting points as fall-backs. Convergence, defined as a maximum absolute residual below 1e-8 volatility units, is recorded together with the Jacobian condition number and the strikes.

Static arbitrage is checked (R2) on log-moneyness grids of ±4 and ±10 ATM standard deviations, using discrete convexity of call prices and the density condition of Gatheral and Jacquier (2014). For the butterfly diagnostic, the calibrated smile predicts the 10Δ risk reversal and butterfly under the same reading, and the absolute prediction errors of the two readings are compared month by month.

## Outcomes

Every month-end smile converged under both readings and passed the static-arbitrage checks on both grids. The smile-strangle reading predicts the held-out 10Δ butterflies better in most month-ends and in every currency, although the absolute differences are small. Because the provider documentation does not define the quote, the pre-registered rule keeps the market-strangle reading as primary and reports E4 alongside the primary results (research log, 24 September 2026).
