# Smile calibration and butterfly diagnostic

Stage 2 panel calibration and Stage 1 check 7, run on the primary sample. This record gives the method and the outcomes that set a method. Calibrated parameters, currency-month errors and all quote values are LSEG-derived and stay in `data/private/results/`. Only pooled counts and summary statistics are published here.

## Reproduction

```bash
.venv/bin/python scripts/calibrate_smiles.py --start 2010-07-30 --end 2026-08-31
```

## Inputs

Month-end inputs are built by `src/qef/data/smile_inputs.py` under the five-business-day rule of the research design (section 4):

- **Spot and forward.** Mid spot, and the one-month outright forward from mid forward points and the pair's pip factor.
- **Discounting.** The one-month quote-currency rate (Fed Funds OIS for USD, deposit rates otherwise) gives D_q over the forward's spot-to-delivery period, on the currency's money-market basis. D_b = F D_q / S.
- **Time to expiry.** Taken as the spot-to-delivery period in calendar days divided by 365, since the spot lag before expiry and before delivery are equal.
- **Quotes.** The mid ATM, 25Δ and 10Δ risk-reversal and butterfly quotes.

## Method

For each currency, month-end and butterfly reading (market strangle or smile strangle):

1. **Calibration.** A SABR smile with β = 1 (Hagan, Kumar, Lesniewski and Woodward, 2002) is calibrated to the ATM volatility, the 25Δ risk reversal and the 25Δ butterfly. The system of theory result R8 is solved by Levenberg–Marquardt in (ln α, atanh ρ, ln ν). Each month starts from the previous month's solution, with three fixed starting points as fall-backs.
2. **Recorded diagnostics.** Convergence (maximum absolute residual below 1e-8 volatility units), the Jacobian condition number and the strikes are recorded.
3. **Arbitrage check.** Static arbitrage is checked (R2) on log-moneyness grids of ±4 and ±10 ATM standard deviations. The check uses discrete convexity of call prices and the density condition of Gatheral and Jacquier (2014).
4. **Butterfly diagnostic.** The calibrated smile predicts the 10Δ risk reversal and butterfly under the same reading. Absolute prediction errors are compared across readings, month by month.

## Outcomes

- **Calibration.** Every month-end smile converged under both readings, and every smile passed the static-arbitrage checks on both grids.
- **Butterfly diagnostic.** The smile-strangle reading predicts the held-out 10Δ butterflies better in most month-ends and in every currency. Absolute differences are small. Provider documentation does not define the quote. The pre-registered rule therefore keeps the market-strangle reading as primary and reports E4 alongside the primary results (research log, 24 September 2026).
