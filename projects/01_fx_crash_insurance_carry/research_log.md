# Research log

## 2026-09-23

- **Design version 1.0 fixed.** No returns, hedge costs or option-implied statistics had been computed.
- **Butterfly reading.** A three-parameter smile reproduces either butterfly reading exactly, so the reading is set by provider documentation, then by the interdealer convention. The 10Δ fit is supporting evidence only.
- **Primary hedge moneyness is 10Δ.** At one month, 25Δ strikes lie about 2% out of the money.
- **Payoff spot.** WM fixings are not licensed. Option payoffs use end-of-day spot, with a sensitivity check.
- **Independent cross-check.** QuantLib replaces the MATLAB implementation first considered, because it is open source and reproducible without a licence.
- **Pricing and smile engine implemented:**
  - Garman–Kohlhagen pricing in four delta conventions.
  - Delta inversion, including the premium-adjusted branch rule.
  - SABR calibration under both butterfly readings.
  - Static-arbitrage checks.
  - Tests: parity, delta-neutral straddle, inversion round trips, and synthetic SABR recovery under both readings and both delta types.
- **QuantLib cross-check** (QuantLib 1.43, independent C++ implementation). Premiums, deltas in the spot, forward and premium-adjusted conventions, delta-to-strike inversion, delta-neutral ATM strikes and Hagan SABR volatilities agree to between 1e-9 and 1e-13 on the test grid (`tests/test_quantlib_crosscheck.py`).

## 2026-09-23 (data audit)

- **Acquired and audited.** The LSEG panel was acquired: spot, forward points, ATM volatility, 25Δ and 10Δ risk reversals and butterflies at 1M and 3M for the composite and Fenics contributors, one-month rates and VIX futures. Audit checks 1 to 5, 8 and 9 were run. Counts and dates are in the private audit.
- **Conventions confirmed.** Every volatility instrument, including NOK and SEK, is quoted against USD, so no triangulation is needed. Provider scaling factors confirm the configured pip factors. Risk reversals are call minus put on the quoted pair, and the 10Δ quotes follow the same convention.
- **Clarification: extended sample.** The rule requires the calibration set only (ATM, 25Δ risk reversal, 25Δ butterfly) plus forward points. 10Δ quotes are used for the fit test where available. This was recorded before any estimate was computed.
- **Extended-sample quality.** Fenics quotes before 2010 update infrequently, and butterflies in particular are often unchanged for weeks. Extended-sample results are therefore secondary evidence. The stale-butterfly robustness run applies to them, and the long ATM sample is used alongside them for 2008.
- **USD deposit rate added.** `USD1MD=` was added to the catalogue so that the forward-point check compares deposit rates with deposit rates. OIS RICs do not exist for EUR and SEK; deposit rates are primary for all currencies, as designed.
- **Remaining for Stage 1:**
  - check 6 (delta and premium conventions from provider documentation);
  - check 7 (the butterfly convention);
  - a second retrieval to check for quote revisions;
  - the snapshot time of daily history.
