# Research design

Version 1.0, 23 September 2026, fixed before any return, hedge cost or option-implied statistic was computed. Changes are logged in [research_log.md](research_log.md) with their date and reason, and results under the original rule remain reported. On 26 September 2026 the text was re-edited for presentation only, with no change to any rule, estimand or threshold; earlier versions are in the Git history.

## 1. Question

Does the ex-ante price of crash insurance in G10 FX option smiles explain the collapse of the currency carry premium after 2008 and its revival in 2022–2026? How much of the realised cost of option-hedging carry is due to skew pricing, and how much to the option's delta and the volatility level?

## 2. Motivation and identification

Carry investors earn the interest differential, and they lose when high-yielding currencies fall sharply in bad states (Brunnermeier, Nagel and Pedersen, 2008). Let M and M* be discount factors that price one-period bonds and FX consistently in the US and abroad. The expected log excess return on foreign currency is then

E_t[s_{t+1} − f_t] = L_t(M) − L_t(M*) − x_t,

where L_t(M) = log E_t M − E_t log M is conditional entropy and x_t is the covered-interest-parity basis (theory result R7). x_t is minus the basis of Du, Tepper and Verdelhan (2018), who quote exchange rates in foreign currency per USD. Entropy is a sum of cumulants of log marginal utility, so crash risk enters through cumulants of order three and above.

Options on S identify the risk-neutral law of Δs = m* − m, but they do not identify how the entropy difference divides across cumulant orders. The project therefore measures the price of the carry position's left-tail exposure; it does not measure the crash share of discount-factor entropy. A residual premium after hedging can reflect diffusive risk, systemic risk that FX options insure cheaply (Caballero and Doyle, 2012), or peso-state pricing (Burnside, Eichenbaum, Kleshchelski and Rebelo, 2011).

The closest work is Jurek (2014, *JFE*), who studies crash-neutral carry with unhedged returns over 1990–2012 and option-hedged returns over 1999–2012 and finds that crash premia account for at most one third of carry returns. Related studies are Burnside et al. (2011, *RFS*), Chernov, Graveline and Zviadadze (2018, *JFQA*), Della Corte, Ramadorai and Sarno (2016, *JFE*), Fan, Londono and Xiao (2022, *JFE*) and Choi and Suh (2022, *JIFMIM*), together with the working papers of Farhi, Fraiberger, Gabaix, Rancière and Verdelhan (NBER WP 15062, 2009; SSRN version of 12 March 2015) and Kutuk and van Wijnbergen (CEPR DP 20745, USD/TRY only).

The intended contribution is a decomposition of G10 crash-insurance cost into delta, volatility-level and skew components across two rate regimes, and a test of whether the ex-ante skew price explains the change in the carry premium between those regimes. A secondary contribution is the bias from misreading the butterfly quote convention.

## 3. Notation and conventions

The currencies are J = {EUR, GBP, AUD, NZD, JPY, CHF, CAD, NOK, SEK}, all against USD. X_t^j is the USD price of one unit of j, so USD-base quotes (USDJPY, USDCHF, USDCAD, USDNOK, USDSEK) are inverted. F_t^j is the one-month outright forward, spot plus points times the pip factor. With s = log X and f = log F, the forward discount is fd_t^j = s_t^j − f_t^j. If an NOK or SEK volatility instrument is quoted against EUR, the currency is triangulated through EURUSD or excluded; that decision is taken in the data audit, before any return is computed.

Excess returns are per USD of forward notional. The arithmetic return is rx_{t+1}^j = X_{t+1}^j / F_t^j − 1 and the log return is s_{t+1}^j − f_t^j. Estimands use arithmetic returns, R7 concerns log returns, and both are reported.

Pricing is Garman–Kohlhagen in forward form. The undiscounted call is F Φ(d₊) − K Φ(d₋), discounted with the quote-currency factor D_q, and the base-currency factor is implied by the forward, D_b = F D_q / S, so no separate covered-interest-parity assumption is made. USD discounting uses one-month Fed Funds OIS, with SOFR OIS as a check from its start; other currencies use one-month deposit rates, or OIS where available.

The quote conventions follow Reiswich and Wystup (2012), who report those of Clark (2011), and are confirmed per instrument in the audit: spot delta up to one year; premium-adjusted delta for USD-base pairs; the ATM strike as the delta-neutral straddle, K = F exp(+σ²τ/2) for pips delta and F exp(−σ²τ/2) for premium-adjusted delta; and the risk reversal as σ(call) − σ(put) on the quoted pair.

A three-parameter smile fits either butterfly reading exactly, so the quotes cannot decide between them. The reading is set, in order, by (i) provider documentation, (ii) the interdealer market-strangle convention and (iii), as supporting evidence only, which reading better predicts the unused 10Δ quotes.

The smile is SABR with β = 1, calibrated to the ATM volatility, the 25Δ risk reversal and the 25Δ market strangle by solving the system of R8 with Levenberg–Marquardt. The 10Δ quotes are held out as a fit test, and vanna–volga is the alternative. Every smile is checked for static arbitrage (R2), and non-convergence and arbitrage flags are counted, not dropped.

## 4. Samples

Samples are defined by rules; the resulting dates are LSEG-derived and are recorded in the private audit. The primary sample uses composite quotes, from the first month-end at which all nine currencies have two-sided one-month ATM, 25Δ and 10Δ risk reversals and butterflies, and forward points, with three long and three short legs. The extended sample uses Fenics (`=FN`) quotes for the currencies available before the primary start, with two long and two short legs; it covers 2008, and a splice test compares Fenics and composite quotes on their overlap. The long ATM sample is ATM-hedged carry (Burnside et al., 2011) from the start of the ATM and forward series.

Sampling is at New York month-ends. A missing two-sided quote is replaced by the last two-sided quote within five business days; otherwise the currency is excluded that month, and both events are counted. A butterfly unchanged for five or more business days is flagged, and treated as missing in a robustness run. Option payoffs use end-of-day spot on the expiry date; WM fixings and the 10:00 New York cut are not available, so same-day and previous-day closes are compared.

## 5. Portfolios and hedges

HML_FX follows Lustig, Roussanov and Verdelhan (2011). At month-end t the currencies are ranked by fd_t^j, and the portfolio goes long 1/3 USD of forward notional in the top three and short 1/3 in the bottom three, so that

HML^U_{t+1} = (1/3) Σ_long rx^j − (1/3) Σ_short rx^j,   FD_t = (1/3) Σ_long fd^j − (1/3) Σ_short fd^j ≥ 0.

The alternatives are dollar carry (Lustig, Roussanov and Verdelhan, 2014) and a ten-currency ranking in which USD can be a leg.

Each leg buys a one-month option that pays when the leg loses: a put on j for a long leg, which is a USD call for USD-base pairs, and a call on j for a short leg. The notional matches the forward, and the strike is the market 10Δ strike in the pair's own convention (R5). 10Δ is primary, because a one-month 25Δ option is only about 2% out of the money, and 25Δ and ATM are comparisons. The hedged return is

HML^H_{t+1} = HML^U_{t+1} + Σ_legs (1/3)[payoff_{t+1} − premium_t (1 + r_t^{USD} τ)],

in USD per USD of notional at t + 1. Estimands use mid prices. The implementable version trades forwards at bid and ask and sets the option half-spread to k × the ATM half-spread in volatility, k ∈ {1, 1.5, 2}, converted to premium by vega, because package quotes do not give a single-option ask.

## 6. Estimands

E1 is the primary estimand and is measured ex ante. For each leg, V^smile is the 10Δ hedge premium from the calibrated smile and V^flat the premium at the same strike with the ATM volatility, both in USD per USD of notional at t + 1. Then

C_t^skew = Σ_legs (1/3)(V^smile − V^flat),   φ_t = C_t^skew / FD_t.

E1 reports φ_t, its means in the zero-rate regime (primary start to December 2021) and the hiking regime (January 2022 onwards), and their difference with HAC and bootstrap 95% intervals. The 25Δ and ATM versions are comparisons.

E2 concerns predictability, through HML^U_{t+1} = a + b φ_t + e_{t+1}. The alternative predictor is the oriented 10Δ risk reversal, the volatility of the protective option minus the opposite option, averaged over legs. The out-of-sample test uses an expanding window with a first forecast for January 2017 and initial training on the extended sample, and the Clark and West (2007) test against the historical mean. Option-implied variance and skewness are secondary predictors; they depend on tail extrapolation and are reported as R4 intervals.

E3 is measured ex post. The realised excess payoff of each option, π = payoff − premium × (1 + rτ), is split into three terms: (i) payoff − V^GK(σ̂^P)(1 + rτ), the delta and diffusive term; (ii) −[V^flat − V^GK(σ̂^P)](1 + rτ), the volatility-level term; and (iii) −[V^smile − V^flat](1 + rτ), the skew term. Here σ̂^P is the annualised realised volatility of daily log spot changes over the previous 21 business days, and R6 gives the diffusive null θ₀ ≈ Φ(−d₁). E3 reports θ_UB = 1 − μ_H/μ_U, which can take any real value. It bounds the crash share from above in the sense of Jurek (2014), but an unlevered hedge also gives up part of the diffusive premium, so it is compared with θ₀ (Farhi et al., 2015, state the leading-order form of θ₀).

E4 recomputes E1 and E3(3) under the other butterfly reading. E5 regresses HML^U and HML^H on two factors: the VIX roll-down return, the negated monthly change of the second-month VIX future with the same contract within the month; and global FX-volatility innovations, the AR(1) residual of the monthly mean absolute daily log change across the nine currencies (Menkhoff, Sarno, Schmeling and Schrimpf, 2012). It reports θ_UB after adjusting for these exposures. A descriptive account covers smiles, φ_t and returns around 5 August 2024.

## 7. Inference

Standard errors are those of Newey and West (1987) with a Bartlett kernel, and the bandwidth follows the automatic procedure of Newey and West (1994): pilot truncation n = ⌊4(T/100)^{2/9}⌋, then bandwidth m = ⌊γ̂ T^{1/3}⌋. The stationary bootstrap of Politis and Romano (1994) uses the block length of Politis and White (2004) as corrected by Patton, Politis and White (2009), with 9,999 draws. θ_UB is obtained by test inversion, CS = {θ₀ : |μ̂_H − (1 − θ₀)μ̂_U| ≤ 1.96 ŝe_HAC}, the construction of Fieller (1954) for a ratio of means; unbounded sets are reported as unbounded (Dufour, 1997), and delta-method intervals are secondary. Stambaugh bias in E2 is assessed by residual bootstrap under the null. E1 is the single primary estimand, so no multiple-testing adjustment is made.

## 8. Robustness

1. 25Δ and ATM hedges.
2. Three-month tenor with quarterly rebalancing.
3. Fenics in place of composite quotes.
4. Stale butterflies treated as missing.
5. Excluding March 2020, and CHF over December 2014 to March 2015.
6. Half-spread k grid.
7. Dollar-carry and ten-currency portfolios.
8. Log returns.
9. Vanna–volga smiles.
10. Previous-day spot for payoffs.

## 9. Interpretation rules

For E1, a regime difference is claimed only if its 95% interval excludes zero, and its size is read relative to FD_t. The statement that the skew price "explains the regime shift" requires all three of: the E1 difference, b > 0 in E2 after the Stambaugh correction, and an out-of-sample Clark–West rejection at 5%; otherwise the result is reported as partial or no support. In E3, θ_UB is not described as a crash share, and a skew term whose interval includes zero is reported as no evidence of priced crash risk beyond the diffusive and volatility-level terms. Null results are reported with the same prominence as positive ones.

## 10. Validation

Pricing and delta conventions are checked against analytic identities: parity, the delta-neutral straddle and inversion round trips. Smile calibration is checked by synthetic recovery and against an independent open-source implementation (QuantLib). Moment and density code is checked against Garman–Kohlhagen, Merton and Heston closed forms; SABR is not a benchmark, because Hagan's formula is an approximation. HML_FX is compared in sign and magnitude with Verdelhan's public portfolios, and θ_UB is compared only with premium shares in the literature, noting differences in hedge moneyness.

## 11. Stages and completion criteria

1. Design and data audit. The design is frozen, the LSEG panel acquired, and coverage, two-sidedness, staleness, underlying pairs and conventions audited. The stage is complete when the private audit and the public audit record exist.
2. Pricing and smile engine. Complete when the analytic checks and synthetic recovery pass, the QuantLib cross-check agrees, and calibration and arbitrage flag rates across the panel are recorded.
3. E1 and E4. Complete when the regime estimates with intervals and the splice test are done.
4. Portfolios, E2, E3 and E5. Complete when every estimand is reported as specified, including null and unbounded results.
5. Robustness, proofs of R1–R9, and the C++ kernel with parity against the Python reference. Complete when every claim traces to a derivation, a citation or a reproducible run.
6. Paper and release. Complete when clean-environment reproduction passes and publication terms for LSEG-derived outputs are confirmed.
