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
- **USD deposit rate added.** `USD1MD=` was added to the catalogue so that the forward-point check compares deposit rates with deposit rates. OIS RICs do not exist for EUR and SEK. (Corrected on 24 September: USD discounting uses Fed Funds OIS as designed; the other currencies use deposit rates. See the entry of that date.)
- **Remaining for Stage 1** at that point: checks 6 and 7, a second retrieval to check for quote revisions, and the snapshot time of daily history.

## 2026-09-24

- **Snapshot time.** Daily history was matched against intraday bars over windows in January, March and August to September 2026. The composite daily value corresponds to about 21:18 UTC in both summer and winter, that is, a fixed UTC time rather than a New York local close. The Fenics daily value corresponds to about 17:16 London time. Option payoffs and month-end alignment use these times; the sensitivity check in the design (same-day and previous-day close) is unchanged.
- **Check 6 (conventions).** The instrument metadata name only the quote type and the delta. LSEG's surface-construction documentation gives defaults for its own surfaces (spot delta, premium-adjusted) but does not define the contributed quotes. Delta and premium conventions therefore follow market convention (Clark, 2011; Reiswich and Wystup, 2010), as configured.
- **Check 7 (butterfly reading).** Provider documentation is silent. Under the pre-registered order of evidence, the market-strangle reading is primary. The supporting diagnostic points the other way: across the primary sample, smiles calibrated under the smile-strangle reading predict the held-out 10Δ butterflies better in most month-ends and every currency, with small absolute differences. This is recorded as evidence against the primary reading. E4, which recomputes E1 and the skew term of E3 under the smile-strangle reading, is therefore reported alongside the primary results, not as a footnote.
- **Stage 2 (panel calibration).** Every month-end smile in the primary sample was calibrated under both readings: all converged, with residuals below 1e-13 and all smiles free of static arbitrage within ±4 and ±10 ATM standard deviations. Counts and diagnostics are in the private results. Stage 2's completion criterion is met.
- **Bug fixed during calibration.** The market-strangle butterfly bracket in `quotes_from_smile` could test a negative strangle volatility for low-volatility pairs, which sent the premium-adjusted delta maximiser into an unbounded loop. Volatility is now validated, bracket expansions are bounded and the butterfly bracket keeps σ_ATM + BF positive. Regression tests were added. The delta-to-strike search now brackets locally before scanning and evaluates the smile on the whole grid at once.
- **Remaining for Stage 1:** a second retrieval to check for quote revisions.
- **Citations verified.** Every cited source was checked for its bibliographic details (Crossref, OpenAlex, JSTOR, Project Euclid, publisher pages). Each specific claim attributed to a source was checked against the text, and [references.md](references.md) records what each source is used for and which version was read. Corrections made:
  - The entropy bound is attributed jointly to Bansal and Lehmann (1997) and Alvarez and Jermann (2005), as stated in Backus, Chernov and Zin (2014).
  - The remark on strict-local-martingale exchange rates now states that the foreign measure exists but is not equivalent (Carr, Fisher and Ruf, 2014).
  - Brandt, Cochrane and Santa-Clara (2006) are cited for a risk-sharing index, not a correlation.
  - Delta conventions are cited to Reiswich and Wystup (2012), whose working-paper version was read, and Clark (2011) only as reported there.
  - Brunnermeier, Nagel and Pedersen is dated 2008 (NBER Macroeconomics Annual 2008).
  - The bootstrap conditions are cited to the technical-report version of Politis and Romano.
- **Deviation (before any inference was run): HAC bandwidth.** The design stated a Newey–West lag of ⌊4(T/100)^{2/9}⌋ and attributed it to Newey and West (1994). In that paper this quantity is the pilot truncation of an automatic procedure whose bandwidth is ⌊γ̂ T^{1/3}⌋; using the pilot value directly as the lag is a software convention. The design now uses the full automatic procedure. No standard error had been computed under the earlier rule.

## 2026-09-24 (Stage 3)

- **Stage 3 run.** E1, E4, the 25Δ comparison, the extended sample and the splice test were estimated as specified (`scripts/estimate_e1.py`; method record in [reports/e1.md](reports/e1.md)). The results are held in `data/private/results/` until the licence terms for publishing LSEG-derived findings are confirmed. The ATM comparison named in the design is identically zero by construction, so it is not reported.
- **Independent review before the results were recorded.** A separate review of the estimand and inference code found no error affecting the primary estimates. An independent QuantLib rebuild of every leg agreed to machine precision, and the standard errors matched the R package sandwich exactly. It found five departures, all corrected before the final run; each had a negligible effect on the primary estimates:
  - The extended sample now requires the calibration set only, as logged on 23 September.
  - Option dates are built on the New York calendar, and volatility time is trade to expiry.
  - Provider days-to-maturity is used where reported, and computed spot-to-delivery days otherwise, instead of a 30-day fallback.
  - Currencies, not whole months, are excluded when inputs are missing; non-converged smiles are used and counted.
  - The splice test compares skew costs leg by leg, including the pre-2013 overlap, instead of rows that held by construction.
- **Deviation: discount rates.** The design allows OIS "where available" for non-USD currencies. The implementation uses one-month deposit rates for every non-USD currency, for consistency across currencies. Rates enter E1 only through the base-currency discount factor in spot delta. The review bounded the effect on that factor at about 0.1 per cent, which is negligible for the strikes.
- **Supplementary analyses (post hoc).** Added after the first E1 results, because the series is highly persistent and its zero-rate distribution has a heavy right tail from months with very small forward discounts. With so few effective observations, the automatic Newey–West bandwidth is capped well below what the persistence would require, so the pre-registered intervals are likely too narrow. The additions, listed in [reports/e1.md](reports/e1.md), are reported alongside the pre-registered estimand and do not replace it.
- **Next:** Stage 4 (portfolio returns, E2, E3 and E5).

## 2026-09-24 (Stage 4, part A)

- **Licence.** The licence holder advised that LSEG data are for individual study and research, may not be redistributed, and may not be placed in a public repository. Publishing the extraction code is acceptable. The repository therefore publishes code, methods and aggregate results, and keeps LSEG data, month-level series and calibrated parameters in `data/private/`.
- **Sources read before the Stage 4 results.** Clark and West (working-paper version, Section 2) for the MSPE-adjusted statistic. Stambaugh (1999, published version, eq. 18) for the bias approximation. Lustig, Roussanov and Verdelhan (2011, published version) for the portfolio construction. Background sources on which no computation depends are deferred to the report; [references.md](references.md) marks them.
- **Correction before recording: the E2 bootstrap.** The first implementation resampled blocks of (return, predictor-innovation) pairs. Adjacent pairs inside a block reintroduced the sample's predictive relation into the null distribution, which inflated the estimated bias about tenfold relative to Stambaugh's approximation. The bootstrap now draws pairs of unrestricted-regression residuals and AR(1) residuals independently, and its bias agrees with the approximation.
- **Correction before recording: the E3 reference θ₀.** It is now the average of the hedges' forward deltas over legs. Leg weights had summed to two.
- **Clarification.** The design names the Clark–West test without specifying the standard error. Clark and West recommend the least-squares standard error for one-step forecasts, so that is primary and Newey–West is reported alongside.
- **Results.** Aggregate results are in [reports/e1.md](reports/e1.md) and [reports/stage4.md](reports/stage4.md). Under the pre-registered rules:
  - The E1 regime difference is not significant.
  - The bias-corrected E2 slope on φ is positive (one-sided p = 0.021).
  - The out-of-sample Clark–West statistic is 1.640, against a 5% critical value of 1.645.
  - The regime-shift claim therefore has partial support only.
  - In E3, the realised hedge cost is almost entirely the skew term, and θ_UB is weakly identified.
- **Remaining for Stage 4:** E5, which needs individual VIX futures contracts, not the continuation series acquired so far.

## 2026-09-24 (Stage 4, part B: E5)

- **VIX data.** Individual VIX futures contracts were taken from Cboe's historical settlement files, one per contract from January 2007 to December 2026, stored privately (`scripts/acquire_cboe_vx.py`). Before 26 March 2007 the archive quotes contracts at ten times the index level. The change of scale is detected in the files themselves and removed, including for contracts that expired before that date. After this, Cboe second-month settlements equal the LSEG second-month continuation at every month-end.
- **VIX roll-down.** Caballero and Doyle (2012) short the VIX future whose expiry matches the one-month forward's maturity and hold it to expiry. The design's factor is a month-end analogue: the second-month contract, held over the carry return window from the month-end to the FX option expiry.
- **Correction before recording: Stage 4 return windows.** The data check exposed an error. The last E1 month-end's return window ends after the last spot observation, and the spot look-up had substituted an earlier quote. Windows that end after the data are now excluded, and substitution is limited to five New York business days as designed. The primary return sample is 159 months. The E2 and E3 estimates changed only in the third significant figure, and [reports/stage4.md](reports/stage4.md) reports the corrected values.
- **E5 result.** Unhedged carry loads strongly on the VIX roll-down. The 10Δ hedge lowers the loading by about a fifth. After adjusting for the VIX and FX-volatility factors, neither portfolio earns a significant α, and the exposure-adjusted hedge-cost ratio is unidentified.
- **Stage 4 complete.** Next: Stage 5.

## 2026-09-24 (Stage 5: proofs)

- **Proofs completed.** The remaining results are now proved in `theory/notes.tex`:
  - R3 (spanning, with the integrability condition for taking expectations term by term);
  - R4 (bounds on the unobserved tails from observed wing prices and a tail exponent);
  - R5(c) (monotonicity of the smile delta, and hence injectivity of delta-to-strike inversion);
  - R8 (smoothness of the Hagan volatility for β = 1 and local well-posedness of the smile calibration by the implicit function theorem);
  - the Garman–Kohlhagen PDE (verification, and uniqueness among solutions of polynomial growth).
- **Supporting checks on the panel.** The R5(c) condition holds over ±4 ATM standard deviations at every calibrated month-end, with a maximum of 0.76 against the bound of 1. Premium-adjusted deltas are monotone on the conventional branch everywhere. The Jacobian condition numbers confirm the non-singularity hypothesis of R8 at every calibration.

## 2026-09-24 (Stage 5: sources and robustness)

- **Sources read.** Brunnermeier, Nagel and Pedersen (2008), Burnside et al. (2011), Chernov, Graveline and Zviadadze (2018), Fan, Londono and Xiao (2022), Lustig, Roussanov and Verdelhan (2014) and Bakshi, Kapadia and Madan (2003) in full; Jurek (2014), Farhi et al., Della Corte, Ramadorai and Sarno (2016) and Du, Tepper and Verdelhan (2018) in working-paper or manuscript versions; Choi and Suh (2022) in abstract and introduction only. [references.md](references.md) records the version and sections of each.
- **Corrections from the reading.**
  - R6: Farhi et al. (version of 12 March 2015, Section 5.1) state the leading-order form of the diffusive null, and Jurek (2014) notes that an unlevered hedge gives up part of the diffusive premium. The notes had described R6 as my own without qualification. They now credit both, and claim only the exact mean-value form and the error bound.
  - Design, section 6: θ_UB is an upper bound in Jurek's sense only after allowing for the diffusive premium, so it is read against θ₀. The wording is corrected; the procedure is unchanged.
  - Jurek's option-hedged sample is 1999–2012, not 1990–2012; 1990–2012 is his unhedged sample.
  - Du, Tepper and Verdelhan quote foreign currency per USD, so x_t in R7 is minus their basis. The design and notes now say so.
  - R3: the moment contracts of Bakshi, Kapadia and Madan expand around the spot price; the notes state that the R3 contracts are centred at the forward.
- **Robustness grid.** `scripts/robustness.py` estimates every variant of design section 8 except the three-month tenor and vanna–volga smiles. The rules the design leaves open (stale-quote rule, exclusion windows, dollar-carry signal, the USD leg, log hedged returns, the out-of-sample training start) were fixed in the script before any variant was estimated. The base variant reproduces E1 and Stage 4. Results are in [reports/robustness.md](reports/robustness.md):
  - φ is lower in the hiking regime in every variant; the difference is significant at 5% only in the ten-currency ranking, and the skew price per month does not change across regimes.
  - The in-sample E2 slope survives every HML variant; the out-of-sample test does not reject in any variant when training starts in 2013.
  - Quoted spreads remove most of the unhedged mean, and previous-day spot lowers the hedged mean by 3.1 bp.
- **Implementation.** `option_dates` now builds the New York holiday table once; dates are unchanged and the tests pass. This reduced the grid's run time from ten minutes to under half a minute.
- **Remaining for Stage 5:** the three-month tenor (needs three-month forwards and rates), vanna–volga smiles, and the R4 moment intervals.

## 2026-09-24 (Stage 5: choices fixed before the remaining variants)

Recorded before any of the estimates below were computed.

- **Three-month data.** Three-month forward points and rates (USD3MOIS=, USD3MD=, USDSROIS3M=, and <CCY>3MD= for the nine currencies) were retrieved on 24 September 2026 into a separate private retrieval folder. Volatility quotes come from the 23 September retrieval.
- **Three-month variant.**
  - Rebalancing and returns: quarterly and non-overlapping, at calendar-quarter month-ends from June 2013.
  - Discount rates: OIS for USD and deposits otherwise, as at one month.
  - Legs and options: 3L/3S on the three-month forward discount, with protective options at the 10Δ strike of the three-month market-reading SABR smile.
  - Estimands: E1, in-sample E2 and E3 as at one month.
  - Supplementary: the two other quarterly phases are reported alongside.
- **Vanna–volga variant.**
  - Pillars: the ATM delta-neutral-straddle strike, and the 25Δ put and call strikes at the pillar volatilities in the pair's delta convention.
  - Butterfly: under the market reading, the smile strangle is solved so that the vanna–volga smile prices the market strangle. Under the smile reading, the quote is used directly.
  - Method: the primary vanna–volga method is chosen from the source reading before estimation. Months where a method is undefined at the 10Δ strikes are counted, not dropped silently.
- **Option-implied moments (R4, secondary E2 predictors).**
  - Measure: moments of log(X/F_X) under the USD forward measure; for USD-base pairs, prices are converted with the change of numeraire of R1.
  - Integration: the calibrated smile between the 10Δ strikes, with R4 tail bounds beyond them.
  - Tail exponents, two settings: (i) the local elasticities of the calibrated smile's density at the 10Δ strikes, which assumes the tails beyond the quotes are no heavier than at the quotes; (ii) γ = η = 2.
  - Predictors: portfolio variance and oriented skewness, averaged over the E1 legs.
  - Regressions: at the lower endpoint, midpoint and upper endpoint of each interval under setting (i). A conclusion is reported only if it holds at all three.

## 2026-09-25 (Stage 5: vanna–volga, three-month tenor, moment intervals)

- **Sources read for vanna–volga.** Castagna and Mercurio (2007; authors' post-review preprint), their working paper (September 2006 revision, compared with the January 2006 version), Bossens et al. (2010; arXiv v3) and Reiswich (2010, doctoral dissertation). Reiswich and Wystup (2012) was checked and does not treat vanna–volga. [references.md](references.md) records versions and locations.
- **Choice fixed from the reading, before estimation.** The primary vanna–volga smile is the implied volatility of the vanna–volga price, not either approximation. The properties Castagna and Mercurio prove belong to the price. The first-order approximation overvalues the wings. The second-order approximation can be undefined. The approximations are computed as checks.
- **Rule added in implementation.** A 10Δ strike on a vanna–volga smile must be connected to the ATM pillar through strikes where the smile is defined. Across an undefined region the delta map is not continuous. The rule is not from a source, and it does not bind for the price smile's strikes used in the variant.
- **Independent review.** Each of the three pieces was reviewed separately after implementation, and confirmed defects were fixed before any result was recorded. None of the fixes changed an estimate.
  - Vanna–volga:
    - missing diagnostic flags had been counted as passes;
    - the held-out comparison had included non-converged smiles;
    - a citation range was wrong;
    - the strike-set invariance test was missing.
  - Three-month tenor:
    - the holiday-table guard did not bound the delivery date;
    - a reproduction command was missing.
  - Moment intervals:
    - a closed-form check had been applied outside its stated domain;
    - a test case did not exercise the edge extremum it claimed to;
    - one docstring statement was wrong.
- **Results.** Vanna–volga and three-month results are in [reports/robustness.md](reports/robustness.md); the moment predictors are in [reports/stage4.md](reports/stage4.md).
  - Vanna–volga: smiles leave E1, E2 and E3 unchanged, and fit the held-out 10Δ quotes better than SABR.
  - Three-month tenor: it repeats the one-month pattern. φ falls in the hiking regime because the carry spread widens while the skew price is stable. The difference is significant under the bootstrap but not under Newey–West.
  - Moment intervals:
    - Under γ = η = 2 no skewness interval excludes zero, so truncated quotes do not identify skewness.
    - Under the boundary-elasticity setting, oriented skewness predicts carry returns with the crash-compensation sign under the bootstrap at every endpoint.
    - That setting is contradicted by the smile's own wings in about three quarters of currency-months, so the result is conditional supporting evidence only.
- **Open.** USD3MOIS= is used as USD OIS at three months; the 24 September retrieval has no provider metadata confirming its description.
- **Remaining:** the C++ kernel, the Stage 1 quote-revision check, and Stage 6.

## 2026-09-25 (licence compliance audit)

- **Trigger.** A review of LSEG's published guidance on derived data. LSEG's data-redistribution page treats an output as derived data if it is unrecognisable, non-reversible and cannot be traced back to the original content without exceptional effort. The London Stock Exchange market data policy (Schedule B, 2026) defines derived data as data that cannot be reverse engineered back to the data or used as a substitute for them. Both documents govern commercial and exchange licences, not the academic Workspace licence; the binding terms remain the licence holder's guidance (entry of 24 September). The principle agrees with the adopted policy: publish statistical summaries only.
- **Audit.** Every tracked file and every commit in the history was examined for LSEG-derived content beyond statistical summaries. No data file and no path under `data/private/` has ever been tracked. Findings and actions:
  - A regression-test fixture in `tests/test_smile.py` used rounded SABR parameters, a forward and a discount factor that match one calibrated currency-month. It is replaced by synthetic values, which reproduce the original failure on the code before the fix.
  - A verbatim instrument-metadata label in this log is paraphrased.
  - A single month-level value in the robustness report is replaced by a pooled statement.
  - **History rewritten.** The public history was rewritten on 25 September 2026 so that no commit contains the fixture values, the metadata label or the month-level value. Each affected commit now carries the replacement wording, the current files are unchanged, and the commits from the Stage 2 calibration onwards have new identifiers.
  - The data plan and the calibration record said coverage facts were kept private, while the reports publish them. The policy is now stated consistently: coverage facts (sample windows, first available dates, counts) are published; quote values, month-level series and calibrated parameters are not.
- **Safeguards.** `.gitignore` now excludes spreadsheet and binary data formats, and `tests/test_repository_hygiene.py` (run in CI) fails if a data extract, anything under `data/private/` or a notebook with stored outputs is tracked.
- **Open.** Confirmation from the licence holder that pooled statistical summaries and coverage facts may be published.

## 2026-09-25 (Stage 6: paper draft)

- **Paper.** I drafted [paper/paper.tex](paper/paper.tex) (23 pages, main text 18) and [paper/summary.md](paper/summary.md). Three independent reviews followed: one for numbers and licence, one for citations and originality, and one referee review of the economics, inference and interpretation rules. Each finding was checked against the record before it was applied; 49 were applied and 9 rejected with reasons.
- **Corrections to the reports found by the reviews.**
  - The primary sample starts at the 31 May 2013 month-end, not June 2013. The 159 return month-ends run from May 2013 to July 2026, and the extended sample starts in January 2007. I corrected e1.md and stage4.md.
  - The split of φ into C_skew and FD is a post hoc supplementary analysis. The E1 table and the reports now label it so.
  - The skew term of E3 is an ex-ante premium, so its standard error measures time variation, not the precision of realised crash compensation. The reports no longer call it precisely estimated.
  - e1.md listed supplementary analyses without their results; it now reports all of them.
- **Added, post hoc.** Realised carry returns by regime. The unhedged mean is 7.4 bp per month before 2022 and 37.0 after, a difference of 29.6 bp (s.e. 23.8) that is not significant; the hedged difference is also not significant. Stage 4 report.
- **Interpretation rule.** I require both the Newey–West and the bootstrap interval to exclude zero before an E1 difference counts as significant. This is my reading of the design's rule, stated as such in the paper.
- **Open.**
  - Attribute the test-inversion set for θ_UB to Fieller's method once the source is read.
  - Clean-environment reproduction.
  - The licence holder's confirmation on publishing pooled summaries.
  - The validation checks not yet run: the published-portfolio comparison and the Merton and Heston checks of the moment code.
