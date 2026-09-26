# Theory

This note lists the results used in the [research design](../research_design.md). Every result marked P or P\* is proved in `notes.tex` (compiled as `notes.pdf`), and R9 is cited with its conditions. A result marked P is proved in the notes, and the part I claim as my own is named in its sources, where the related literature is cited. A result marked P\* is a known result, and the notes write out the published or standard proof for completeness. A result marked C is cited with its hypotheses.

## R1 (P\*). Domestic and foreign risk-neutral measures

Let (Ω, 𝔽, (𝔽_t), Q^d) satisfy the usual conditions, let B^d and B^f be strictly positive money-market accounts with B^d_0 = B^f_0 = 1, and let S be strictly positive and càdlàg. If Z_t = S_t B^f_t / (S_0 B^d_t) is a true Q^d-martingale, then dQ^f/dQ^d = Z_T defines a measure Q^f equivalent to Q^d, with E^{Q^f}[X | 𝔽_t] = E^{Q^d}[Z_T X | 𝔽_t] / Z_t for X ≥ 0 or X ∈ L¹(Q^f). Foreign prices discounted by B^f are then Q^f-martingales exactly when their domestic values discounted by B^d are Q^d-martingales, and under Garman–Kohlhagen E^{Q^d} S_T = F and E^{Q^f}(1/S_T) = 1/F, which resolves Siegel's paradox.

If Z is only a strict local martingale, then E Z_T < 1 and Z_T does not define an equivalent measure. Carr, Fisher and Ruf (2014) construct the foreign measure in that case as a Föllmer measure that is not equivalent to Q^d, and restore the domestic-foreign symmetry with a modified pricing operator. The result is the change of numeraire of Geman, El Karoui and Rochet (1995).

## R2 (P\*). Arbitrage-free call prices

A function c satisfies c(K) = E(X − K)⁺ for some X ≥ 0 with E X = F if and only if c is convex, c(0) = F, c′₊(0) ≥ −1 and c(K) → 0 as K → ∞. The law of X is then unique, with distribution function G = 1 + c′₊ and P(X = 0) = 1 + c′₊(0); extending c by c(K) = F − K for K < 0, the law equals c″ in D′(ℝ). On a strike grid the tests used are slopes in [−1, 0] that are nondecreasing, and on total variance w(k) the condition g(k) ≥ 0 of Gatheral and Jacquier (2014). The link between call prices and state prices is due to Breeden and Litzenberger (1978), and the characterisation and its proof through the right derivative and Fubini's theorem are in Föllmer and Schied (2004, Lemma 7.23).

## R3 (P\*). Spanning

If f′ is locally absolutely continuous on (0, ∞), Q(X > 0) = 1 and ∫₀^F |f″|P dK + ∫_F^∞ |f″|C dK < ∞, then E f(X) = f(F) + ∫₀^F f″(K)P(K) dK + ∫_F^∞ f″(K)C(K) dK. The formula is that of Carr and Madan (1998), and the proof is the standard argument by Taylor's theorem with integral remainder and Fubini's theorem, written out with explicit integrability hypotheses. The moment contracts of Bakshi, Kapadia and Madan (2003) are of this form.

## R4 (P). Partial identification from truncated quotes

If Q(S_T ≤ K)/K^γ is nondecreasing on (0, K_min], then P(K) ≤ P(K_min)(K/K_min)^{γ+1} there, and if K^η Q(S_T > K) is nonincreasing on [K_max, ∞) with η > 1, then C(K) ≤ C(K_max)(K/K_max)^{1−η} there; for example, ∫₀^{K_min} P K⁻² dK ≤ P(K_min)/(γ K_min). The same inequalities bound the tails of any weight, including the moment contracts adapted from Bakshi, Kapadia and Madan (2003), with sign-changing weights split into positive and negative parts, so that variance and skewness are reported as intervals over stated tail exponents when the quotes stop at the 10Δ strikes. Between the quoted strikes the integrals use quadrature with F as a node, because the integrand has a kink there; the trapezium error is O(h²) and the Simpson error O(h⁴) on each side.

Jiang and Tian (2005, Proposition 2 and Appendix) bound the truncation error of model-free implied variance beyond the quoted strikes, and Lee (2004) links finite moments to the wings of implied variance. The power-type tail condition, the bounds for general and sign-changing weights and the resulting variance and skewness intervals are my own. The implementation is `src/qef/fx/moments.py`, with R1 for USD-base pairs, and the results are in `reports/stage4.md`.

## R5. Delta-to-strike inversion

For pips spot delta the strike has the closed form K = F exp(σ²τ/2 − φσ√τ Φ⁻¹(|Δ|/D_b)).

(b), marked P\*. The premium-adjusted call delta h(K) = (K/F) D_b Φ(d₋) tends to zero at both ends, and h′ has the sign of Φ(d₋) − φ(d₋)/(σ√τ). Since φ/Φ is strictly decreasing, h is strictly increasing and then strictly decreasing, with a unique maximiser K*, so each attainable delta has exactly two strikes and the convention takes the root in [K*, ∞).

(c), marked P. For pips delta Δ_φ(K) = φ D Φ(φ d₊(K, σ(K))), if |∂σ/∂k| √τ |d₋| < 1 on an interval, with k = ln(K/F), then Δ_φ is strictly decreasing there and the smile delta-to-strike map is injective. The condition holds over ±4 ATM standard deviations at every calibrated month-end, and premium-adjusted deltas are checked numerically.

Reiswich and Wystup (2012; working-paper version of 20 March 2010, Sections 1.3 and 1.5) give the delta definitions, the non-monotonicity of the premium-adjusted call delta, the right-branch convention, the bound by the simple delta and the first-order condition σ√τ Φ(d₋) = φ(d₋) for its maximum. Jäckel (2020, Section 2, eq. (20)) gives the single maximum through the inverse Mills ratio and the choice of the larger root, which is the argument written out for (b); that the inverse Mills ratio is decreasing follows from the log-concavity of Φ (Bagnoli and Bergstrom, 2005, Section 2.1, Theorem 1) and is classical (Sampford, 1953). Reiswich (2010, Section 3.8) shows that a smile delta can fail to be monotone. Fukasawa (2012, Lemma 2.1 and Propositions 2.1–2.2) shows by the same derivative calculation that d₊ and d₋ are strictly decreasing for every arbitrage-free smile, and Mingone (2023, Lemma 4.1) links this to a unique strike for each delta; my own part of (c) is only the explicit slope bound, checkable on a fitted smile not assumed to be free of arbitrage, for pips spot and forward deltas. The implementation is `src/qef/fx/gk.py`.

## R6 (P). Diffusive null of the hedge cost

Under Garman–Kohlhagen with equal P and Q volatility and no jumps, where the P-drift of log S exceeds the Q-drift by λ, E^P(K − S_T)⁺ − E^Q(K − S_T)⁺ = −FλτΦ(−d₁) + o(λτ), and hence θ₀ = Φ(−d₁) + O(λτ). The proof uses the mean value theorem with ∂_G p(G) = −Φ(−d₊(G)) for the undiscounted put. The leading-order form is stated by Farhi et al. (2015, Section 5.1, eqs. (8)–(9)), and Jurek (2014) notes that an unlevered hedge gives up part of the diffusive premium; the exact mean-value form and the error bound, proved in the notes, are my own.

## R7 (P\*). Entropy decomposition

Let M, M* > 0 lie in L²(𝔽_{t+1}) with E_t|log M| < ∞ and E_t|log M*| < ∞, E_t M = e^{−r_t}, E_t M* = e^{−r*_t} and M S_{t+1}/S_t ∈ L². If one-period markets are complete in L², then M* = M S_{t+1}/S_t, and E_t[Δs + r* − r] = L_t(M) − L_t(M*) and E_t[s_{t+1} − f_t] = L_t(M) − L_t(M*) − x_t. If the cumulant generating function of m has a Maclaurin series with radius of convergence greater than 1, then L_t(M) = Σ_{j≥2} κ_j/j!. The entropy bound E L_t(M) ≥ E[log R − log R_f] for funded returns R > 0 is due to Bansal and Lehmann (1997) and Alvarez and Jermann (2005), as stated in Backus, Chernov and Zin (2014).

Backus, Foresi and Telmer (2001, Proposition 1 and eqs. (9)–(13)) give the decomposition of the currency premium into the difference between the two kernels' log-expectation and expectation-of-log terms, and its cumulant expansion. Backus, Chernov and Zin (2014) supply conditional entropy and the bound, Lustig and Verdelhan (2019) the decomposition with a wedge, and Jiang, Krishnamurthy and Lustig (2021) the convenience yield inferred from the Treasury basis; Brandt, Cochrane and Santa-Clara (2006) is cited for the remark on risk sharing. Writing the covered-interest-parity basis into the identity is an elementary step of my own.

## R8 (P). Local well-posedness of the smile calibration

The calibration solves three equations in (α, ρ, ν) at fixed β: the ATM delta-neutral-straddle volatility, the 25Δ risk reversal at the smile's own 25Δ strikes, and the equality of the market-strangle premium under the smile and under the flat volatility σ_ATM + BF. If the Hagan volatility is smooth, which a lemma establishes for β = 1, the strike fixed points are non-degenerate (R5(c)) and the Jacobian is non-singular at a solution, then the solution is locally unique and C¹ in the quotes, by the implicit function theorem applied first to the strikes and then to the system.

The implicit-function-theorem argument follows Reiswich (2010, Theorems 2 and 3, Appendices A and B), who uses it for the volatility-strike function and for the market-strangle calibration of a simplified parabolic smile. Le Floc'h and Kennedy (2014, Section 4, eqs. (30)–(31)) invert the leading-order expansion of the Hagan formula in closed form, recovering (α, ρ, ν) from the ATM level, skew and curvature, which is an approximate, global counterpart of R8. The three equations are the conditions of Reiswich and Wystup (2012, working-paper Section 2, eqs. (17), (18) and (25)), to which Reiswich (2010, Section 3.2) calibrates SABR, and Henrard (2011) applies the implicit function theorem to any calibration with as many instruments as parameters. In the lemma, the value of z/x(z) at z = 0 is known (Choi and Wu, 2021) and the integral representation is Hadamard's lemma (Milnor, 1963, Lemma 2.1). My own part is the explicit statement that the β = 1 formula is jointly smooth on the stated domain, including ν = 0, and the local well-posedness statement for SABR under conditions (i) to (iv), written for FX delta conventions. The implementation is `src/qef/fx/smile.py`, in the variables (ln α, atanh ρ, ln ν).

## R9 (C). Inference

The central limit theorem for the sample mean holds under strict stationarity, E|X|^{2+δ} < ∞, Σ α(n)^{δ/(2+δ)} < ∞ and a positive long-run variance (Ibragimov, 1962, Theorem 1.7). The stationary bootstrap is that of Politis and Romano (1994), valid under the moment and mixing conditions of Theorem 2 of their technical-report version (Politis and Romano, 1991); Gonçalves and de Jong (2003) weaken the conditions. For the ratio θ, confidence procedures that are bounded with probability one have zero worst-case coverage (Gleser and Hwang, 1987; Dufour, 1997), so test inversion is used; the set is Fieller's (1954, Sections 2–4) for a ratio of means with a Newey–West variance, and it is bounded only when the denominator is significantly different from zero.

## Supporting derivation: the Garman–Kohlhagen PDE (P\*)

The Garman–Kohlhagen price solves its PDE with terminal condition (S − K)⁺, which is verified directly (Garman and Kohlhagen, 1983; Shreve, 2004). It is the unique C^{1,2} solution of polynomial growth that is continuous up to T, by the Feynman–Kac argument with localisation, and no boundary condition at S = 0 is needed, because the exchange rate never reaches zero.
