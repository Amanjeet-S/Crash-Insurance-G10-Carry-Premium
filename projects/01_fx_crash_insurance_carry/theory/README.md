# Theory

Results used in the [research design](../research_design.md). Proofs are written in `notes.tex`.

**Labels:**

- **P:** proved here.
- **P\*:** a standard result, proved here for completeness.
- **C:** cited, with hypotheses stated.

## R1 (P\*). Domestic and foreign risk-neutral measures

**Setting.**
- (Ω, 𝔽, (𝔽_t), Q^d) satisfies the usual conditions.
- B^d and B^f are strictly positive money-market accounts with B^d_0 = B^f_0 = 1.
- S is strictly positive and càdlàg.
- Z_t = S_t B^f_t / (S_0 B^d_t) is a true Q^d-martingale.
- Define dQ^f/dQ^d = Z_T.

**Claims.**
- (a) Q^f ~ Q^d.
- (b) E^{Q^f}[X | 𝔽_t] = E^{Q^d}[Z_T X | 𝔽_t] / Z_t, for X ≥ 0 or X ∈ L¹(Q^f).
- (c) Foreign prices discounted by B^f are Q^f-martingales exactly when their domestic values discounted by B^d are Q^d-martingales.
- (d) Under Garman–Kohlhagen, E^{Q^d} S_T = F and E^{Q^f}(1/S_T) = 1/F. This is Siegel's paradox resolved.

**Remark.** If Z is a strict local martingale, then E Z_T < 1 and the symmetry fails (Carr, Fisher and Ruf, 2014).

**Source.** Geman, El Karoui and Rochet (1995).

## R2 (P). Arbitrage-free call prices

**Statement.** c(K) = E(X − K)⁺ for some X ≥ 0 with E X = F if and only if:
- c is convex,
- c(0) = F,
- c′₊(0) ≥ −1, and
- c(K) → 0 as K → ∞.

**Consequences.**
- The law of X is unique, with distribution function G = 1 + c′₊ and P(X = 0) = 1 + c′₊(0).
- With c(K) = F − K for K < 0, the law equals c″ in D′(ℝ).

**Tests used.**
- On a grid: slopes in [−1, 0] and nondecreasing.
- On total variance w(k): g(k) ≥ 0 (Gatheral and Jacquier, 2014).

## R3 (P). Spanning

**Hypotheses.**
- f′ is locally absolutely continuous on (0, ∞).
- Q(X > 0) = 1.
- ∫₀^F |f″|P dK + ∫_F^∞ |f″|C dK < ∞.

**Statement.** E f(X) = f(F) + ∫₀^F f″(K)P(K) dK + ∫_F^∞ f″(K)C(K) dK.

**Proof.** Taylor's theorem with integral remainder, then Fubini.

**Sources.** Carr and Madan (1998); Bakshi, Kapadia and Madan (2003).

## R4 (P). Partial identification from truncated quotes

**Tail bound.**
- If Q(S_T ≤ K)/K^α is nondecreasing on (0, K_min], then P(K) ≤ P(K_min)(K/K_min)^{α+1} there.
- Hence ∫₀^{K_min} P K⁻² dK ≤ P(K_min)/(α K_min).

**Further bounds.** Analogous bounds for the Bakshi–Kapadia–Madan weights, and for the upper tail.

**Output.** Variance and skewness are reported as intervals over a stated range of α, truncating at the 10Δ strikes.

**Quadrature.** F is a node, because the integrand has a kink there. Trapezoid error is O(h²) and Simpson error O(h⁴) on each side.

## R5. Delta-to-strike inversion

**(a) Closed form.** For pips spot delta: K = F exp(σ²τ/2 − φσ√τ Φ⁻¹(|Δ|/D_b)).

**(b) (P) Premium-adjusted call delta** h(K) = (K/F) D_b Φ(d₋):
- h → 0 at both ends.
- h′ has the sign of Φ(d₋) − φ(d₋)/(σ√τ).
- Since φ/Φ is strictly decreasing, h is strictly increasing then strictly decreasing, with a unique maximiser K*.
- Each attainable Δ has exactly two strikes; the convention takes the root in [K*, ∞).

**(c)** Conditions on σ(·) under which K ↦ Δ(K, σ(K)) is monotone; checked numerically on the panel.

**Implementation.** `src/qef/fx/gk.py`.

## R6 (P). Diffusive null of the hedge cost

**Setting.** Garman–Kohlhagen with equal P and Q volatility and no jumps; the P-drift of log S exceeds the Q-drift by λ.

**Statement.**
- E^P(K − S_T)⁺ − E^Q(K − S_T)⁺ = −FλτΦ(−d₁) + o(λτ).
- Hence θ₀ = Φ(−d₁) + O(λτ).

**Proof.** Dominated convergence, together with ∂_F E(K − S_T)⁺ = −Φ(−d₁).

## R7 (P). Entropy decomposition

**Hypotheses.**
- M, M* > 0 in L²(𝔽_{t+1}).
- E_t|log M| < ∞ and E_t|log M*| < ∞.
- E_t M = e^{−r_t} and E_t M* = e^{−r*_t}.
- M S_{t+1}/S_t ∈ L².

**Lemma.** If one-period markets are complete in L², then M* = M S_{t+1}/S_t.

**Claims.**
- (i) E_t[Δs + r* − r] = L_t(M) − L_t(M*), and E_t[s_{t+1} − f_t] = L_t(M) − L_t(M*) − x_t.
- (ii) If the cumulant generating function of m has a Maclaurin series with radius of convergence greater than 1, then L_t(M) = Σ_{j≥2} κ_j/j!.

**Entropy bound** (Bansal and Lehmann, 1997). E L_t(M) ≥ E[log R − log R_f] for funded returns R > 0.

**Sources.** Backus, Foresi and Telmer (2001); Backus, Chernov and Zin (2014); Brandt, Cochrane and Santa-Clara (2006).

## R8 (P). Local well-posedness of the smile calibration

**The system.** Three equations in (α, ρ, ν) at fixed β:
- the ATM delta-neutral-straddle volatility;
- the 25Δ risk reversal at the smile's own 25Δ strikes;
- equality of the market-strangle premium under the smile and under the flat volatility σ_ATM + BF.

**Statement.** If the Jacobian is non-singular at a solution and R5(c) holds at the three strikes, the solution is locally unique and C¹ in the quotes (implicit function theorem).

**Implementation.** `src/qef/fx/smile.py`, in the variables (ln α, atanh ρ, ln ν).

## R9 (C). Inference

**Central limit theorem.**
- Strict stationarity.
- E|X|^{2+δ} < ∞.
- Σ α(n)^{δ/(2+δ)} < ∞.
- A positive long-run variance.

**Stationary bootstrap.** Politis and Romano (1994, Thm 1), or Gonçalves and de Jong (2003).

**The ratio θ.** Confidence procedures that are bounded with probability one have zero worst-case coverage (Gleser and Hwang, 1987; Dufour, 1997), so test inversion is used.

## Supporting derivation: the Garman–Kohlhagen PDE

- Via Feynman–Kac, with terminal condition (S − K)⁺.
- Unique among C^{1,2} solutions of polynomial growth.
- No boundary condition at S = 0.
