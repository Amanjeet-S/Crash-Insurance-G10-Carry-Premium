# Research roadmap

## 1. Current state

The paper is a draft. Its theory notes are complete, as are the estimation record of E1 to E5, the robustness grid with the vanna-volga smiles, the three-month tenor, the model validation of the moment code against Merton and Heston models (a design item, computed after the main results), the post hoc analyses (the regime comparison of skew price and carry spread, the identification of option-implied skewness and its out-of-sample test with 5-delta quotes), the comparison with published currency portfolios, the check of the bootstrap against arch, the long at-the-money sample, the design checks and the C++ kernel.

On 2 October 2026 I computed the items that the design or the research log had committed to and that no earlier script had computed: E3 by rate regime, the stale-butterfly check on the extended sample, and the robustness grid, the three-month tenor and the secondary moment predictors with the design's 9,999 bootstrap draws in place of the 1,999 first used (`scripts/estimate_outstanding_items.py`); and a descriptive account of smiles, φ and carry returns around 5 August 2024 (`scripts/describe_august_2024.py`). The paper and the Stage 4 report also compare the hedge-cost ratio θ<sub>UB</sub> with the crash-premium shares reported in the literature. On 3 October 2026 I added the payoff and volatility-level terms and θ₀ to the three-month variant, the test-inversion sets of θ<sub>UB</sub> to the robustness grid and static-arbitrage checks to the daily smiles of the 5 August 2024 account, and ran the full pipeline again from a fresh clone in a new environment.

The public reproduction covers only the results that the published series support; the [reproduction instructions](reproducing_paper.md) list them and the results that need the reader's own LSEG licence. The [research log](../research_log.md) records the order of work, including the analyses added after the main results.

## 2. Limitations that further data would address

The paper states its limitations. The hiking regime is short (56 month-ends), and a longer one would add effective observations. A bounded confidence set for the hedge-cost ratio θ<sub>UB</sub> at 10 delta needs a significant unhedged carry premium, about 24 years of 10-delta-hedged returns at the primary sample's ratio of mean to standard deviation, while the one-month risk-reversal quotes that such a hedge needs begin in 2007 in the licensed panel. Option-implied skewness is identified only under a strong tail assumption, which quotes beyond 10 delta over a longer span would test further. Executable single-option quotes would replace the spread rule, and fixings at the 10:00 New York cut would settle payoff timing.

## 3. Later projects

Each later project has its own repository, states a focused question and fixes its design before estimation.
