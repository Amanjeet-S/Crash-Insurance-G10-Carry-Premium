# Research roadmap

## 1. Current state

The paper, its theory notes, the estimation record, the robustness grid, the long at-the-money sample, the design checks, the C++ kernel and the public reproduction of the portfolio-level results (apart from E4 and the 25-delta comparisons, whose series are withheld) are complete. The [research log](../research_log.md) records the order of work, including the analyses added after the main results.

## 2. Limitations that further data would address

The paper states its limitations. The hiking regime is short. The share of the carry premium paid for crash insurance at 10 delta needs a significant unhedged premium over a sample of about 24 years of option-hedged returns, while the licensed risk-reversal quotes begin in 2007. Option-implied skewness is identified only under a strong tail assumption, which quotes beyond 10 delta over a longer span would test further. Executable quotes and fixings at the 10:00 New York cut would replace the spread rule and settle payoff timing.

## 3. Later projects

Each later project has its own repository, states a focused question and fixes its design before estimation.
