# Replication guide

The results of this project can be reproduced at two levels.

| Level | What it reproduces | Access needed |
| --- | --- | --- |
| Public | The portfolio-level results of the paper (E1 at 10 delta, E2, E3, E5, the hedge-cost ratio and its confidence sets, the long at-the-money sample, the regime comparison and the design checks that use portfolio series), recomputed from the [published series](../data/public/fx_carry_portfolio_series/README.md) | None; the Cboe and Verdelhan files for E5 and the portfolio comparison are public downloads |
| Licensed | The whole pipeline from the quotes, including E4 and the 25-delta comparisons: audit, smile calibration, the series themselves, the per-currency moment and identification results, the robustness grid and the three-month tenor | The reader's own LSEG Workspace licence |

The public route is described in the [reproduction instructions](../docs/reproducing_paper.md). The licensed route is described in [licensed reproduction](licensed_reproduction.md).

Both routes use the project's own code without modification. The theory results and the model validation of the moment code are synthetic and are reproduced by the test suite and `scripts/validate_moments_models.py`.
