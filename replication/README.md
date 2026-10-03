# Replication guide

The results of this project can be reproduced at two levels.

| Level | What it reproduces | Access needed |
| --- | --- | --- |
| Public | The results that the [published series](../data/public/fx_carry_portfolio_series/README.md) support: E1 under the market reading with 10-delta hedges, E2 with φ as predictor, the secondary moment predictors, E3 and its split by rate regime, E5, the hedge-cost ratio and its confidence sets, the long at-the-money sample, the post hoc regime comparison, the comparison with published currency portfolios, and the base and at-the-money rows of the robustness grid with the returns and hedge-cost ratio of its 25-delta row | None; the Cboe and Verdelhan files for E5 and the portfolio comparison are public downloads |
| Licensed | The whole pipeline from the quotes: the audit, smile calibration and the series themselves, and the results the public series cannot reproduce, among them E4 and the 25-delta skew prices of E1, every result of the oriented 10-delta risk reversal, the other robustness variants, the three-month tenor and the per-currency analyses | The reader's own LSEG Workspace licence |

The public route is described in the [reproduction instructions](../docs/reproducing_paper.md), which list the results each level reproduces. The licensed route, with its steps, inputs and outputs, is described in [licensed reproduction](licensed_reproduction.md).

Both routes use the project's own code without modification. The theory results and the model validation of the moment code are synthetic and are reproduced by the test suite and `scripts/validate_moments_models.py`.
