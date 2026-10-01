# Crash insurance and the G10 carry premium

Does the ex-ante price of crash insurance in G10 FX option smiles explain the collapse of the currency carry premium after 2008 and its revival in 2022–2026?

The data are daily dealer-contributed FX option quotes from LSEG Workspace for nine G10 currencies against the US dollar: at-the-money volatility, 25Δ and 10Δ risk reversals and butterflies, forward points and money-market rates. Each month-end smile is converted into strikes and prices under market quoting conventions. The project then measures the skew price of crash protection per unit of carry, tests whether it explains carry returns in and out of sample, and decomposes the realised cost of hedging into payoff, volatility-level and skew components. The [paper](paper/paper.pdf) and its [extended abstract](paper/summary.md) report the results.

| Document | Content |
| --- | --- |
| [research_design.md](research_design.md) | Question, estimands, samples, inference, stages (pre-registered) |
| [theory/](theory/README.md) | Results R1–R10 and proofs |
| [data_plan.md](data_plan.md) | Instruments, storage and audit checks |
| [research_log.md](research_log.md) | Dated decisions and deviations |
| [references.md](references.md) | Every source, what it is used for and the version consulted |
| [reports/](reports/) | Audit record and results |

## Status

Stages 1 to 5 are complete. Stage 5 comprises the proofs, the full robustness grid (including the three-month tenor and vanna–volga smiles), the R4 moment intervals and the C++ kernel of the pricing, delta, SABR and moment code, with parity tests against the Python reference ([cpp/](../../cpp/README.md)). Results: [E1 and E4](reports/e1.md); [portfolio returns, E2 (with the secondary moment predictors), E3 and E5](reports/stage4.md); [robustness](reports/robustness.md). For Stage 6 the [paper draft](paper/paper.pdf) and [summary](paper/summary.md) are written and reviewed. A clean-environment run of the whole pipeline reproduces every reported number (research log, 26 September 2026). A second retrieval found no quote revisions ([data audit](reports/data_audit.md)). After the pre-registered results, I added three post hoc analyses: the admissible tail exponents and sharp identified sets of option-implied moments (theory result R10), which show how strong a tail assumption the sign of one-month skewness needs; a validation of the moment code against Merton and Heston models ([model validation](reports/model_validation.md)); and a direct comparison of the regime changes in the skew price and the carry spread ([E1 record](reports/e1.md)). An out-of-sample test with one-month 5Δ quotes, planned in the research log before their retrieval, checks the tail hypotheses behind the identification results ([Stage 4 record](reports/stage4.md)). Five items that the design names were computed after the main results and are reported with that label: the long at-the-money sample from 1995, which bounds the hedge-cost ratio but only for an at-the-money hedge ([Stage 4 record](reports/stage4.md)); the SOFR OIS and named-broker checks ([E1 record](reports/e1.md)); the delta-method intervals; and the closed-form checks of the moment code. HML_FX agrees in sign and magnitude with the published portfolios of Lustig, Roussanov and Verdelhan, and the stationary bootstrap agrees with the arch package. The licence holder has confirmed that manipulated or transformed data, which include these pooled results, may be published and that the raw data may not (research log, 1 October 2026).

## Reproduction

```bash
python -m venv .venv && .venv/bin/pip install -r requirements.lock && .venv/bin/pip install -e . --no-deps
.venv/bin/python -m pytest -q
```

Rerunning the LSEG acquisition requires the reader's own LSEG Workspace licence, not mine, and uses a separate environment:

```bash
python -m venv .venv-lseg && .venv-lseg/bin/pip install -r requirements-lseg.lock
.venv-lseg/bin/python scripts/acquire_lseg_fx.py
.venv/bin/python scripts/audit_fx_panel.py
```

The script reads the App Key of whoever runs it from `~/.lseg/app_key`. Given the private data, `scripts/reproduce.sh` builds a new environment from `requirements.lock`, runs the tests and then every estimation step in order; its header states the data layout it expects and the outcome of my own run. I obtained my data under an LSEG Workspace student licence provided by my university, which covers only me, permits individual study and research and does not permit redistribution. The repository therefore holds code, methods and aggregate results (means, standard errors, test statistics and intervals). It holds no LSEG data and no month-level series or calibrated parameters from which quotes could be reconstructed; these stay in `data/private/`, which Git ignores.
