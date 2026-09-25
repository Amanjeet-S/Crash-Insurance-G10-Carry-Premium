# Crash insurance and the G10 carry premium

Does the ex-ante price of crash insurance in G10 FX option smiles explain the collapse of the currency carry premium after 2008 and its revival in 2022–2026?

The data are daily dealer-contributed FX option quotes from LSEG Workspace, for nine G10 currencies against the US dollar:

- at-the-money volatility;
- 25Δ and 10Δ risk reversals;
- 25Δ and 10Δ butterflies;
- forward points and money-market rates.

Each month-end smile is converted into strikes and prices under market quoting conventions. The project then:

1. measures the skew price of crash protection per unit of carry;
2. tests whether it explains carry returns in and out of sample;
3. decomposes the realised cost of hedging into delta, volatility-level and skew components.

| Document | Content |
| --- | --- |
| [research_design.md](research_design.md) | Question, estimands, samples, inference, stages (pre-registered) |
| [theory/](theory/README.md) | Results R1–R9 and proofs |
| [data_plan.md](data_plan.md) | Instruments, storage and audit checks |
| [research_log.md](research_log.md) | Dated decisions and deviations |
| [references.md](references.md) | Every source, what it is used for and the version consulted |
| [reports/](reports/) | Audit record and results |

**Status:** Stages 1 to 4 complete. Stage 5: proofs, the full robustness grid (including the three-month tenor and vanna–volga smiles) and the R4 moment intervals complete. Results: [E1 and E4](reports/e1.md); [portfolio returns, E2 (with the secondary moment predictors), E3 and E5](reports/stage4.md); [robustness](reports/robustness.md). Stage 6: [paper draft](paper/paper.pdf) and [summary](paper/summary.md) written and reviewed. Remaining: the C++ kernel, the Stage 1 quote-revision check, clean-environment reproduction, and the licence holder's confirmation that pooled summaries may be published.

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

The script reads the App Key of whoever runs it from `~/.lseg/app_key`. I obtained my data under an LSEG Workspace student licence provided by my university, which covers only me, permits individual study and research and does not permit redistribution. The repository therefore holds code, methods and aggregate results (means, standard errors, test statistics and intervals). It holds no LSEG data and no month-level series or calibrated parameters from which quotes could be reconstructed; these stay in `data/private/`, which Git ignores.
