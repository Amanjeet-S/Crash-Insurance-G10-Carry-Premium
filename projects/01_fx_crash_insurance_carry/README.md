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
| [reports/](reports/) | Audit record and results |

**Status:** Stage 1. The LSEG panel is acquired and audited (checks 1 to 5, 8 and 9; see [reports/data_audit.md](reports/data_audit.md)). Checks 6 and 7 and a revision check remain. The pricing and smile engine is implemented and cross-checked against QuantLib. No empirical results yet.

## Reproduction

```bash
python -m venv .venv && .venv/bin/pip install -r requirements.lock && .venv/bin/pip install -e . --no-deps
.venv/bin/python -m pytest -q
```

LSEG acquisition needs a Workspace licence and uses a separate environment:

```bash
python -m venv .venv-lseg && .venv-lseg/bin/pip install -r requirements-lseg.lock
.venv-lseg/bin/python scripts/acquire_lseg_fx.py
.venv/bin/python scripts/audit_fx_panel.py
```

The App Key is read from `~/.lseg/app_key`. LSEG data and every value derived from them stay in `data/private/`, which Git ignores, until the licence terms for publication are confirmed.
