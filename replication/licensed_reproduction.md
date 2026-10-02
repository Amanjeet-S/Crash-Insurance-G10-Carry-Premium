# Licensed reproduction

A reader with their own LSEG Workspace licence can rebuild the private data layer and rerun every step. The acquisition runs in a separate environment and reads the App Key of whoever runs it from `~/.lseg/app_key`:

```bash
python -m venv .venv-lseg && .venv-lseg/bin/pip install -r requirements-lseg.lock
.venv-lseg/bin/python scripts/acquire_lseg_fx.py --retrieval-date 2026-09-23
python scripts/acquire_cboe_vx.py --retrieval-date 2026-09-24
python scripts/acquire_verdelhan.py --retrieval-date 2026-09-27
```

The three-month forwards and rates, the one-month 5-delta quotes and the pre-1995 coverage check come from further retrievals with the options `--only`, `--start` and `--end` of the acquisition script; the header of `scripts/reproduce.sh` lists them. Given the private inputs, `scripts/reproduce.sh` builds a new environment from `requirements.lock`, builds the C++ kernel, runs the tests and then every step in order: the audit, the one-month composite and Fenics calibrations, E1 and E4, Stage 4 (returns, E2 and E3), E5, the moment intervals, the vanna-volga calibration, the robustness grid, the three-month tenor, the post hoc analyses, the long at-the-money sample and the design checks. `scripts/export_public_series.py` then rebuilds the published series, and `scripts/check_public_reversibility.py` repeats the check that they cannot be traced back to the quotes.

Results can differ slightly from those reported if the provider has revised quotes since my retrievals; the [data audit](../reports/data_audit.md) records a revision check over the days between my two retrievals.
