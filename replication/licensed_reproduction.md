# Licensed reproduction

A reader with their own LSEG Workspace licence can rebuild the private data layer and rerun every estimation step. The acquisition runs in a separate environment and reads the App Key of whoever runs it from `~/.lseg/app_key`:

```bash
python -m venv .venv-lseg && .venv-lseg/bin/pip install -r requirements-lseg.lock
.venv-lseg/bin/python scripts/acquire_lseg_fx.py --retrieval-date 2026-09-23
python scripts/acquire_cboe_vx.py --retrieval-date 2026-09-24
python scripts/acquire_verdelhan.py --retrieval-date 2026-09-27
```

The three-month forwards and rates, the one-month 5-delta quotes and the pre-1995 coverage check come from further retrievals with the options `--only`, `--start`, `--end` and `--skip-search` of the acquisition script (the further retrievals made no metadata search); the header of `scripts/reproduce.sh` lists them.

## Inputs

`scripts/reproduce.sh` expects, under `data/private/` or in a directory given as its argument: the LSEG retrievals of 23 September 2026 (`raw/`, `manifest.json`, `metadata.csv` and `requests.jsonl`), 24 September 2026 (three-month forwards and rates) and 27 September 2026 (one-month 5-delta quotes) under `lseg/`, with the optional coverage retrieval of 1 October 2026; the Cboe settlement files under `cboe/2026-09-24/raw/`; and Verdelhan's `CurrencyPortfolios.xls` with its manifest under `verdelhan/2026-09-27/`. It checks these inputs before it changes anything.

## Steps

The script builds a new environment from `requirements.lock`, builds the C++ kernel with `cpp/build.sh` (which needs a C++17 compiler) and runs the tests. It then runs every estimation step, in this order:

1. the audit of the 23 September 2026 retrieval (`scripts/audit_fx_panel.py`);
2. the one-month composite and Fenics smile calibrations (`scripts/calibrate_smiles.py`);
3. E1 and E4 (`scripts/estimate_e1.py`);
4. Stage 4: returns, E2 and E3 (`scripts/estimate_stage4.py`);
5. E5 (`scripts/estimate_e5.py`);
6. the moment intervals (`scripts/estimate_moments.py`);
7. the vanna-volga calibration (`scripts/calibrate_vv.py`);
8. the robustness grid (`scripts/robustness.py`);
9. the three-month calibration and estimation (`scripts/calibrate_smiles.py --tenor 3M` and `scripts/estimate_3m.py`);
10. the regime attribution (post hoc), the model validation of the moment code (a design item computed after the main results), and the post hoc sharp identification of option-implied moments and the out-of-sample test with 5-delta quotes (`scripts/estimate_regime_attribution.py`, `scripts/validate_moments_models.py`, `scripts/estimate_identification.py` and `scripts/estimate_wing_test.py`);
11. the comparison with published currency portfolios (`scripts/compare_verdelhan.py`), which reads Verdelhan's file with `xlrd`, pinned in `requirements.lock` (xlrd 2.0.2);
12. the E1 bootstrap intervals recomputed with arch (`scripts/check_bootstrap_arch.py`);
13. the long at-the-money sample (`scripts/estimate_long_atm.py`);
14. the design checks (`scripts/estimate_design_checks.py`);
15. the design items completed on 2 October 2026: E3 by regime, the stale-butterfly check on the extended sample, and the comparison of 1,999 and 9,999 bootstrap draws for the robustness grid, the three-month tenor and the secondary moment predictors (`scripts/estimate_outstanding_items.py`);
16. the descriptive account of smiles, φ and returns around 5 August 2024 (`scripts/describe_august_2024.py`);
17. the checks behind statements of the paper and the theory notes that no earlier step produced: the hypotheses of R5(c) and R8 on the calibrated smiles, a denser grid of means in the identified sets, and the 5-delta breakdown comparison on the same currency-months (`scripts/check_reported_conditions.py`).

It runs no acquisition, and it does not run the kernel benchmark (`scripts/benchmark_kernel.py`), which rewrites the timing block of `cpp/README.md`. The header of `scripts/reproduce.sh` names the report, or the section of the paper, that records each step. It stops at the first step that fails and refuses to overwrite an existing environment, output or link.

## Outputs

Outputs go to `data/private/audit/2026-09-23/`, `data/private/results/2026-09-23/` and, for the three-month estimation, `data/private/results/2026-09-24/`. The calibration logs are written next to the results; the other step logs, the interpreter version and the installed package versions go to `data/private/reproduce/`. Every output is LSEG-derived and stays under `data/private/`.

`scripts/reproduce.sh` does not touch the published series. After it, `scripts/export_public_series.py` rebuilds them from the private results, and `scripts/check_public_reversibility.py` repeats the check that they cannot be traced back to the quotes.

Results can differ slightly from those reported if the provider has revised quotes since my retrievals; the [data audit](../reports/data_audit.md) records a revision check over the days between my two retrievals.
