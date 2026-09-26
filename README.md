# Quantitative Economics and Finance

Research projects that combine economic theory, financial market data and mathematical methods. Each project states a question, fixes its estimands before estimation, and reports results with their uncertainty and limitations.

| Project | Question | Status |
| --- | --- | --- |
| [01 Crash insurance and the G10 carry premium](projects/01_fx_crash_insurance_carry/README.md) | Does the ex-ante price of crash insurance in FX option smiles explain the post-2008 collapse and 2022–2026 revival of the currency carry premium? | Estimation and robustness complete; paper in preparation |

## Layout

```
projects/<project>/   design, theory, data plan, log, reports
src/qef/              reusable code (FX conventions, pricing, smiles, data handling)
scripts/              acquisition and audit entry points
tests/                analytic and synthetic checks
cpp/                  compiled kernels (planned)
data/public/          attributed public snapshots (none yet)
data/private/         licensed data and derived values (not tracked)
```

## Environment

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.lock && .venv/bin/pip install -e . --no-deps
.venv/bin/python -m pytest -q
```

## Data licence

I obtained the LSEG data used in this research under an LSEG Workspace student licence provided by my university. That licence covers only me: it permits individual study and research and does not permit redistribution, so this repository contains no LSEG data. Anyone else who wants to rerun the acquisition and estimation needs their own LSEG Workspace licence; the data are acquired in a separate environment (`requirements-lseg.lock`).

The repository holds code, methods and aggregate results (means, standard errors, test statistics and intervals). It holds no month-level series or calibrated parameters from which quotes could be reconstructed; these stay in `data/private/`, which Git ignores. Coverage facts (sample windows, counts) are published because they reveal no values. `tests/test_repository_hygiene.py` fails if a data extract or anything under `data/private/` is ever tracked.
