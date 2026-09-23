# Quantitative Economics and Finance

Research projects that combine economic theory, financial market data and mathematical methods. Each project states a question, fixes its estimands before estimation, and reports results with their uncertainty and limitations.

| Project | Question | Status |
| --- | --- | --- |
| [01 Crash insurance and the G10 carry premium](projects/01_fx_crash_insurance_carry/README.md) | Does the ex-ante price of crash insurance in FX option smiles explain the post-2008 collapse and 2022–2026 revival of the currency carry premium? | Data audit |

## Layout

```
projects/<project>/   design, theory, data plan, log, reports
src/qef/              reusable code (FX conventions, pricing, smiles, data handling)
scripts/              acquisition and audit entry points
tests/                analytic and synthetic checks
cpp/                  compiled kernels
data/public/          attributed public snapshots
data/private/         licensed data and derived values (not tracked)
```

## Environment

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.lock && .venv/bin/pip install -e . --no-deps
.venv/bin/python -m pytest -q
```

LSEG Workspace data are acquired in a separate environment (`requirements-lseg.lock`) with a user's own licence. Licensed data and values derived from them are kept in `data/private/` until the licence terms for publication are confirmed.
