# Reproducing the crash-insurance paper

## 1. What can be checked publicly

The [paper](../paper/README.md) studies whether the price of crash insurance in G10 currency option smiles accounts for the carry premium. Its inputs are LSEG data that I may not redistribute, so the public reproduction starts from the [portfolio-level series](../data/public/fx_carry_portfolio_series/README.md) derived from them.

| Material | Public reproduction | Interpretation |
| --- | --- | --- |
| E1 (market reading, 10 delta), E2, E3, E5, the hedge-cost ratio and its confidence sets, the long at-the-money sample and the regime comparison | Recomputed from the public series by the project's own estimation functions and compared with the numbers printed in the paper | Checks those results; does not recompute the series from quotes |
| VIX roll-down factor (E5) | Rebuilt from Cboe's public settlement files, downloaded by `scripts/acquire_cboe_vx.py` | Optional; the E5 rows are skipped without them |
| Comparison with published currency portfolios | Uses Verdelhan's public file, downloaded by `scripts/acquire_verdelhan.py` (reading it needs `xlrd`) | Optional |
| E4, the 25-delta comparisons of E1, and per-currency results: smile calibration, option-implied moment intervals, identification, the 5-delta test, robustness variants and the three-month tenor | Not reproducible from public data | Requires the reader's own LSEG licence (section 4) |
| Theory results and the model validation of the moment code | Synthetic; reproduced by the tests and `scripts/validate_moments_models.py` | No market data involved |

## 2. Install once

```bash
git clone https://github.com/Amanjeet-S/Crash-Insurance-G10-Carry-Premium.git
cd Crash-Insurance-G10-Carry-Premium
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock
python -m pip install -e . --no-deps
```

## 3. Run the public reproduction

```bash
python scripts/reproduce_from_public.py --out reproduction_report.md
```

The command checks the SHA-256 hashes in the public manifest, recomputes each result, compares it with the paper at the printed precision and lists the results it cannot reproduce with the reason. It exits with a non-zero status on any mismatch. With the optional inputs:

```bash
python scripts/acquire_cboe_vx.py --retrieval-date 2026-09-24
python scripts/acquire_verdelhan.py --retrieval-date 2026-09-27
python scripts/reproduce_from_public.py --out reproduction_report.md \
  --cboe-dir data/private/cboe/2026-09-24/raw \
  --verdelhan-file data/private/verdelhan/2026-09-27/CurrencyPortfolios.xls
```

A reader should check the number of skipped results as well as the exit status.

## 4. Licensed reproduction

With their own LSEG Workspace licence, a reader can rebuild the private layer with `scripts/acquire_lseg_fx.py` and run the whole pipeline with `scripts/reproduce.sh`, which builds a new environment, runs the tests and every estimation step in order. The [replication guide](../replication/README.md) lists the steps, their inputs and outputs.
