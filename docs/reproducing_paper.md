# Reproducing the crash-insurance paper

## 1. What can be checked publicly

The [paper](../paper/README.md) studies whether the price of crash insurance in G10 currency option smiles accounts for the carry premium. Its inputs are LSEG data that I may not redistribute, so the public reproduction starts from the [portfolio-level series](../data/public/fx_carry_portfolio_series/README.md) derived from them.

| Material | Public reproduction | Interpretation |
| --- | --- | --- |
| E1 under the market reading with 10-delta hedges; E2 with φ as predictor; the secondary moment predictors of E2; E3 and its split by rate regime; the hedge-cost ratio θ<sub>UB</sub> and its confidence sets, including the delta-method intervals; the long at-the-money sample; the post hoc regime comparison of the skew price and the carry spread; and, in the robustness grid, the base and at-the-money rows and the returns and θ<sub>UB</sub> of the 25-delta row | Recomputed from the public series by the project's own estimation functions and compared with the numbers printed in the paper | Checks those results; does not recompute the series from quotes |
| E5, with its VIX roll-down factor | Factor rebuilt from Cboe's public settlement files, downloaded by `scripts/acquire_cboe_vx.py` | Optional; the E5 rows are skipped without them |
| Comparison with Verdelhan's published currency portfolios | Uses Verdelhan's public file, downloaded by `scripts/acquire_verdelhan.py`; reading it needs `xlrd`, which `requirements.lock` pins | Optional; skipped without the file |
| E4 and the 25-delta skew prices of E1; every result of the oriented 10-delta risk reversal (the secondary E2 predictor and its regime change); the other robustness variants, vanna-volga smiles included; the three-month tenor; the per-currency analyses (smile calibration and its diagnostics, per-currency moment intervals, identification of option-implied skewness, the 5-delta test, the splice of Fenics and composite quotes, quote-level audit counts and coverage, and the design checks that recalibrate smiles or test the moment code); the check of the Cboe settlements against LSEG's VIX futures continuation; the stale-butterfly check on the extended sample; the descriptive account around 5 August 2024; and two details of the long sample (the at-the-money payoff and premium in 2008, and θ<sub>0</sub> of the at-the-money hedge) | Not reproducible from public data; `scripts/reproduce_from_public.py` lists each with its reason | Requires the reader's own LSEG licence (section 4) |
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

With their own LSEG Workspace licence, a reader can rebuild the private layer with `scripts/acquire_lseg_fx.py` and run the pipeline with `scripts/reproduce.sh`, which builds a new environment from `requirements.lock`, runs the tests and then every step of the pipeline in order, except the data acquisitions (LSEG, Cboe and Verdelhan) and the kernel benchmark (`scripts/benchmark_kernel.py`, which rewrites the timing block of `cpp/README.md`). The comparison with published portfolios reads Verdelhan's file with `xlrd`, which `requirements.lock` pins (xlrd 2.0.2). [Licensed reproduction](../replication/licensed_reproduction.md) lists the steps, their inputs and outputs.
