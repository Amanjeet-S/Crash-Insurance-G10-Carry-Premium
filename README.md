# Crash Insurance and the G10 Carry Premium

## 1. Research question and findings

I ask whether the price of crash insurance in G10 currency option smiles accounts for the carry premium and for its change between the zero-rate regime (May 2013 to December 2021) and the hiking regime from January 2022. Month-end SABR smiles for nine currencies against the US dollar, calibrated to dealer quotes from LSEG Workspace, price one-month protective options on a three-long, three-short carry portfolio. The primary estimand is the skew price of that protection per unit of the portfolio's forward discount; the design also tests whether it predicts carry returns, decomposes the realised hedge cost into payoff, volatility-level and skew terms, and measures exposure to systemic risk.

The rise in the carry spread between the regimes was not matched by a rise in the price of protection, so the regime-shift explanation has partial support only. The share of the carry premium paid for crash insurance is not identified at one month in a sample of this length; over 1995 to 2026 with at-the-money protection the confidence set is bounded but contains the diffusive null. The sign of one-month option-implied skewness is identified only under a strong and explicit assumption about the tails beyond the quoted strikes, which I characterise exactly and test out of sample with 5-delta quotes.

The design was fixed on 23 September 2026 before any return or option-implied statistic was computed. Every later change is recorded in the [research log](research_log.md) with its date and reason, and analyses added after the main results are labelled post hoc in the paper and the reports.

## 2. Reading and execution

Read the [paper (PDF)](paper/Crash_Insurance_and_the_G10_Carry_Premium.pdf), with its [LaTeX source](paper/main.tex) and [extended abstract](paper/summary.md).

| Document | Content |
| --- | --- |
| [research_design.md](research_design.md) | Question, estimands, samples, inference and stages, fixed before estimation |
| [theory/](theory/README.md) | Results R1 to R10 and their proofs |
| [data_plan.md](data_plan.md) | Instruments, storage and audit checks |
| [research_log.md](research_log.md) | Dated decisions and deviations |
| [references.md](references.md) | Every source, what it is used for and the version consulted |
| [reports/](reports/) | Data audit, smile calibration, E1 and E4, returns and E2 to E5, robustness and model validation |
| [replication/](replication/README.md) | Public and licensed reproduction |

Python implements the data audit, smile calibration, estimation and inference; a [C++ kernel](cpp/README.md) of the pricing, delta, SABR and moment code is checked against the Python reference, and QuantLib, arch and statsmodels serve as independent cross-checks.

The [reproduction instructions](docs/reproducing_paper.md) give one command that recomputes the paper's portfolio-level results from the [published series](data/public/fx_carry_portfolio_series/README.md) and compares each with the printed number; it needs no data licence. E4 and the 25-delta comparisons are the exception, because their series are withheld. Recomputing the series and the per-currency results needs the reader's own LSEG Workspace licence. The [development guide](docs/development.md) covers installation and the automated checks.

## 3. Repository structure

```text
Crash-Insurance-G10-Carry-Premium/
├── .github/workflows/   Automated tests
├── cpp/                 C++ kernel with parity tests
├── data/public/         Published portfolio-level series
├── docs/                Development, reproduction and publication policies
├── paper/               Paper, LaTeX source and extended abstract
├── replication/         Public and licensed reproduction guides
├── reports/             Audit, calibration and estimation records
├── scripts/             Acquisition, estimation and verification entry points
├── src/qef/             Reusable Python implementation
├── tests/               Analytic, numerical and data-hygiene checks
└── theory/              Theory notes with proofs
```

## 4. Evidence and licensing

I distinguish observed market data, derived portfolio-level series, synthetic validation models and stated mathematical results. They answer different questions and are labelled accordingly. Established methods are cited where they are used; claims of originality are made only after a literature check.

The LSEG data were obtained under a university student licence that permits individual study and research and does not permit redistribution. The licence holder confirmed on 1 October 2026 that manipulated or transformed data may be published and that the raw data may not. Released results are research findings, original code and writing, and portfolio-level series that cannot be traced back to the quotes. Raw provider data, per-currency series and calibrated parameters remain private. The [publication policy](docs/publication_policy.md) specifies the release boundary.

Released original software uses the [MIT licence](LICENSE), and released original research writing and series use [CC BY 4.0](LICENSES/CC-BY-4.0.txt). The [licensing scope](LICENSING.md), [provenance record](PROVENANCE.md) and [third-party notices](THIRD_PARTY_NOTICES.md) preserve external rights and exclude restricted observations. To cite this work, use the [citation file](CITATION.cff).
