# Source and provenance record

This repository distinguishes established methods from its implementation, its empirical results and its own mathematical contributions. Citations appear close to the relevant method or claim, and each project keeps a fuller source register.

## Current release

| Material | Origin | Attribution record |
| --- | --- | --- |
| Asset-pricing, option-pricing and econometric methods | Established research literature | [Project references](references.md), with the version of each source consulted |
| Theory results R1 to R10, with proofs of all but R9, which is cited | Known results credited where they are used, and derivations prepared for this project | [Theory notes](theory/README.md) |
| Pricing, calibration, estimation and inference code, including the C++ kernel | Prepared for this repository | Source files, tests and [kernel record](cpp/README.md) |
| Numerical linear algebra, optimisation and random generation | Installed third-party libraries | Dependency locks and [third-party notices](THIRD_PARTY_NOTICES.md) |
| Independent cross-checks | QuantLib (pricing, deltas, strikes and Hagan volatilities), the R package sandwich (Newey-West standard errors, checked once) and the arch package (stationary bootstrap and block length), used as separate implementations | [Validation section of the paper](paper/main.tex) (Section 4.6) and the cross-check tests |
| FX spot, forwards, option quotes and rates | LSEG Workspace, under a university student licence | Acquisition code and the [data audit](reports/data_audit.md); observations held privately |
| VIX futures settlement prices | Cboe public historical data | Acquisition code; files held privately because their redistribution terms are not established |
| Currency portfolios of Lustig, Roussanov and Verdelhan | Adrien Verdelhan's public data page | Acquisition code with URL, retrieval time and hash; file held privately because no terms are stated |
| Portfolio-level month series | Derived in this project from the private data | [Public series and their definitions](data/public/fx_carry_portfolio_series/README.md) |
| Model validation of the moment code | Synthetic Merton and Heston models | [Model validation record](reports/model_validation.md) |

No external code file, paper table or figure is incorporated into the implementation. The [public reproduction](docs/reproducing_paper.md) recomputes the paper's results that the public series support, using the public Cboe and Verdelhan files for E5 and the comparison with published portfolios, and lists the results it cannot reproduce with the reason; it does not read or export LSEG observations.

## Rules for additions

- Write explanations from an understood mathematical argument and cite the source of the underlying result.
- Identify quotations as quotations, with a precise source. Changing a few words is not a substitute for attribution.
- Record the origin, licence and modifications of any reused code before adding it.
- Record the origin, transformations and redistribution conditions of any external dataset.
- Generate figures from documented code, or label and attribute any permitted external figure.
- Reserve claims of originality for results checked against the literature, and record post hoc analyses as post hoc.
- Keep numerical assumptions, unsuccessful checks and material limitations visible.

This record documents the material used and the attribution practice. It is not an automated similarity certificate.
