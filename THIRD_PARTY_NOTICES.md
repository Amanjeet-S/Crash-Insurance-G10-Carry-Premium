# Third-party material and attribution

The [project references](references.md) identify the mathematical, econometric and empirical methods used. Citations acknowledge prior research and do not imply that this repository owns those methods or that its implementation is an official replication.

## Data read by the code

- **LSEG Workspace:** spot, forward points, FX option quotes and money-market rates, acquired under a university student licence that permits individual study and research and does not permit redistribution. No observation is included. A reader needs their own LSEG Workspace licence to rerun the acquisition.
- **Cboe:** VIX futures daily settlement files from Cboe's public historical data, downloaded by `scripts/acquire_cboe_vx.py`. They are not included, because their redistribution terms are not established; the script retrieves them.
- **Verdelhan's currency portfolios:** the file `CurrencyPortfolios.xls` of Lustig, Roussanov and Verdelhan (2011), from Adrien Verdelhan's data page, downloaded by `scripts/acquire_verdelhan.py`. It is not included, because the page states no terms of use.

## Subscription data and derived series

The licence holder confirmed on 1 October 2026 that manipulated or transformed data may be published and that the raw data may not. Raw observations, cleaned subsets and reversible transformations remain excluded. The portfolio-level series in `data/public/` are original derived results, released after the review described in the [publication policy](docs/publication_policy.md).

## Code that follows third-party implementations

Two routines in this repository were written to follow an existing implementation step by step, so that their results agree with it exactly. They are not original MIT material in that respect, and the notices of the two projects apply to them:

- **NumPy:** `pairwise_sum` in `cpp/qef_kernel.cpp` follows NumPy's pairwise summation routine (`numpy/_core/src/umath/loops_utils.h.src`), so that the kernel's sums agree bitwise with `numpy.sum`. NumPy is distributed under the BSD 3-Clause licence, Copyright (c) 2005-2025, NumPy Developers; the text is in [LICENSES/NumPy-BSD-3-Clause.txt](LICENSES/NumPy-BSD-3-Clause.txt).
- **arch:** `_block_length` in `src/qef/stats/bootstrap.py` follows the steps of `optimal_block_length` in the arch package (Kevin Sheppard), including its caps, with two stated differences of convention and one degenerate case. arch is distributed under the University of Illinois/NCSA Open Source Licence, Copyright (c) 2017 Kevin Sheppard; the text is in [LICENSES/arch-NCSA.txt](LICENSES/arch-NCSA.txt).

Both licence texts are copied unchanged from the installed distributions (numpy 2.3.5 and arch 8.0.0).

## Software and publications

Dependencies retain their own licences and are installed through the lock files; their source is not vendored. Research papers are linked to their academic or publisher sources. Released original software is licensed under [MIT](LICENSE), and released original research writing and series under [CC BY 4.0](LICENSES/CC-BY-4.0.txt), as specified in the [licensing scope](LICENSING.md).
