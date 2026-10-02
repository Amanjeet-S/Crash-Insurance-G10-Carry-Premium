# Third-party material and attribution

The [project references](references.md) identify the mathematical, econometric and empirical methods used. Citations acknowledge prior research and do not imply that this repository owns those methods or that its implementation is an official replication.

## Data read by the code

- **LSEG Workspace:** spot, forward points, FX option quotes and money-market rates, acquired under a university student licence that permits individual study and research and does not permit redistribution. No observation is included. A reader needs their own LSEG Workspace licence to rerun the acquisition.
- **Cboe:** VIX futures daily settlement files from Cboe's public historical data, downloaded by `scripts/acquire_cboe_vx.py`. They are not included, because their redistribution terms are not established; the script retrieves them.
- **Verdelhan's currency portfolios:** the file `CurrencyPortfolios.xls` of Lustig, Roussanov and Verdelhan (2011), from Adrien Verdelhan's data page, downloaded by `scripts/acquire_verdelhan.py`. It is not included, because the page states no terms of use.

## Subscription data and derived series

The licence holder confirmed on 1 October 2026 that manipulated or transformed data may be published and that the raw data may not. Raw observations, cleaned subsets and reversible transformations remain excluded. The portfolio-level series in `data/public/` are original derived results, released after the review described in the [publication policy](docs/publication_policy.md).

## Software and publications

Dependencies retain their own licences and are installed through the lock files; their source is not vendored. Research papers are linked to their academic or publisher sources. Released original software is licensed under [MIT](LICENSE), and released original research writing and series under [CC BY 4.0](LICENSES/CC-BY-4.0.txt), as specified in the [licensing scope](LICENSING.md).
