# Data organisation

`public/` holds research data that may be published, each folder with a README stating its origin, definitions and licence, and a manifest of SHA-256 hashes.

- [Portfolio-level series of the crash-insurance project](public/fx_carry_portfolio_series/README.md)

`private/` is ignored by Git. Licensed provider exports, cleaned observations, calibrated parameters, month-level per-currency series, third-party files whose redistribution terms are not established, and all intermediate results belong there. Continuous integration does not need them.

The LSEG data were obtained under a university student licence, which permits individual study and research and does not permit redistribution. The licence holder confirmed on 1 October 2026 that manipulated or transformed data may be published and that the raw data may not; the [publication policy](../docs/publication_policy.md) sets out how a transformed series is judged publishable. Reproduction of the licensed-data analysis requires the reader's own authorised access; the repository does not offer these data on request.

Cleaned data are generated from an identified retrieval by documented code, with units, transformations, date alignment, substitutions and exclusions recorded. A missing observation never silently becomes a zero.

See the [third-party notices](../THIRD_PARTY_NOTICES.md) and the [provenance record](../PROVENANCE.md).
