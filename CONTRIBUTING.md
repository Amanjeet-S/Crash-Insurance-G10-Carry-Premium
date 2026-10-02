# Development and contributions

Read the [development guide](docs/development.md), [research roadmap](docs/roadmap.md) and [source record](PROVENANCE.md) before proposing a change.

Explain the research question or defect addressed, state the assumptions affected, and include the smallest reproducible example that supports the change. Distinguish implementation corrections from changes to the research design: the design of each project is fixed before estimation, and a later change is recorded in the project's research log with its date and reason. Update documentation when conventions or results change.

Run the [public reproduction](docs/reproducing_paper.md) and the automated tests before submitting a change. Keep empirical observations separate from synthetic validation inputs. Do not select specifications solely because they produce a preferred sign or significance level. Record unsuccessful checks and material limitations. Credit ideas, datasets and adapted code at the point of use.

Keep credentials, account exports, private correspondence and licensed provider files outside the tracked repository. `tests/test_repository_hygiene.py` fails if a data extract or anything under `data/private/` is tracked.

## Licensing contributions

Original software contributions are offered under [MIT](LICENSE), and original research-writing contributions under [CC BY 4.0](LICENSES/CC-BY-4.0.txt), unless separately agreed and clearly documented. Contributors retain their copyright. Follow the [licensing scope](LICENSING.md), identify third-party material and its terms, and do not submit restricted data.
