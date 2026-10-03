# Development and reproduction

For the paper's public material, start with the [reproduction instructions](reproducing_paper.md). The commands below cover development.

## Installation

Clone the repository and run these commands from its root. Python 3.13 is the reference version; continuous integration also checks Python 3.11.

```bash
git clone https://github.com/Amanjeet-S/Crash-Insurance-G10-Carry-Premium.git
cd Crash-Insurance-G10-Carry-Premium
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock
python -m pip install -e . --no-deps
python -m pytest -q
```

On Windows, activate with `.venv\Scripts\Activate.ps1` in PowerShell. The LSEG acquisition uses a separate environment built from `requirements-lseg.lock` and needs the reader's own LSEG Workspace licence; nothing else requires it.

## C++ kernel

The [C++ kernel](../cpp/README.md) of the pricing, delta, SABR and moment code needs a C++17 compiler. `bash cpp/build.sh` builds it; its parity tests build it when needed and skip, with the reason, if no compiler is installed. `scripts/benchmark_kernel.py` rewrites the timing section of `cpp/README.md`.

## Cross-checks

`requirements.lock` pins the two libraries that serve as independent implementations: QuantLib, against which `tests/test_quantlib_crosscheck.py` checks premia, deltas, strike inversion, delta-neutral strikes and Hagan volatilities, and arch, against which `tests/test_bootstrap_arch.py` and `scripts/check_bootstrap_arch.py` check the stationary bootstrap and its block length. Both tests skip, with the reason, if their library is missing; the model validation of the moment code (`tests/test_moments_models.py` and `scripts/validate_moments_models.py`) needs QuantLib. statsmodels and patsy are pinned only because arch depends on them; no project code or test imports statsmodels. The Newey-West standard errors were checked once against the R package sandwich, which the project does not install. mpmath is not in `requirements.lock`: the 50-digit check of the inverse normal in `tests/test_cpp_kernel.py` skips without it, and `scripts/benchmark_kernel.py` reports that check only when it is installed (`pip install -e .[crosscheck]`). matplotlib and its dependencies are in the lock file but no project code imports them.

## Continuous integration

GitHub Actions installs the frozen dependencies on Python 3.11 and 3.13 and runs the test suite, including the parity tests of the C++ kernel, the checks of the public portfolio series and `tests/test_repository_hygiene.py`, which fails if a data extract or anything under `data/private/` is tracked.
