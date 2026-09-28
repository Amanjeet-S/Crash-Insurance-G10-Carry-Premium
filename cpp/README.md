# C++ kernel

This directory holds the C++ kernel of Stage 5 of the research design for
project 01 (crash insurance and the G10 carry premium): the inner loops of the
pricing, delta, SABR and moment code, with parity against the Python reference
in `src/qef/fx`. The Python reference remains the code of record; the kernel
is a faster implementation of the same arithmetic, and the parity tests show
that the two agree to the stated tolerances.

| File | Contents |
|---|---|
| `qef_kernel.cpp` | the kernel, with the formulas, the sources and the numerical choices in its header comment |
| `qef_kernel.h` | the C interface (`extern "C"`), its array conventions and status codes |
| `build.sh` | builds `build/libqef_kernel.dylib` (macOS) or `build/libqef_kernel.so` (elsewhere) |
| `../src/qef/kernel.py` | ctypes wrappers with the signatures of the Python functions they mirror |
| `../tests/test_cpp_kernel.py` | parity tests; builds the library if needed and skips without a compiler |
| `../scripts/benchmark_kernel.py` | timings and parity on realistic synthetic workloads; writes the section below |

## Build and use

The kernel is plain C++17 with no dependencies beyond the standard library,
compiled by the system compiler (`$CXX`, default `c++`):

```
bash cpp/build.sh
```

The build directory is not tracked (`.gitignore`). From Python:

```python
from qef import kernel
from qef.fx.gk import CALL, DeltaConvention
from qef.fx.smile import SabrSmile

kernel.forward_premium(1.10, [1.05, 1.15], 0.08, 1 / 12, CALL)
kernel.strike_from_delta(0.25, 150.0, 0.10, 1 / 12, CALL, DeltaConvention(spot=True, premium_adjusted=True))
smile = SabrSmile(1.10, 1 / 12, 0.08, -0.2, 1.0)
kernel.middle_contracts(1.10, 1.05, 1.16, 1 / 12, smile.vol, False)  # the dict of qef.fx.moments.middle_contracts
kernel.middle_contracts_sabr(F, K_min, K_max, tau, alpha, rho, nu, usd_base)  # arrays, one element per smile
```

If the library is absent, importing `qef.kernel` still works and the first
call raises `KernelUnavailableError` with the build command. The tests and
the benchmark run with:

```
.venv/bin/python -m pytest tests/test_cpp_kernel.py
.venv/bin/python scripts/benchmark_kernel.py
```

## Design

Each kernel routine mirrors one Python function operation by operation: the
same expressions, in the same order, with every product and sum rounded
separately (`-ffp-contract=off`, never `-ffast-math`). Nodes are placed as
`numpy.linspace` places them, and Simpson's rule sums with NumPy's partial
pairwise summation, so that the two implementations differ only where a
library function differs. On this machine NumPy's `log`, `exp` and `sqrt`
agree bitwise with the C library, and the SABR volatilities agree bitwise with
the reference.

The one material difference is the normal distribution function. The kernel
uses Φ(x) = erfc(−x/√2)/2 with the C library's `erfc`; SciPy's `ndtr` forms
the same argument but uses its own `erfc`. Measured on this machine, the two
differ by about 20 units in the last place for \|x\| < 6 (19 on 2,000,001
evenly spaced points of [−6, 6], 21 on 5 million uniform draws; sample
maxima, not a bound), most for x < −4 and by at most 2 for x ≥ 0. In the
lower tail neither is accurate to better than tens of units in the last
place, because the condition number of Φ grows like x²: at the 200 points of
largest disagreement the kernel's error against 40-digit values is at most 54
units and `ndtr`'s at most 65. Where a premium is a difference of two nearly equal terms (out of
the money, or at the money with a small σ√τ) that difference is amplified by
the condition number of the subtraction, so the premium and the cubic moment
contract are compared normwise, relative to the size of the terms, as
explained in `tests/test_cpp_kernel.py`. At the points of largest
disagreement, a 50-digit evaluation of the same formula at the same d± (which
both implementations form identically) shows how much of the disagreement is
each implementation's rounding; the benchmark section reports it.

The inverse normal is computed by Newton's method on ln Φ, which converges
monotonically from an explicit starting point (the proof is in the header of
`qef_kernel.cpp`), and the premium-adjusted strikes by Newton's method
safeguarded by bisection on the reference's own bracket, with the reference's
checks. On one sample, the first 150 premium-adjusted targets of each
premium-adjusted convention that `test_strike_from_delta_parity` retains (300
in all), the kernel's strikes are at most 4.6 units in the last place
(`np.spacing` of the root) from 40-digit roots, with median 0.6, and the
reference's at most 36, with median 0.5. These are sample maxima, not
bounds; both lie well inside the reference's guarantee from
`scipy.optimize.brentq`, 1e-14 F plus 1e-14 relative, and at 108 of the 300
points the kernel's root is the further from the exact one. In the adaptive middle part the kernel reuses the integrand at the
previous level's nodes, which are bitwise the even nodes of the next level.

Differences from the reference, all outside its documented use, are listed in
the docstring of `src/qef/kernel.py`: array inputs where the reference is
scalar, β = 1 only for SABR, SABR smiles only for the middle part, and the
reference's `ZeroDivisionError` for a zero passed as a Python float.

## Benchmark

<!-- benchmark:begin -->

Generated by `scripts/benchmark_kernel.py` on 28 September 2026; do not edit by hand.

### Build

Compiler: Apple clang version 21.0.0 (clang-2100.1.1.101). Command, as run by `cpp/build.sh`:

```
c++ -O2 -std=c++17 -fPIC -shared -ffp-contract=off -Wall -Wextra -o cpp/build/libqef_kernel.dylib cpp/qef_kernel.cpp
```

### Coverage

| Kernel wrapper (`qef.kernel`) | Python reference | Branches and conventions covered | Parity tests |
|---|---|---|---|
| `forward_premium` | `qef.fx.gk.forward_premium` | calls and puts; any broadcast shape | `test_forward_premium_*` |
| `delta` | `qef.fx.gk.delta` | spot and forward, pips and premium-adjusted; the nine G10 conventions | `test_delta_*` |
| `sabr_vol` | `qef.fx.sabr.sabr_vol`, β = 1 | both z/x(z) branches and the branch point, ν = 0, ρ near ±1, deep out-of-the-money strikes, NaN and infinite inputs | `test_sabr_vol_*` |
| `pa_call_delta_maximiser` | `qef.fx.gk.pa_call_delta_maximiser` | bracket expansion and root | `test_strike_near_maximal_*` |
| `strike_from_delta` | `qef.fx.gk.strike_from_delta` | pips closed form; premium-adjusted put root; premium-adjusted call root on [K*, ∞); every check and error of the reference | `test_strike_from_delta_*`, `test_strike_near_maximal_*` |
| `middle_contracts`, `middle_contracts_fixed`, `middle_contracts_sabr` | `qef.fx.moments.middle_contracts` for a SABR smile | k = 1, 2, 3; EURUSD-type and USD-base pairs; fixed n and the adaptive doubling; many smiles in one call | `test_middle_contracts_*` |

### Parity tolerances

From `qef.kernel.PARITY_TOLERANCES`; `tests/test_cpp_kernel.py` states why each looser bound holds.

| Quantity | Criterion |
|---|---|
| normal distribution function, deltas, SABR volatilities, pips strikes | relative error ≤ 1e-13 |
| inverse normal | absolute error ≤ 1e-13 max(1, \|x\|) |
| premium | \|ΔC\| ≤ 1e-13 [F Φ(φd+) + K Φ(φd−)]; relative ≤ 1e-13 where that sum over \|C\| is at most 10 |
| premium-adjusted strikes | \|ΔK\| ≤ 1e-14 F + 1e-14 \|K\| + 1e-15 \|K\| (the reference's brentq tolerance plus 4.5 ulp); near the maximal call delta, the reference's delta at the kernel's strike within 1e-14 of the target |
| middle part | \|ΔI_k\| ≤ 1e-13 A_k, A_k Simpson's rule for \|h_k'' P\|; relative ≤ 1e-13 where A_k/\|I_k\| ≤ 10 (always for k = 1, 2) |

### Timings

Machine: Apple M2, Darwin 25.6.0; Python 3.13.9, NumPy 2.3.5, SciPy 1.16.3; one thread. Median of 5 repetitions (3 for the reference's Python loops) after a warm-up call. In the adaptive middle part the reference converged at n from 256 to 4096 subintervals per side (median 1024).

| Kernel | Workload | Python reference | C++ kernel | Speed-up |
|---|---|---:|---:|---:|
| `forward_premium` | 1,000,000 options (both vectorised) | 61.6 ms | 41.9 ms | 1.5× |
| `delta` (spot, pips) | 1,000,000 options (both vectorised) | 46.4 ms | 33.6 ms | 1.4× |
| `delta` (spot, premium-adjusted) | 1,000,000 options (both vectorised) | 44.4 ms | 34.4 ms | 1.3× |
| `delta` (forward, pips) | 1,000,000 options (both vectorised) | 59.6 ms | 50.6 ms | 1.2× |
| `delta` (forward, premium-adjusted) | 1,000,000 options (both vectorised) | 60.0 ms | 34.5 ms | 1.7× |
| `sabr_vol` | one smile, 1,000,000 strikes (both vectorised) | 52.2 ms | 19.0 ms | 2.7× |
| `sabr_vol` | 1,440 smiles × 513 strikes (reference loops over smiles) | 45.9 ms | 14.8 ms | 3.1× |
| `strike_from_delta` (pips) | 1,280 strikes (reference loops) | 40.8 ms | 586 µs | 69.6× |
| `strike_from_delta` (premium-adjusted) | 1,600 strikes (reference loops) | 1,199.2 ms | 3.4 ms | 351.3× |
| middle part, adaptive | 1,440 SABR smiles, one kernel call | 2,165.9 ms | 232.4 ms | 9.3× |
| middle part, adaptive | 1,440 SABR smiles, one kernel call per smile | 2,165.9 ms | 253.6 ms | 8.5× |
| middle part, fixed n | 1,440 SABR smiles at n = 1024 (reference integrate(n) loop) | 557.3 ms | 187.5 ms | 3.0× |

### Parity on the benchmark workloads

| Kernel | Results | Largest relative error | Other criteria | Bitwise equal |
|---|---:|---|---|---:|
| `forward_premium` | 1,000,000 | 3.8e-13 | normwise 1.2e-15; relative ≤ 1e-13 for 99.87%, median 0; at the 20 largest disagreements the kernel is closer to a 50-digit value at 18, with largest relative errors 1.8e-13 (kernel) and 3.9e-13 (reference) | 56.6% |
| `delta` (spot, pips) | 1,000,000 | 1.8e-15 |  | 66.6% |
| `delta` (spot, premium-adjusted) | 1,000,000 | 1.7e-15 |  | 67.5% |
| `delta` (forward, pips) | 1,000,000 | 1.7e-15 |  | 66.3% |
| `delta` (forward, premium-adjusted) | 1,000,000 | 1.7e-15 |  | 67.5% |
| `sabr_vol` | 1,738,720 | 0 |  | 100.0% |
| `strike_from_delta` (pips) | 1,280 | 0 |  | 100.0% |
| `strike_from_delta` (premium-adjusted) | 1,600 | 5.2e-15 | largest ΔK as a share of the root bound 0.25 | 23.1% |
| middle part, adaptive | 1,440 smiles × 3 contracts | k = 1, 2: 1.5e-15; k = 3: 7.6e-12 | normwise 2.6e-15; same n in 100% and same convergence flag in 100% of smiles | 20.8% |
| middle part, n = 1024 | 1,440 smiles × 3 contracts | k = 1, 2: 1.6e-15; k = 3: 3.3e-12 | normwise 2.3e-15 | 22.1% |

<!-- benchmark:end -->
