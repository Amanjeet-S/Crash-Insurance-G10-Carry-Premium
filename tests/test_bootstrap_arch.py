"""Parity of the stationary bootstrap and its block length with arch 8.0.0 (Sheppard)."""
import numpy as np
import pytest
from scipy.stats import ks_2samp

from qef.stats.bootstrap import _block_length, bootstrap_distribution, optimal_block_length

arch = pytest.importorskip("arch.bootstrap")


def ar1(T, rho, rng, df=5):
    e = rng.standard_t(df, size=T + 100)
    x = np.empty_like(e)
    x[0] = e[0]
    for t in range(1, len(e)):
        x[t] = rho * x[t - 1] + e[t]
    return x[100:]


@pytest.mark.parametrize("T", [56, 104, 160, 500, 2000])
def test_block_length_equals_arch_under_its_two_conventions(T):
    rng = np.random.default_rng(T)
    checked = 0
    for rho in (0.0, 0.3, 0.6, 0.85, 0.95):
        for _ in range(10):
            x = ar1(T, rho, rng)
            ref = float(arch.optimal_block_length(x)["stationary"].iloc[0])
            mine = _block_length(x, arch_conventions=True)
            if ref < 1.0 or mine == 1.0:
                continue  # degenerate cases differ by design (module docstring)
            assert mine == pytest.approx(ref, rel=1e-12)
            checked += 1
    assert checked >= 40


def test_degenerate_long_run_variance_differs_from_arch_by_design():
    # MA(1) with coefficient -0.99: under arch's lag window the kernel estimate of the long-run
    # variance is not positive, the block-length formula does not apply, and this module falls back
    # to the i.i.d. bootstrap (b = 1) while arch squares the estimate and returns its cap (module
    # docstring). Under this module's own window the estimate is positive for this series.
    e = np.random.default_rng(0).standard_normal(105)
    x = e[1:] - 0.99 * e[:-1]
    assert _block_length(x, arch_conventions=True) == 1.0 and optimal_block_length(x) > 1.0
    assert float(arch.optimal_block_length(x)["stationary"].iloc[0]) == np.ceil(min(3 * np.sqrt(104), 104 / 3))


def test_resampler_matches_arch_in_distribution():
    rng = np.random.default_rng(3)
    x = ar1(104, 0.8, rng, df=1000)
    b = 6.0
    mine = bootstrap_distribution(np.mean, [x], B=20000, seed=5, block=[b])
    theirs = arch.StationaryBootstrap(b, x, seed=7).apply(np.mean, reps=20000).ravel()
    # Standard deviation of the bootstrap mean within 3% (Monte Carlo s.e. of the ratio ≈ 0.7%).
    assert mine.std() / theirs.std() == pytest.approx(1.0, abs=0.03)
    assert ks_2samp(mine, theirs).pvalue > 1e-3
