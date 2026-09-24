import numpy as np
import pytest

from qef.stats.bootstrap import bootstrap_distribution, optimal_block_length, stationary_indices
from qef.stats.hac import long_run_variance, mean_and_se, nw_bandwidth, ols_hac


def _ar1(T, rho, seed):
    rng = np.random.default_rng(seed)
    e = rng.standard_normal(T)
    x = np.empty(T)
    x[0] = e[0]
    for t in range(1, T):
        x[t] = rho * x[t - 1] + e[t]
    return x


def test_hac_se_close_to_iid_se_for_white_noise():
    x = np.random.default_rng(1).standard_normal(4000)
    r = mean_and_se(x)
    assert r["se"] == pytest.approx(x.std() / np.sqrt(x.size), rel=0.1)


def test_hac_lrv_approaches_ar1_long_run_variance():
    # Long-run variance of a unit-innovation AR(1) is 1 / (1 - rho)^2.
    rho = 0.5
    x = _ar1(20000, rho, 2)
    v, lag = long_run_variance(x - x.mean())
    assert lag > 0
    assert v == pytest.approx(1 / (1 - rho) ** 2, rel=0.15)


def test_bandwidth_grows_with_persistence():
    assert nw_bandwidth(_ar1(500, 0.8, 3)) > nw_bandwidth(_ar1(500, 0.0, 3))


def test_ols_hac_matches_least_squares_point_estimates():
    rng = np.random.default_rng(4)
    X = np.column_stack([np.ones(300), rng.standard_normal(300)])
    y = X @ np.array([0.5, -1.0]) + _ar1(300, 0.4, 5)
    r = ols_hac(y, X)
    np.testing.assert_allclose(r["beta"], np.linalg.lstsq(X, y, rcond=None)[0], rtol=1e-12)
    assert np.all(r["se"] > 0)


def test_block_length_larger_for_persistent_series():
    assert optimal_block_length(_ar1(1000, 0.8, 6)) > optimal_block_length(_ar1(1000, 0.0, 6))


def test_stationary_indices_mean_block_length():
    rng = np.random.default_rng(7)
    idx = stationary_indices(100000, 5.0, rng)
    breaks = np.sum(np.diff(idx) % 100000 != 1) + 1
    assert idx.min() >= 0 and idx.max() < 100000
    assert 100000 / breaks == pytest.approx(5.0, rel=0.05)


def test_bootstrap_distribution_of_mean_is_centred():
    x = np.random.default_rng(8).standard_normal(400) + 1.0
    d = bootstrap_distribution(np.mean, [x], B=999)
    assert d.mean() == pytest.approx(x.mean(), abs=0.02)
