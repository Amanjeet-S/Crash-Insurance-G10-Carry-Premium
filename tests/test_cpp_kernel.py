"""Parity of the C++ kernel (cpp/qef_kernel.cpp, through qef.kernel) with the Python reference.

The module builds the library with cpp/build.sh if it is absent or older than
its sources, and skips, with the reason, if no working C++ compiler is found.
All inputs are synthetic and seeded; no licensed data are read.

Rules fixed on 27 September 2026 before the parity results below were
computed. A pilot on other seeds (one-off scripts, not kept) measured the
disagreement of the two normal distribution functions and motivated the
normwise criteria; the grids, seeds, tolerances and exclusions were then
fixed as stated here and not changed after the tests were run.

Grids (seed 20260927 plus an offset per test). Time to expiry: half of the
points at exactly 1/12 or 1/4, half at 28 to 31 or 89 to 92 days over 365.
Volatility (for SABR, alpha): uniform on [0.03, 0.4]. Forward: a rounded
public level per currency (EUR 1.10, GBP 1.30, AUD 0.70, NZD 0.65, JPY 130,
CHF 0.95, CAD 1.30, NOK 10, SEK 10) times exp(U(-0.3, 0.3)). Base-currency
discount factor: uniform on [0.97, 1.01]. Strikes: from 5-delta to 95-delta,
K = F exp(s²/2 − s Φ⁻¹(Δ)) with Δ uniform on [0.05, 0.95] and s the
volatility times √τ. Option type ±1 with equal probability. Conventions: the
four DeltaConventions and the nine G10 conventions of qef.fx.conventions.
SABR: ρ uniform on (−0.9, 0.9), ν uniform on [0, 4] with one point in ten
exactly 0. Strike targets: |Δ| uniform on [0.05, 0.95] with the sign of the
option type. Middle part: K_min and K_max are the flat-volatility 10-delta
put and call strikes at σ = alpha in the currency's convention, and n is 32,
64, 256 or 1,024 subintervals per side.

Criteria. Relative error |kernel − reference| / |reference| ≤ 1e-13 for the
normal distribution function (x in [−20, 8]), deltas, SABR volatilities and
pips strikes; |ΔΦ⁻¹| ≤ 1e-13 max(1, |x|) for the inverse. Looser or
differently normalised bounds, with their justification:

- Premium: |ΔC| ≤ 1e-13 [F Φ(φd+) + K Φ(φd−)] everywhere, and the relative
  bound where the condition number κ = [F Φ(φd+) + K Φ(φd−)] / |C| is at
  most 10. The kernel's Φ (system erfc) and scipy.special.ndtr differ by
  about 20 units in the last place for |x| < 6 (at most 19 on 2,000,001
  evenly spaced points of [−6, 6] and 21 on 5 million uniform draws, largest
  for x < −4 and at most 2 for x ≥ 0; sample maxima, not a bound, and in
  the lower tail neither function is accurate to better than tens of units
  in the last place), and an out-of-the-money premium is a difference of two
  such terms, so the relative error is that difference times κ, which
  reaches several hundred at 5-delta and low volatility. The normwise
  criterion bounds the error by the size of the terms, which is what either
  implementation can guarantee. A 50-digit evaluation (mpmath) at the points
  of largest disagreement checks that the kernel is accurate there: at the
  20 points of largest relative disagreement its error is at most
  1e-14 [F Φ(φd+) + K Φ(φd−)]. That check does not show how the
  disagreement divides between the two implementations; the benchmark in
  cpp/README.md reports both errors. Amendment after the first run: the
  first version of the test also asserted that more than half of the grid
  has κ ≤ 10. That expectation was wrong, because the premium is a
  difference even at the money, where κ ≈ 2.5/s; only about 3% of the grid
  has κ ≤ 10. I replaced it by the requirement that this set is not empty
  and left every tolerance unchanged; cpp/README.md reports the plain
  relative errors of the whole grid.
- Premium-adjusted strikes (root finding): |ΔK| ≤ 1e-14 F + 1e-14 |K| +
  1e-15 |K|. scipy.optimize.brentq guarantees its root x0 only to
  |x − x0| ≤ xtol + rtol |x0| with the reference's xtol = 1e-14 F and
  rtol = 1e-14; the last term allows 4.5 units in the last place for the
  kernel's own root and the conditioning of the delta equation. Targets
  within 2% below the maximal premium-adjusted call delta are excluded,
  because the equation is ill-conditioned next to the maximiser; there the
  test checks the backward error instead: the reference's delta at the
  kernel's strike is within 1e-14 |target| of the target, on the
  conventional branch K ≥ K*. The maximiser K* itself must agree within
  1e-12 relative: the reference solves for d* to 1e-14 absolute plus 1e-14
  relative, an error of s times that in ln K*. Every input must give an
  error in both implementations or in neither, and invalid inputs must raise
  the reference's exception type.
- Middle part: |ΔI_k| ≤ 1e-13 A_k, with A_k Simpson's rule applied to
  |h_k'' P| on the same nodes, and the relative bound where A_k/|I_k| ≤ 10.
  For k = 1, 2 the weight has one sign on [K_min, K_max], so A_k = |I_k| and
  this is the relative bound; the cubic weight changes sign at K = F, and
  its middle part is a difference of two similar integrals. The adaptive
  rule must stop at the same n with the same converged flag.

Amendment of 27 September 2026, after an independent review. Some rules
above were stated but not asserted everywhere they apply, and one statement
was wrong. (1) The backward-error check on targets within 2% below the
maximal premium-adjusted call delta was made only in
test_strike_near_maximal_premium_adjusted_call_delta; _check_strikes now
applies it, as stated above, to the targets it excludes from the parity
grids (the reference's delta at the kernel's strike within 1e-14 |target|
of the target, and K ≥ K* (1 − 1e-12)). (2) The relative bound of the middle
part where A_k/|I_k| ≤ 10 was asserted only at fixed n; the adaptive test
now asserts it too. (3) The 1e-14 bound of the 50-digit premium check was in
the test from the start but not stated here; it is stated above, and the
test, formerly test_forward_premium_disagreement_is_the_reference_error, is
renamed for what it asserts. (4) The difference of the two Φ was given as up
to about 13 units in the last place, which a finer grid did not reproduce;
the measured figures above replace it. (5) The kernel's bracket check
underflowed for two tiny end values of one sign and is corrected; the case
"PA put product underflow" is added to the invalid inputs. (6) A float n or
n0 in the middle part now raises TypeError, as the reference does, where
the test expected ValueError for n = 2.5, and a real n_max is accepted and
stops the doubling where the reference's does; the test checks both, with
a non-integer n_max below the converged n. (7) Subnormal probabilities were
not among the inputs of the inverse-normal test, and the kernel was wrong
there; the corrected kernel is now tested on 2,004 subnormal p against
ndtri at the same relative tolerance. (8) The module now skips, with the
reason, when the library is older than its sources and no compiler can
rebuild it, and fails with the compiler's output when the build fails. No tolerance, grid, seed or
exclusion was changed, and none was tuned to results.
"""

from __future__ import annotations

import math
import os
import shutil
import subprocess
import warnings

import numpy as np
import pytest
from scipy import special

from qef import kernel as kk
from qef.fx import gk, moments, sabr
from qef.fx.conventions import G10
from qef.fx.smile import SabrSmile

BUILD = kk.ROOT / "cpp" / "build.sh"
SOURCES = [kk.ROOT / "cpp" / "qef_kernel.cpp", kk.ROOT / "cpp" / "qef_kernel.h", BUILD]


def _ensure_library():
    lib = kk.LIBRARY
    stale = not lib.exists() or lib.stat().st_mtime < max(p.stat().st_mtime for p in SOURCES)
    if stale:
        cxx = os.environ.get("CXX", "c++")
        works = shutil.which(cxx) is not None and subprocess.run(
            [cxx, "--version"], capture_output=True).returncode == 0
        if works:
            try:
                subprocess.run(["bash", str(BUILD)], check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                pytest.fail(f"cpp/build.sh failed:\n{e.stdout}{e.stderr}", pytrace=False)
        elif not lib.exists():
            pytest.skip(f"no working C++ compiler ({cxx}) to build the kernel with cpp/build.sh",
                        allow_module_level=True)
        else:
            pytest.skip(f"the kernel library is older than its sources and no working C++ compiler ({cxx}) "
                        "can rebuild it; run cpp/build.sh", allow_module_level=True)
    if not kk.available():
        pytest.fail(f"the kernel library at {lib} does not load", pytrace=False)


_ensure_library()

TOL = kk.PARITY_TOLERANCES
RTOL = TOL["relative"]
SEED = 20260927
CONVENTIONS = [gk.DeltaConvention(s, p) for s in (True, False) for p in (False, True)]
LEVELS = {"EUR": 1.10, "GBP": 1.30, "AUD": 0.70, "NZD": 0.65, "JPY": 130.0, "CHF": 0.95, "CAD": 1.30,
          "NOK": 10.0, "SEK": 10.0}


def _rel(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.abs(a - b) / np.abs(b)
    return np.where(a == b, 0.0, r)


def _grid(rng, n, ccy=None):
    """τ, σ, F, D_b, φ and 5- to 95-delta strikes (module docstring)."""
    tau = np.where(rng.random(n) < 0.5, rng.choice([1 / 12, 0.25], n),
                   np.where(rng.random(n) < 0.5, rng.integers(28, 32, n), rng.integers(89, 93, n)) / 365.0)
    sigma = rng.uniform(0.03, 0.4, n)
    ccys = np.array([ccy] * n) if ccy else rng.choice(list(LEVELS), n)
    F = np.array([LEVELS[c] for c in ccys]) * np.exp(rng.uniform(-0.3, 0.3, n))
    s = sigma * np.sqrt(tau)
    K = F * np.exp(0.5 * s**2 - s * special.ndtri(rng.uniform(0.05, 0.95, n)))
    return {"tau": tau, "sigma": sigma, "F": F, "K": K, "df": rng.uniform(0.97, 1.01, n),
            "phi": rng.choice([-1.0, 1.0], n), "ccy": ccys}


# Normal distribution


def test_library_loads_with_matching_interface():
    assert kk.available()
    assert kk._lib().qef_abi_version() == kk.ABI_VERSION


def test_normal_distribution_function_and_inverse_match_scipy():
    x = np.linspace(-20.0, 8.0, 280_001)
    assert _rel(kk._norm_cdf(x), special.ndtr(x)).max() <= RTOL
    rng = np.random.default_rng(SEED)
    p = np.concatenate([np.exp(rng.uniform(np.log(1e-300), np.log(0.5), 20_000)), rng.uniform(0.0, 1.0, 20_000),
                        1.0 - np.exp(rng.uniform(np.log(1e-16), np.log(0.5), 5_000)), [0.5, 1e-300]])
    ref = special.ndtri(p)
    assert np.max(np.abs(kk._norm_ppf(p) - ref) / np.maximum(1.0, np.abs(ref))) <= RTOL
    assert kk._norm_ppf(0.0) == -np.inf and kk._norm_ppf(1.0) == np.inf and np.isnan(kk._norm_ppf(1.5))
    # Subnormal p (added after the review): scipy's ndtri is accurate there to about 1e-14 absolute.
    tiny = np.finfo(float).tiny
    p = np.concatenate([[5e-324, 1e-323, np.nextafter(tiny, 0.0), tiny],
                        np.exp(rng.uniform(np.log(5e-324), np.log(tiny), 2_000))])
    ref = special.ndtri(p)
    assert np.max(np.abs(kk._norm_ppf(p) - ref) / np.abs(ref)) <= RTOL


# Garman–Kohlhagen premium and delta


def _premium_scale(F, K, sigma, tau, phi):
    d1, d2 = gk.d_plus_minus(F, K, sigma, tau)
    return F * special.ndtr(phi * d1) + K * special.ndtr(phi * d2)


def test_forward_premium_parity():
    g = _grid(np.random.default_rng(SEED + 1), 20_000)
    args = (g["F"], g["K"], g["sigma"], g["tau"], g["phi"])
    ref, out = gk.forward_premium(*args), kk.forward_premium(*args)
    scale = _premium_scale(*args)
    assert np.max(np.abs(out - ref) / scale) <= TOL["premium_normwise"]
    well = scale / np.abs(ref) <= TOL["premium_condition_for_relative"]
    assert well.any()
    assert _rel(out, ref)[well].max() <= RTOL


def test_forward_premium_kernel_accurate_at_largest_disagreements():
    mpmath = pytest.importorskip("mpmath")
    mpmath.mp.dps = 50
    g = _grid(np.random.default_rng(SEED + 1), 20_000)
    args = (g["F"], g["K"], g["sigma"], g["tau"], g["phi"])
    ref, out = gk.forward_premium(*args), kk.forward_premium(*args)
    d1, d2 = gk.d_plus_minus(*args[:4])
    scale = _premium_scale(*args)
    for i in np.argsort(-_rel(out, ref))[:20]:
        # Both implementations form d± identically in double precision; only Φ differs.
        ph = float(g["phi"][i])
        cdf = lambda d: mpmath.ncdf(mpmath.mpf(ph) * mpmath.mpf(float(d)))
        exact = ph * (mpmath.mpf(float(g["F"][i])) * cdf(d1[i]) - mpmath.mpf(float(g["K"][i])) * cdf(d2[i]))
        assert abs(float(mpmath.mpf(float(out[i])) - exact)) <= 1e-14 * scale[i]


def test_forward_premium_deep_out_of_the_money_and_shapes():
    F, tau, sigma = 1.1, 0.25, 0.08
    s = sigma * math.sqrt(tau)
    K = F * np.exp(s * np.array([-10, -8, -6, -4, 4, 6, 8, 10, *special.ndtri([1e-4, 1e-3, 0.999, 0.9999])]))
    for phi in (gk.CALL, gk.PUT):
        ref, out = gk.forward_premium(F, K, sigma, tau, phi), kk.forward_premium(F, K, sigma, tau, phi)
        assert np.max(np.abs(out - ref) / _premium_scale(F, K, sigma, tau, phi)) <= TOL["premium_normwise"]
    scalar = kk.forward_premium(1.1, 1.2, 0.1, 1 / 12, gk.CALL)
    assert type(scalar) is type(gk.forward_premium(1.1, 1.2, 0.1, 1 / 12, gk.CALL)) is np.float64
    K2 = np.linspace(1.0, 1.2, 4)[::2]  # a non-contiguous view
    phis = np.array([[gk.CALL], [gk.PUT], [gk.CALL]])
    ref2, out2 = gk.forward_premium(1.1, K2, 0.1, 0.25, phis), kk.forward_premium(1.1, K2, 0.1, 0.25, phis)
    assert out2.shape == ref2.shape == (3, 2)
    assert _rel(out2, ref2).max() <= RTOL
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bad = (np.array([1.1, 1.1, np.nan]), np.array([0.0, 1.2, 1.2]), np.array([0.1, 0.0, 0.1]), 0.25, gk.CALL)
        np.testing.assert_array_equal(kk.forward_premium(*bad), gk.forward_premium(*bad))
    with pytest.raises(TypeError):
        gk.forward_premium("a", 1.0, 0.1, 0.25, gk.CALL)
    with pytest.raises(TypeError):
        kk.forward_premium("a", 1.0, 0.1, 0.25, gk.CALL)


@pytest.mark.parametrize("convention", CONVENTIONS, ids=lambda c: f"spot={c.spot}-pa={c.premium_adjusted}")
def test_delta_parity(convention):
    g = _grid(np.random.default_rng(SEED + 2), 20_000)
    args = (g["F"], g["K"], g["sigma"], g["tau"], g["phi"], convention, g["df"])
    assert _rel(kk.delta(*args), gk.delta(*args)).max() <= RTOL
    F, tau, sigma = 130.0, 1 / 12, 0.12
    s = sigma * math.sqrt(tau)
    K = F * np.exp(s * np.array([-8.0, -6.0, -4.0, 4.0, 6.0, 8.0]))
    for phi in (gk.CALL, gk.PUT):
        deep = (F, K, sigma, tau, phi, convention, 0.995)
        assert _rel(kk.delta(*deep), gk.delta(*deep)).max() <= RTOL
    fwd = kk.delta(1.1, np.array([1.0, 1.2]), 0.1, 0.25, gk.CALL, convention, np.ones((3, 1)))
    assert fwd.shape == np.shape(gk.delta(1.1, np.array([1.0, 1.2]), 0.1, 0.25, gk.CALL, convention, np.ones((3, 1))))


@pytest.mark.parametrize("ccy", sorted(G10))
def test_delta_parity_by_currency(ccy):
    g = _grid(np.random.default_rng(SEED + 3), 4_000, ccy)
    args = (g["F"], g["K"], g["sigma"], g["tau"], g["phi"], G10[ccy].delta, g["df"])
    assert _rel(kk.delta(*args), gk.delta(*args)).max() <= RTOL


def test_delta_convention_errors_match():
    with pytest.raises(AttributeError):
        gk.delta(1.1, 1.2, 0.1, 0.25, gk.CALL, None)
    with pytest.raises(AttributeError):
        kk.delta(1.1, 1.2, 0.1, 0.25, gk.CALL, None)


# SABR


def _sabr_reference(K, F, tau, alpha, rho, nu):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.array([sabr.sabr_vol(k, f, t, a, r, v) for k, f, t, a, r, v in
                         np.broadcast(K, F, tau, alpha, rho, nu)]).reshape(np.broadcast(K, F, tau, alpha, rho, nu).shape)


def test_sabr_vol_parity():
    rng = np.random.default_rng(SEED + 4)
    g = _grid(rng, 20_000)
    alpha = g["sigma"]
    s = alpha * np.sqrt(g["tau"])
    K = g["F"] * np.exp(0.5 * s**2 - s * special.ndtri(rng.uniform(0.05, 0.95, alpha.size)))
    rho = rng.uniform(-0.9, 0.9, alpha.size)
    nu = np.where(rng.random(alpha.size) < 0.1, 0.0, rng.uniform(0.0, 4.0, alpha.size))
    args = (K, g["F"], g["tau"], alpha, rho, nu)
    assert _rel(kk.sabr_vol(*args), _sabr_reference(*args)).max() <= RTOL
    one = (K[:1000], 1.1, 0.25, 0.1, -0.3, 1.5)  # one smile, array strikes, as the reference is used
    assert _rel(kk.sabr_vol(*one), sabr.sabr_vol(*one)).max() <= RTOL


def test_sabr_vol_edge_cases():
    F, tau, alpha, rho, nu = 1.1, 1 / 12, 0.1, -0.4, 2.0
    c = nu / alpha
    # z = c ln(F/K) on both sides of the branch point |z| = 1e-7, at 0 and across scales.
    zs = np.concatenate([[0.0], np.outer([1.0, -1.0], [1e-12, 1e-9, 5e-8, 1e-7 * (1 - 1e-9), 1e-7,
                                                       1e-7 * (1 + 1e-9), 2e-7, 1e-5, 1e-3]).ravel()])
    K = F * np.exp(-zs / c)
    z_ref = c * np.log(F / K)
    assert (np.abs(z_ref) < 1e-7).any() and (np.abs(z_ref) >= 1e-7).any()
    for r in (rho, 0.0, 0.999, -0.999):
        assert _rel(kk.sabr_vol(K, F, tau, alpha, r, nu), sabr.sabr_vol(K, F, tau, alpha, r, nu)).max() <= RTOL
    K2 = F * np.array([1e-3, 0.2, 0.5, 2.0, 5.0, 1e3])
    for v in (0.0, 1e-300, 4.0):
        assert _rel(kk.sabr_vol(K2, F, tau, alpha, 0.5, v), sabr.sabr_vol(K2, F, tau, alpha, 0.5, v)).max() <= RTOL
    np.testing.assert_array_equal(kk.sabr_vol(K2, F, tau, alpha, 0.5, 0.0), alpha)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bad = np.array([0.0, -1.0, np.nan, np.inf, 1.0])
        np.testing.assert_array_equal(kk.sabr_vol(bad, F, tau, alpha, 0.5, nu), sabr.sabr_vol(bad, F, tau, alpha, 0.5, nu))
        zero = np.float64(0.0)  # a Python-float zero makes the reference raise ZeroDivisionError (kernel.py)
        np.testing.assert_array_equal(kk.sabr_vol(K2, F, tau, zero, 0.5, nu), sabr.sabr_vol(K2, F, tau, zero, 0.5, nu))
    assert type(kk.sabr_vol(1.2, F, tau, alpha, rho, nu)) is type(sabr.sabr_vol(1.2, F, tau, alpha, rho, nu))
    with pytest.raises(ValueError):
        sabr.sabr_vol("x", F, tau, alpha, rho, nu)
    with pytest.raises(ValueError):
        kk.sabr_vol("x", F, tau, alpha, rho, nu)
    with pytest.raises(NotImplementedError):
        kk.sabr_vol(1.2, F, tau, alpha, rho, nu, beta=0.5)


# Strike from delta


def _strike_reference(t, F, sigma, tau, phi, convention, df):
    K, err = np.full(t.size, np.nan), np.zeros(t.size, dtype=bool)
    for i in range(t.size):
        try:
            K[i] = gk.strike_from_delta(t[i], F[i], sigma[i], tau[i], phi[i], convention, df[i])
        except ValueError:
            err[i] = True
    return K, err


def _strike_case(seed, n, convention, ccy=None):
    rng = np.random.default_rng(seed)
    g = _grid(rng, n, ccy)
    t = g["phi"] * rng.uniform(0.05, 0.95, n)
    args = (t, g["F"], g["sigma"], g["tau"], g["phi"], convention, g["df"])
    ref, err = _strike_reference(*args)
    out, status = kk.strike_from_delta_status(*args)
    np.testing.assert_array_equal(status != 0, err)
    assert np.isnan(out[err]).all()
    return g, t, ref, out, ~err


def _root_bound(F, K_ref):
    return TOL["root_xtol_per_unit_forward"] * F + (TOL["root_rtol"] + TOL["root_slack"]) * np.abs(K_ref)


def _check_strikes(convention, g, t, ref, out, ok):
    if not convention.premium_adjusted:
        assert _rel(out[ok], ref[ok]).max() <= RTOL
        return
    # Exclude targets within 2% below the maximal premium-adjusted call delta.
    calls = ok & (g["phi"] > 0)
    far = ok.copy()
    k_star = np.full(t.size, np.nan)
    for i in np.flatnonzero(calls):
        k_star[i] = gk.pa_call_delta_maximiser(g["F"][i], g["sigma"][i], g["tau"][i])
        top = gk.delta(g["F"][i], k_star[i], g["sigma"][i], g["tau"][i], gk.CALL, convention, g["df"][i])
        far[i] = t[i] <= 0.98 * top
    assert far.sum() > 0.8 * ok.sum()
    assert np.all(np.abs(out[far] - ref[far]) <= _root_bound(g["F"][far], ref[far]))
    assert _rel(out[far], ref[far]).max() <= RTOL
    # The excluded targets: backward error and the conventional branch (amendment (1) of the docstring).
    for i in np.flatnonzero(ok & ~far):
        d = gk.delta(g["F"][i], out[i], g["sigma"][i], g["tau"][i], gk.CALL, convention, g["df"][i])
        assert abs(d - t[i]) <= 1e-14 * abs(t[i])
        assert out[i] >= k_star[i] * (1 - 1e-12)


@pytest.mark.parametrize("convention", CONVENTIONS, ids=lambda c: f"spot={c.spot}-pa={c.premium_adjusted}")
def test_strike_from_delta_parity(convention):
    n = 1_500 if convention.premium_adjusted else 5_000
    g, t, ref, out, ok = _strike_case(SEED + 5, n, convention)
    if convention.premium_adjusted:
        assert (~ok).sum() > 0  # targets above the maximal call delta are in the grid
    _check_strikes(convention, g, t, ref, out, ok)


@pytest.mark.parametrize("ccy", sorted(G10))
def test_strike_from_delta_parity_by_currency(ccy):
    conv = G10[ccy].delta
    g, t, ref, out, ok = _strike_case(SEED + 6, 400 if conv.premium_adjusted else 1_000, conv, ccy)
    _check_strikes(conv, g, t, ref, out, ok)


@pytest.mark.parametrize("convention", [c for c in CONVENTIONS if c.premium_adjusted],
                         ids=lambda c: f"spot={c.spot}")
def test_strike_near_maximal_premium_adjusted_call_delta(convention):
    for F, sigma, tau, df in ((1.1, 0.4, 0.25, 0.98), (130.0, 0.05, 1 / 12, 1.0), (10.0, 0.2, 0.25, 1.005)):
        k_star = gk.pa_call_delta_maximiser(F, sigma, tau)
        assert _rel(kk.pa_call_delta_maximiser(F, sigma, tau), k_star) <= 1e-12
        top = gk.delta(F, k_star, sigma, tau, gk.CALL, convention, df)
        for j in range(2, 7):
            target = (1.0 - 10.0**-j) * top
            K = kk.strike_from_delta(target, F, sigma, tau, gk.CALL, convention, df)
            assert abs(gk.delta(F, K, sigma, tau, gk.CALL, convention, df) - target) <= 1e-14 * target
            assert K >= k_star * (1 - 1e-12)
        above = 1.0001 * top
        with pytest.raises(ValueError):
            gk.strike_from_delta(above, F, sigma, tau, gk.CALL, convention, df)
        with pytest.raises(ValueError):
            kk.strike_from_delta(above, F, sigma, tau, gk.CALL, convention, df)


PIPS, PA = gk.DeltaConvention(True, False), gk.DeltaConvention(True, True)
INVALID = {
    "sign mismatch": (0.25, 1.1, 0.1, 0.25, gk.PUT, PIPS, 1.0),
    "zero target": (0.0, 1.1, 0.1, 0.25, gk.CALL, PA, 1.0),
    "NaN target": (np.nan, 1.1, 0.1, 0.25, gk.CALL, PIPS, 1.0),
    "phi 2": (0.25, 1.1, 0.1, 0.25, 2, PIPS, 1.0),
    "zero volatility": (0.25, 1.1, 0.0, 0.25, gk.CALL, PIPS, 1.0),
    "negative volatility": (-0.25, 1.1, -0.1, 0.25, gk.PUT, PA, 1.0),
    "infinite volatility": (0.25, 1.1, np.inf, 0.25, gk.CALL, PA, 1.0),
    "NaN expiry": (0.25, 1.1, 0.1, np.nan, gk.CALL, PA, 1.0),
    "zero expiry": (0.25, 1.1, 0.1, 0.0, gk.CALL, PIPS, 1.0),
    "pips delta of one": (1.0, 1.1, 0.1, 0.25, gk.CALL, PIPS, 1.0),
    "pips delta above D_b": (0.99, 1.1, 0.1, 0.25, gk.CALL, PIPS, 0.98),
    "pips NaN D_b": (0.25, 1.1, 0.1, 0.25, gk.CALL, PIPS, np.nan),
    "pips zero D_b (NumPy)": (np.float64(0.25), 1.1, 0.1, 0.25, gk.CALL, PIPS, np.float64(0.0)),
    "PA call above maximum": (0.999, 1.1, 0.1, 0.25, gk.CALL, PA, 1.0),
    "PA call infinite target": (np.inf, 1.1, 0.1, 0.25, gk.CALL, PA, 1.0),
    "PA put below bracket": (-1e-30, 1.1, 0.1, 0.25, gk.PUT, PA, 1.0),
    "PA put infinite target": (-np.inf, 1.1, 0.1, 0.25, gk.PUT, PA, 1.0),
    "PA put NaN D_b": (-0.25, 1.1, 0.1, 0.25, gk.PUT, PA, np.nan),
    "PA call NaN D_b": (0.25, 1.1, 0.1, 0.25, gk.CALL, PA, np.nan),
    "PA put negative D_b": (-0.25, 1.1, 0.1, 0.25, gk.PUT, PA, -1.0),
    "PA call negative D_b": (0.25, 1.1, 0.1, 0.25, gk.CALL, PA, -1.0),
    "PA put NaN forward": (-0.25, np.nan, 0.1, 0.25, gk.PUT, PA, 1.0),
    "PA call NaN forward": (0.25, np.nan, 0.1, 0.25, gk.CALL, PA, 1.0),
    "PA put negative forward": (-0.25, -1.0, 0.1, 0.25, gk.PUT, PA, 1.0),
    "PA call negative forward": (0.25, -1.0, 0.1, 0.25, gk.CALL, PA, 1.0),
    "PA put zero forward (NumPy)": (-0.25, np.float64(0.0), 0.1, 0.25, gk.PUT, PA, 1.0),
    "PA call huge volatility": (0.25, 1.0, 50.0, 1.0, gk.CALL, PA, 1.0),
    "PA put product underflow": (-1e-300, 1.1, 0.1, 0.25, gk.PUT, PA, 1e-152),
    "string target": ("a", 1.1, 0.1, 0.25, gk.CALL, PIPS, 1.0),
}


@pytest.mark.parametrize("case", list(INVALID))
def test_strike_from_delta_invalid_inputs_raise_like_the_reference(case):
    args = INVALID[case]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(Exception) as ref:
            gk.strike_from_delta(*args)
        with pytest.raises(Exception) as out:
            kk.strike_from_delta(*args)
    # The reference raises NumPy's UFuncTypeError, a TypeError, on non-numeric input; the wrapper raises TypeError.
    assert out.type is ref.type or (out.type is TypeError and issubclass(ref.type, TypeError)), (ref, out)


def test_strike_from_delta_arrays_raise_at_the_first_invalid_element():
    t = np.array([0.25, 0.5, -0.1, 0.25])
    K, status = kk.strike_from_delta_status(t, 1.1, 0.1, 0.25, gk.CALL, PIPS)
    np.testing.assert_array_equal(status, [0, 0, 1, 0])
    assert np.isnan(K[2]) and np.isfinite(K[[0, 1, 3]]).all()
    with pytest.raises(ValueError, match="sign of the option type"):
        kk.strike_from_delta(t, 1.1, 0.1, 0.25, gk.CALL, PIPS)
    with pytest.raises(ValueError):
        kk.pa_call_delta_maximiser(1.1, np.array([0.1, 0.0]), 0.25)


# Middle part of the moment contracts


def _abs_simpson(F, K_min, K_max, tau, vol, usd, n):
    """A_k: Simpson's rule for |h_k'' P| on the reference's nodes."""
    total = np.zeros(len(moments.CONTRACTS))
    for a, b in ((K_min, F), (F, K_max)):
        K = np.linspace(a, b, n + 1)
        p = moments.otm_price(K, F, tau, vol)
        total += [moments._simpson(np.abs(moments.contract_weight(k, K, F, usd) * p), (b - a) / n)
                  for k in moments.CONTRACTS]
    return total


def _smiles(seed, m):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(m):
        ccy = str(rng.choice(list(LEVELS)))
        conv = G10[ccy]
        tau = float(rng.choice([1 / 12, 0.25])) if rng.random() < 0.5 else float(
            rng.integers(28, 32) if rng.random() < 0.5 else rng.integers(89, 93)) / 365.0
        F = LEVELS[ccy] * math.exp(rng.uniform(-0.3, 0.3))
        alpha, rho = rng.uniform(0.03, 0.4), rng.uniform(-0.9, 0.9)
        nu = 0.0 if rng.random() < 0.1 else rng.uniform(0.0, 4.0)
        df = rng.uniform(0.97, 1.01)
        K_min = gk.strike_from_delta(-0.10, F, alpha, tau, gk.PUT, conv.delta, df)
        K_max = gk.strike_from_delta(0.10, F, alpha, tau, gk.CALL, conv.delta, df)
        out.append((F, K_min, K_max, tau, SabrSmile(F, tau, alpha, rho, nu), conv.usd_base, int(rng.choice([32, 64, 256, 1024]))))
    return out


def test_middle_contracts_fixed_n_parity():
    worst_rel_well, worst_normwise, n_usd = 0.0, 0.0, 0
    for F, K_min, K_max, tau, smile, usd, n in _smiles(SEED + 7, 300):
        ref = moments.middle_contracts(F, K_min, K_max, tau, smile.vol, usd, n0=n // 2, n_max=n)
        assert ref["n"] == n
        out = kk.middle_contracts_fixed(F, K_min, K_max, tau, smile.vol, usd, n)
        A = _abs_simpson(F, K_min, K_max, tau, smile.vol, usd, n)
        diff = np.abs(out - ref["values"])
        worst_normwise = max(worst_normwise, float(np.max(diff / A)))
        well = A / np.abs(ref["values"]) <= TOL["premium_condition_for_relative"]
        assert well[:2].all()  # k = 1, 2: weights of one sign
        worst_rel_well = max(worst_rel_well, float(_rel(out, ref["values"])[well].max()))
        n_usd += usd
    assert worst_normwise <= TOL["middle_normwise"]
    assert worst_rel_well <= RTOL
    assert 0 < n_usd < 300


def test_middle_contracts_adaptive_parity_and_node_reuse():
    for F, K_min, K_max, tau, smile, usd, _ in _smiles(SEED + 8, 60):
        ref = moments.middle_contracts(F, K_min, K_max, tau, smile.vol, usd)
        out = kk.middle_contracts(F, K_min, K_max, tau, smile.vol, usd)
        assert out["n"] == ref["n"] and out["converged"] == ref["converged"]
        assert type(out["n"]) is int and type(out["rel_change"]) is float and type(out["converged"]) is bool
        A = _abs_simpson(F, K_min, K_max, tau, smile.vol, usd, ref["n"])
        assert np.max(np.abs(out["values"] - ref["values"]) / A) <= TOL["middle_normwise"]
        well = A / np.abs(ref["values"]) <= TOL["premium_condition_for_relative"]
        assert well[:2].all() and _rel(out["values"], ref["values"])[well].max() <= RTOL
        kappa = float(np.max(A / np.abs(ref["values"])))
        assert abs(out["rel_change"] - ref["rel_change"]) <= 3 * TOL["middle_normwise"] * kappa
        # The adaptive rule reuses the previous level's nodes; it must equal a fresh fixed-n evaluation bitwise.
        np.testing.assert_array_equal(out["values"], kk.middle_contracts_fixed(F, K_min, K_max, tau, smile, usd, out["n"]))


def test_middle_contracts_batch_equals_single_calls():
    smiles = _smiles(SEED + 9, 40)
    cols = list(zip(*[(F, K_min, K_max, tau, s.alpha, s.rho, s.nu, usd) for F, K_min, K_max, tau, s, usd, _ in smiles]))
    batch = kk.middle_contracts_sabr(*map(np.array, cols))
    fixed = kk.middle_contracts_sabr(*map(np.array, cols), n=128)
    for i, (F, K_min, K_max, tau, s, usd, _) in enumerate(smiles):
        single = kk.middle_contracts(F, K_min, K_max, tau, s.vol, usd)
        np.testing.assert_array_equal(batch["values"][i], single["values"])
        assert batch["n"][i] == single["n"]
        np.testing.assert_array_equal(fixed["values"][i], kk.middle_contracts_fixed(F, K_min, K_max, tau, s, usd, 128))


def test_middle_contracts_invalid_inputs():
    smile = SabrSmile(1.1, 1 / 12, 0.1, -0.2, 1.0)
    for bounds in ((1.2, 1.3), (1.0, 1.1), (1.1, 1.2), (np.nan, 1.2)):
        with pytest.raises(ValueError):
            moments.middle_contracts(1.1, *bounds, 1 / 12, smile.vol, False)
        with pytest.raises(ValueError):
            kk.middle_contracts(1.1, *bounds, 1 / 12, smile.vol, False)
    with pytest.raises(ValueError):
        kk.middle_contracts_sabr([1.1, 1.1], [1.0, 1.2], 1.2, 1 / 12, 0.1, -0.2, 1.0, False)
    with pytest.raises(TypeError):
        kk.middle_contracts(1.1, 1.0, 1.2, 1 / 12, lambda K: 0.1 + 0 * K, False)
    with pytest.raises(NotImplementedError):
        kk.middle_contracts(1.1, 1.0, 1.2, 1 / 12, SabrSmile(1.1, 1 / 12, 0.1, -0.2, 1.0, beta=0.5), False)
    for bad in (0, -4, 2**22 + 1):
        with pytest.raises(ValueError):
            kk.middle_contracts_fixed(1.1, 1.0, 1.2, 1 / 12, smile, False, bad)
    with pytest.raises(TypeError):
        kk.middle_contracts_fixed(1.1, 1.0, 1.2, 1 / 12, smile, False, 2.5)
    for kw in ({"n0": 2.5}, {"n0": 32.0}):
        with pytest.raises(TypeError):
            moments.middle_contracts(1.1, 1.0, 1.2, 1 / 12, smile.vol, False, **kw)
        with pytest.raises(TypeError):
            kk.middle_contracts(1.1, 1.0, 1.2, 1 / 12, smile, False, **kw)
    for n_max in (65536.0, 1000.5, np.int64(512), -5, 128.5, 64.5):
        ref = moments.middle_contracts(1.1, 1.0, 1.2, 1 / 12, smile.vol, False, n_max=n_max)
        out = kk.middle_contracts(1.1, 1.0, 1.2, 1 / 12, smile, False, n_max=n_max)
        assert out["n"] == ref["n"] and out["converged"] == ref["converged"]
    # A non-integer n_max below the converged n must round up, not down: the reference stops at 256 here.
    assert moments.middle_contracts(1.1, 1.0, 1.2, 1 / 12, smile.vol, False, n_max=128.5)["n"] == 256
    for n_max in (np.nan, np.inf, 2**22 + 1):
        with pytest.raises(ValueError):
            kk.middle_contracts(1.1, 1.0, 1.2, 1 / 12, smile, False, n_max=n_max)
    with pytest.raises(ValueError):
        kk.middle_contracts(1.1, 1.0, 1.2, 1 / 12, smile, False, n0=0)
    # n_max below the first doubled n stops there, as in the reference.
    ref = moments.middle_contracts(1.1, 1.0, 1.2, 1 / 12, smile.vol, False, n0=16, n_max=0)
    out = kk.middle_contracts(1.1, 1.0, 1.2, 1 / 12, smile, False, n0=16, n_max=0)
    assert out["n"] == ref["n"] == 32
