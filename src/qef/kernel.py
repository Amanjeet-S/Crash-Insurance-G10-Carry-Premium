"""ctypes bindings to the C++ kernel in cpp/qef_kernel.cpp (research design, Stage 5).

Each wrapper takes the arguments of the Python function it mirrors, in the
same order and with the same defaults, and returns what that function
returns: an array of the broadcast shape, or a NumPy scalar for scalar input.

==========================  =========================================
wrapper                     Python reference
==========================  =========================================
forward_premium             qef.fx.gk.forward_premium
delta                       qef.fx.gk.delta (every DeltaConvention)
sabr_vol                    qef.fx.sabr.sabr_vol (beta = 1 only)
pa_call_delta_maximiser     qef.fx.gk.pa_call_delta_maximiser
strike_from_delta           qef.fx.gk.strike_from_delta
middle_contracts            qef.fx.moments.middle_contracts (SABR smiles)
middle_contracts_fixed      its inner integrate(n) at a fixed n
middle_contracts_sabr       both, for many smiles in one call
==========================  =========================================

The library is built by ``bash cpp/build.sh``. Importing this module never
fails; the first call raises KernelUnavailableError, with the build command,
if the library is absent or was built from another version of the interface.

Errors. The kernel returns a status per element, and the wrapper raises the
exception type that the reference raises for that input, at the first failing
element in C order: ValueError for every invalid input the reference checks
(signs, volatility and expiry, the pips range, the maximal premium-adjusted
call delta, the bracket and NaN checks of scipy.optimize.brentq, and
K_min < F < K_max), RuntimeError if the root search fails to converge (as
brentq does) and TypeError for non-numeric input (where the reference's
arithmetic raises a TypeError). The kernel's own limits raise ValueError
(subintervals, below) or MemoryError (node buffers).

Differences from the reference, all outside its documented use:

- strike_from_delta and pa_call_delta_maximiser accept arrays and work
  element by element; sabr_vol accepts arrays of parameters.
- sabr_vol raises NotImplementedError for beta != 1, and the middle-part
  wrappers accept only a qef.fx.smile.SabrSmile with beta = 1 (or its ``vol``
  method) as ``vol_of_strike``; other callables raise TypeError.
- A zero passed as a Python float where the reference divides by it (the
  forward in delta and strike_from_delta, the discount factor in the pips
  strike, alpha in sabr_vol) makes the reference raise ZeroDivisionError
  from float division; with NumPy inputs the reference, like the kernel,
  returns inf or NaN or raises ValueError.
- middle_contracts raises ValueError for n0 < 1 (the reference raises
  ValueError from numpy.linspace for a negative n0, and for n0 = 0 raises
  ZeroDivisionError with Python-float bounds or never terminates with NumPy
  bounds), for n0 or n_max above 2**22, which bounds the kernel's buffers
  (the reference's default n_max is 2**16), and for a NaN or infinite n_max
  (the reference's test n >= n_max is never true for a NaN n_max, so it
  stops only on convergence, and always true for n_max = -inf). As in the
  reference, a non-integer n or n0 raises TypeError (from numpy.linspace
  there) and any other real n_max is accepted; the wrapper passes
  max(ceil(n_max), 0) to the kernel, which stops the doubling where the
  reference's n >= n_max does. (Amended on 27 September 2026 after an
  independent review: the wrapper accepted only integer n_max, raised
  ValueError for a non-integer n or n0, and let ctypes wrap an n_max below
  -2**63.)

PARITY_TOLERANCES records the parity criteria of tests/test_cpp_kernel.py
and scripts/benchmark_kernel.py; the test module's docstring justifies them.
"""

from __future__ import annotations

import ctypes
import math
import numbers
import operator
import sys
from pathlib import Path

import numpy as np

from .fx.moments import CONTRACTS, SIMPSON_RTOL

ABI_VERSION = 1
ROOT = Path(__file__).resolve().parents[2]
LIBRARY = ROOT / "cpp" / "build" / ("libqef_kernel.dylib" if sys.platform == "darwin" else "libqef_kernel.so")
BUILD_COMMAND = "bash cpp/build.sh"
MAX_SUBINTERVALS = 2**22  # bounds n, n0 and n_max, so a level has at most 2 max(n0, n_max) <= 2**23 subintervals

PARITY_TOLERANCES = {
    "relative": 1e-13,  # |kernel - reference| / |reference|, every kernel unless stated below
    "premium_normwise": 1e-13,  # |ΔC| / [F Φ(φd+) + K Φ(φd−)], all premia
    "premium_condition_for_relative": 10.0,  # relative bound asserted where [F Φ(φd+) + K Φ(φd−)] / |C| ≤ 10
    "root_xtol_per_unit_forward": 1e-14,  # premium-adjusted strikes: |ΔK| ≤ 1e-14 F + 1e-14 |K| + 1e-15 |K|
    "root_rtol": 1e-14,
    "root_slack": 1e-15,
    "middle_normwise": 1e-13,  # |ΔI_k| / Simpson integral of |h_k'' P|, all contracts
}

_ERRORS = {
    1: (ValueError, "target delta must be non-zero with the sign of the option type"),
    2: (ValueError, "volatility and time to expiry must be positive and finite"),
    3: (ValueError, "|delta| must lie in (0, D_b) for pips spot delta"),
    4: (ValueError, "target exceeds the maximal premium-adjusted call delta"),
    5: (ValueError, "f(a) and f(b) must have different signs"),
    6: (ValueError, "a function value in the root search is NaN; solver cannot continue"),
    7: (ValueError, "xtol too small: the forward must be positive"),
    8: (RuntimeError, "the root search failed to converge"),
    9: (ValueError, "the middle part needs K_min < F < K_max"),
    10: (ValueError, "n and n0 must be integers in [1, 2**22] and n_max a finite number at most 2**22"),
    11: (MemoryError, "the kernel could not allocate its node buffers"),
}


class KernelUnavailableError(ImportError):
    """The shared library is missing or does not match this module."""


_lib_cache = None
_D = ctypes.POINTER(ctypes.c_double)
_I32 = ctypes.POINTER(ctypes.c_int32)
_I64 = ctypes.POINTER(ctypes.c_int64)
_N = ctypes.c_int64


def _signatures(lib):
    sig = {
        "qef_norm_cdf": [_N, _D, _I64, _D],
        "qef_norm_ppf": [_N, _D, _I64, _D],
        "qef_forward_premium": [_N, _D, _D, _D, _D, _D, _I64, _D],
        "qef_delta": [_N, _D, _D, _D, _D, _D, _D, _I64, ctypes.c_int32, ctypes.c_int32, _D],
        "qef_sabr_vol": [_N, _D, _D, _D, _D, _D, _D, _I64, _D],
        "qef_pa_call_delta_maximiser": [_N, _D, _D, _D, _I64, _D, _I32],
        "qef_strike_from_delta": [_N, _D, _D, _D, _D, _D, _D, _I64, ctypes.c_int32, ctypes.c_int32, _D, _I32],
        "qef_middle_contracts_fixed": [_N, *[_D] * 10, _I64, _N, _D, _I32],
        "qef_middle_contracts_adaptive": [_N, *[_D] * 10, _I64, _N, ctypes.c_double, _N, _D, _I64, _D, _I32, _I32],
    }
    for name, args in sig.items():
        fn = getattr(lib, name)
        fn.argtypes, fn.restype = args, None
    lib.qef_abi_version.argtypes, lib.qef_abi_version.restype = [], ctypes.c_int32


def _lib():
    global _lib_cache
    if _lib_cache is None:
        if not LIBRARY.exists():
            raise KernelUnavailableError(f"the C++ kernel is not built: run `{BUILD_COMMAND}` (expected {LIBRARY})")
        lib = ctypes.CDLL(str(LIBRARY))
        lib.qef_abi_version.restype = ctypes.c_int32
        if lib.qef_abi_version() != ABI_VERSION:
            raise KernelUnavailableError(f"{LIBRARY} was built for another interface version: run `{BUILD_COMMAND}`")
        _signatures(lib)
        _lib_cache = lib
    return _lib_cache


def available() -> bool:
    """Whether the shared library is built and matches this module."""
    try:
        _lib()
    except (KernelUnavailableError, OSError):
        return False
    return True


def _inputs(*xs):
    """Broadcast float64 inputs: (shape, size, arrays, increments), scalars passed with increment 0."""
    arrs = []
    for x in xs:
        a = np.asarray(x)
        if a.dtype.kind not in "biuf":  # the reference's arithmetic raises TypeError on such input
            raise TypeError(f"unsupported input of dtype {a.dtype}: the kernel takes real numbers")
        arrs.append(a.astype(np.float64, copy=False))
    shape = np.broadcast_shapes(*(a.shape for a in arrs))
    size = math.prod(shape)
    out, inc = [], []
    for a in arrs:
        if a.size == 1:
            out.append(np.ascontiguousarray(a.reshape(1)))
            inc.append(0)
        else:
            out.append(np.ascontiguousarray(np.broadcast_to(a, shape)).reshape(-1))
            inc.append(1)
    return shape, size, out, np.array(inc, dtype=np.int64)


def _p(a):
    return a.ctypes.data_as(_I32 if a.dtype == np.int32 else _I64 if a.dtype == np.int64 else _D)


def _result(out, shape):
    out = out.reshape(shape)
    return out[()] if out.ndim == 0 else out


def _raise_first(status):
    bad = np.flatnonzero(status)
    if bad.size:
        exc, msg = _ERRORS[int(status.flat[bad[0]])]
        raise exc(msg)


def _norm_cdf(x):
    """Φ(x) = erfc(−x/√2)/2 (the kernel's normal distribution function; scipy.stats.norm.cdf in the reference)."""
    shape, n, (a,), inc = _inputs(x)
    out = np.empty(n)
    _lib().qef_norm_cdf(n, _p(a), _p(inc), _p(out))
    return _result(out, shape)


def _norm_ppf(p):
    """Φ⁻¹(p) (the kernel's inverse normal; scipy.stats.norm.ppf in the reference)."""
    shape, n, (a,), inc = _inputs(p)
    out = np.empty(n)
    _lib().qef_norm_ppf(n, _p(a), _p(inc), _p(out))
    return _result(out, shape)


def forward_premium(F, K, sigma, tau, phi):
    """Undiscounted Garman–Kohlhagen premium φ[F Φ(φd+) − K Φ(φd−)] (qef.fx.gk.forward_premium)."""
    shape, n, arrs, inc = _inputs(F, K, sigma, tau, phi)
    out = np.empty(n)
    _lib().qef_forward_premium(n, *map(_p, arrs), _p(inc), _p(out))
    return _result(out, shape)


def delta(F, K, sigma, tau, phi, convention, df_base=1.0):
    """Delta in the stated convention (qef.fx.gk.delta); ``df_base`` is ignored for forward deltas."""
    spot, pa = bool(convention.spot), bool(convention.premium_adjusted)
    shape, n, arrs, inc = _inputs(F, K, sigma, tau, phi, df_base if spot else 1.0)
    out = np.empty(n)
    _lib().qef_delta(n, *map(_p, arrs), _p(inc), int(spot), int(pa), _p(out))
    return _result(out, shape)


def sabr_vol(K, F, tau, alpha, rho, nu, beta=1.0):
    """Black implied volatility of the Hagan SABR expansion at beta = 1 (qef.fx.sabr.sabr_vol)."""
    if np.any(np.asarray(beta, dtype=float) != 1.0):
        raise NotImplementedError("the C++ kernel implements beta = 1 only; use qef.fx.sabr.sabr_vol")
    K = np.asarray(K, dtype=float)  # as the reference converts K, raising ValueError for non-numeric strikes
    shape, n, arrs, inc = _inputs(K, F, tau, alpha, rho, nu)
    out = np.empty(n)
    _lib().qef_sabr_vol(n, *map(_p, arrs), _p(inc), _p(out))
    return _result(out, shape)


def pa_call_delta_maximiser(F, sigma, tau):
    """Strike K* maximising the premium-adjusted call delta (qef.fx.gk.pa_call_delta_maximiser)."""
    shape, n, arrs, inc = _inputs(F, sigma, tau)
    out, status = np.empty(n), np.empty(n, dtype=np.int32)
    _lib().qef_pa_call_delta_maximiser(n, *map(_p, arrs), _p(inc), _p(out), _p(status))
    _raise_first(status)
    return _result(out, shape)


def strike_from_delta_status(target, F, sigma, tau, phi, convention, df_base=1.0):
    """strike_from_delta without raising: (strikes, status), NaN where status is non-zero.

    Status codes are those of cpp/qef_kernel.h; ``STATUS_MESSAGES`` maps
    each to the exception type and message the wrapper raises.
    """
    spot, pa = bool(convention.spot), bool(convention.premium_adjusted)
    shape, n, arrs, inc = _inputs(target, F, sigma, tau, phi, df_base if spot else 1.0)
    out, status = np.empty(n), np.empty(n, dtype=np.int32)
    _lib().qef_strike_from_delta(n, *map(_p, arrs), _p(inc), int(spot), int(pa), _p(out), _p(status))
    return _result(out, shape), status.reshape(shape)


def strike_from_delta(target, F, sigma, tau, phi, convention, df_base=1.0):
    """Strike with the given delta under a flat volatility (qef.fx.gk.strike_from_delta)."""
    K, status = strike_from_delta_status(target, F, sigma, tau, phi, convention, df_base)
    _raise_first(status)
    return K


STATUS_MESSAGES = dict(_ERRORS)


def _smile_parameters(vol_of_strike):
    """(F, tau, alpha, rho, nu) of a SabrSmile with beta = 1, given the smile or its ``vol`` method."""
    from .fx.smile import SabrSmile

    smile = getattr(vol_of_strike, "__self__", vol_of_strike)
    if not isinstance(smile, SabrSmile):
        raise TypeError("the C++ kernel integrates SABR smiles only: pass a qef.fx.smile.SabrSmile or its vol method")
    if smile.beta != 1.0:
        raise NotImplementedError("the C++ kernel implements beta = 1 only")
    return smile.F, smile.tau, smile.alpha, smile.rho, smile.nu


def _check_subintervals(n):
    """n (or n0) as an int in [1, MAX_SUBINTERVALS]; TypeError for a non-integer, as numpy.linspace raises."""
    n = operator.index(n)
    if not 1 <= n <= MAX_SUBINTERVALS:
        raise ValueError(_ERRORS[10][1])
    return n


def _check_n_max(n_max):
    """max(ceil(n_max), 0) for a finite real n_max at most MAX_SUBINTERVALS (module docstring).

    For integer n >= 1, n >= n_max holds exactly when n >= ceil(n_max), and
    always when n_max <= 0, so the kernel's integer test against the returned
    value stops where the reference's does. Non-real input raises TypeError,
    as the reference's comparison does.
    """
    if isinstance(n_max, numbers.Integral):  # exact at any size, so no int64 wrap-around
        n_max = int(n_max)
    elif math.isfinite(n_max):
        n_max = math.ceil(n_max)
    else:
        raise ValueError(_ERRORS[10][1])
    if n_max > MAX_SUBINTERVALS:
        raise ValueError(_ERRORS[10][1])
    return max(n_max, 0)


def middle_contracts_sabr(F, K_min, K_max, tau, alpha, rho, nu, usd_base, n=None, n0=32, rtol=SIMPSON_RTOL,
                          n_max=2**16, smile_F=None, smile_tau=None) -> dict:
    """Middle parts of qef.fx.moments.middle_contracts for many SABR (beta = 1) smiles in one call.

    Arguments broadcast to one dimension, one element per smile; the smile is
    sabr_vol(K, smile_F, smile_tau, alpha, rho, nu), with smile_F = F and
    smile_tau = tau unless given. With ``n`` the rule is applied at n
    subintervals per side (the reference's integrate(n)); otherwise n is
    doubled from ``n0`` as in the reference. Returns arrays: ``values`` (m, 3)
    and, per smile, ``n``, ``rel_change`` and ``converged`` (NaN, 0 and False
    at a fixed n).
    """
    smile_F = F if smile_F is None else smile_F
    smile_tau = tau if smile_tau is None else smile_tau
    usd = np.asarray(usd_base, dtype=bool).astype(np.float64)
    shape, m, arrs, inc = _inputs(F, K_min, K_max, tau, usd, smile_F, smile_tau, alpha, rho, nu)
    if len(shape) > 1:
        raise ValueError("middle_contracts_sabr expects scalars or one-dimensional arrays")
    values, status = np.empty((m, len(CONTRACTS))), np.empty(m, dtype=np.int32)
    lib = _lib()
    if n is not None:
        n = _check_subintervals(n)
        lib.qef_middle_contracts_fixed(m, *map(_p, arrs), _p(inc), n, _p(values), _p(status))
        _raise_first(status)
        return {"values": values, "n": np.full(m, n), "rel_change": np.full(m, np.nan),
                "converged": np.zeros(m, dtype=bool)}
    n0, n_max = _check_subintervals(n0), _check_n_max(n_max)
    n_out, change, conv = np.empty(m, dtype=np.int64), np.empty(m), np.empty(m, dtype=np.int32)
    lib.qef_middle_contracts_adaptive(m, *map(_p, arrs), _p(inc), n0, float(rtol), n_max, _p(values),
                                      _p(n_out), _p(change), _p(conv), _p(status))
    _raise_first(status)
    return {"values": values, "n": n_out, "rel_change": change, "converged": conv.astype(bool)}


def middle_contracts_fixed(F, K_min, K_max, tau, vol_of_strike, usd_base: bool, n):
    """The reference's integrate(n): Simpson middle parts of the contracts k = 1, 2, 3 at n subintervals per side."""
    sF, st, a, r, v = _smile_parameters(vol_of_strike)
    out = middle_contracts_sabr(F, K_min, K_max, tau, a, r, v, usd_base, n=n, smile_F=sF, smile_tau=st)
    return out["values"][0]


def middle_contracts(F, K_min, K_max, tau, vol_of_strike, usd_base: bool, n0=32, rtol=SIMPSON_RTOL, n_max=2**16):
    """Middle parts ∫_{K_min}^F h_k'' P dK + ∫_F^{K_max} h_k'' C dK, k = 1, 2, 3 (qef.fx.moments.middle_contracts).

    ``vol_of_strike`` must be a SabrSmile with beta = 1 or its ``vol`` method.
    Returns the reference's dict: values, n, rel_change and converged.
    """
    if not K_min < F < K_max:
        raise ValueError(_ERRORS[9][1])
    sF, st, a, r, v = _smile_parameters(vol_of_strike)
    out = middle_contracts_sabr(F, K_min, K_max, tau, a, r, v, usd_base, n0=n0, rtol=rtol, n_max=n_max,
                                smile_F=sF, smile_tau=st)
    change = float(out["rel_change"][0])
    return {"values": out["values"][0], "n": int(out["n"][0]), "rel_change": change, "converged": change < rtol}
