"""Cross-check of the pricing and delta conventions against QuantLib.

QuantLib is an independent C++ implementation, used here through its Python
bindings. Tests are skipped if it is not installed.
"""

import numpy as np
import pytest

ql = pytest.importorskip("QuantLib")

from qef.fx.gk import CALL, PUT, DeltaConvention, atm_dns_strike, delta, premium, strike_from_delta
from qef.fx.sabr import sabr_vol

TYPES = {
    "Spot": (ql.DeltaVolQuote.Spot, DeltaConvention(spot=True, premium_adjusted=False)),
    "Fwd": (ql.DeltaVolQuote.Fwd, DeltaConvention(spot=False, premium_adjusted=False)),
    "PaSpot": (ql.DeltaVolQuote.PaSpot, DeltaConvention(spot=True, premium_adjusted=True)),
    "PaFwd": (ql.DeltaVolQuote.PaFwd, DeltaConvention(spot=False, premium_adjusted=True)),
}
CASES = [  # spot, quote-currency DF, base-currency DF, vol, tenor
    (1.10, 0.9965, 0.9983, 0.08, 1 / 12),
    (150.0, 0.9995, 0.9960, 0.11, 1 / 12),
    (0.65, 0.9900, 0.9870, 0.14, 0.25),
    (1.35, 0.9600, 0.9750, 0.25, 2.0),
]
QL_TYPE = {CALL: ql.Option.Call, PUT: ql.Option.Put}


def _calc(phi, dt, S, dq, db, sig, tau):
    return ql.BlackDeltaCalculator(QL_TYPE[phi], dt, S, dq, db, sig * np.sqrt(tau))


@pytest.mark.parametrize("S,dq,db,sig,tau", CASES)
def test_premium_matches_black_formula(S, dq, db, sig, tau):
    F = S * db / dq
    for K in F * np.exp(np.linspace(-0.3, 0.3, 7) * sig * np.sqrt(tau) * 3):
        for phi in (CALL, PUT):
            ref = ql.blackFormula(QL_TYPE[phi], K, F, sig * np.sqrt(tau), dq)
            assert premium(F, K, sig, tau, dq, phi) == pytest.approx(ref, rel=1e-12, abs=1e-14)


@pytest.mark.parametrize("name", TYPES)
@pytest.mark.parametrize("S,dq,db,sig,tau", CASES)
def test_delta_from_strike(name, S, dq, db, sig, tau):
    dt, conv = TYPES[name]
    F = S * db / dq
    for K in F * np.exp(np.linspace(-2, 2, 9) * sig * np.sqrt(tau)):
        for phi in (CALL, PUT):
            ref = _calc(phi, dt, S, dq, db, sig, tau).deltaFromStrike(K)
            assert delta(F, K, sig, tau, phi, conv, db) == pytest.approx(ref, abs=1e-13)


@pytest.mark.parametrize("name", TYPES)
@pytest.mark.parametrize("S,dq,db,sig,tau", CASES)
@pytest.mark.parametrize("target", [0.10, 0.25])
def test_strike_from_delta(name, S, dq, db, sig, tau, target):
    dt, conv = TYPES[name]
    F = S * db / dq
    for phi in (CALL, PUT):
        ref = _calc(phi, dt, S, dq, db, sig, tau).strikeFromDelta(phi * target)
        assert strike_from_delta(phi * target, F, sig, tau, phi, conv, db) == pytest.approx(ref, rel=1e-9)


@pytest.mark.parametrize("name", TYPES)
@pytest.mark.parametrize("S,dq,db,sig,tau", CASES)
def test_atm_delta_neutral_strike(name, S, dq, db, sig, tau):
    dt, conv = TYPES[name]
    F = S * db / dq
    ref = _calc(CALL, dt, S, dq, db, sig, tau).atmStrike(ql.DeltaVolQuote.AtmDeltaNeutral)
    assert atm_dns_strike(F, sig, tau, conv) == pytest.approx(ref, rel=1e-12)


@pytest.mark.parametrize("beta", [1.0, 0.5])
@pytest.mark.parametrize("alpha,rho,nu", [(0.10, -0.3, 1.2), (0.08, 0.4, 2.5), (0.12, 0.0, 0.6)])
def test_sabr_matches_hagan_in_quantlib(beta, alpha, rho, nu):
    F, tau = 1.10, 0.25
    a = alpha * F ** (1 - beta)  # keep ATM vol comparable across beta
    for K in F * np.exp(np.linspace(-0.25, 0.25, 11)):
        ref = ql.sabrVolatility(K, F, tau, a, beta, nu, rho)
        assert float(sabr_vol(K, F, tau, a, rho, nu, beta)) == pytest.approx(ref, rel=1e-10)
