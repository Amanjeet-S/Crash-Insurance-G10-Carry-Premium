"""The published portfolio-level series in data/public/fx_carry_portfolio_series/.

These checks need no private data. They check the files against their
manifest, the publication rules of scripts/export_public_series.py (allowed
columns only, none of the columns removed after the derived-data check, at
most six significant figures, sorted New York month-ends, no leg identity or currency
anywhere) and a few key numbers of the paper recomputed
from the files at the paper's printed precision with the project's own
functions. scripts/reproduce_from_public.py recomputes the portfolio-level
results; these tests keep a fast subset in the suite.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
PUB = ROOT / "data" / "public" / "fx_carry_portfolio_series"
sys.path.insert(0, str(ROOT / "scripts"))
import estimate_long_atm  # noqa: E402
import export_public_series  # noqa: E402
from estimate_e1 import REGIME_BREAK, regime_summary  # noqa: E402

from qef.data.panel import ny_month_ends  # noqa: E402
from qef.fx.conventions import G10  # noqa: E402
from qef.stats.hac import mean_and_se  # noqa: E402

E1 = ["C_skew_market_10d", "phi_market_10d"]
WITHHELD = [f"{q}_{r}" for r in ("market_25d", "smile_10d", "smile_25d") for q in ("C_skew", "phi")]
RET = ["U", "H10", "H25", "Hatm", "c1", "c2", "c3"]
COLUMNS = {
    "primary.csv": ["date", "regime", "status", *E1, "FD", *RET, "var_lo", "var_hi", "oskew_lo", "oskew_hi", "sigma_fx"],
    "extended.csv": ["date", "regime", "status", "C_skew", "FD", "phi", *RET, "sigma_fx"],
    "long_atm.csv": ["date", "status", "n_legs", "U", "Hatm"],
}
REMOVED = {"rr_or", "payoff", "premium", *WITHHELD}  # removed after the derived-data check; never to be published
WINDOWS = {"primary.csv": ("2013-05-31", "2026-08-31", 160), "extended.csv": ("2007-01-31", "2013-04-30", 76),
           "long_atm.csv": ("1995-01-31", "2026-07-31", 379)}
STATUSES = {"ok", "missing_spot", "too_few_currencies"}
CODES = [*G10, "USD"]


def within(x, printed: str) -> bool:
    """Whether x rounds to the printed number at its printed precision."""
    dec = len(printed.split(".")[1]) if "." in printed else 0
    return abs(x - float(printed)) <= 0.5 * 10.0 ** (-dec) * (1 + 1e-9)


@pytest.fixture(scope="module")
def manifest():
    path = PUB / "manifest.json"
    assert path.exists(), f"{path} is missing; the public series are part of the repository"
    return json.loads(path.read_text())


@pytest.fixture(scope="module")
def frames():
    return {name: pd.read_csv(PUB / name, parse_dates=["date"]) for name in COLUMNS}


@pytest.fixture(scope="module")
def raw_rows():
    out = {}
    for name in COLUMNS:
        with open(PUB / name, newline="") as f:
            out[name] = list(csv.reader(f))
    return out


def test_columns_are_exactly_the_allowed_ones(frames):
    assert export_public_series.ALLOWED == COLUMNS
    assert REMOVED <= export_public_series.FORBIDDEN
    for name, df in frames.items():
        assert list(df.columns) == COLUMNS[name]
        assert not set(df.columns) & export_public_series.FORBIDDEN


def test_manifest_matches_the_files(manifest, frames):
    assert set(manifest["files"]) == set(COLUMNS)
    for name, meta in manifest["files"].items():
        data = (PUB / name).read_bytes()
        assert hashlib.sha256(data).hexdigest() == meta["sha256"]
        assert meta["bytes"] == len(data)
        df = frames[name]
        assert meta["rows"] == len(df) and meta["columns"] == list(df.columns)
        assert meta["first_date"] == f"{df.date.iloc[0]:%Y-%m-%d}" and meta["last_date"] == f"{df.date.iloc[-1]:%Y-%m-%d}"


def test_rows_and_dates(frames):
    for name, (first, last, n) in WINDOWS.items():
        df = frames[name]
        assert len(df) == n
        assert df.date.is_unique and df.date.is_monotonic_increasing
        # one row per New York month-end in the window, none missing
        assert list(df.date) == list(ny_month_ends(first, last))


def test_numbers_have_at_most_six_significant_figures(raw_rows):
    num = re.compile(r"-?(\d+)(?:\.(\d+))?(?:e[-+]\d+)?")
    for name, rows in raw_rows.items():
        header = rows[0]
        text_cols = {header.index(c) for c in ("date", "regime", "status") if c in header}
        for row in rows[1:]:
            assert len(row) == len(header)
            for i, cell in enumerate(row):
                if i in text_cols or cell == "":
                    continue
                m = num.fullmatch(cell)
                assert m, f"{name}: {cell!r} is not a plain number"
                digits = (m.group(1) + (m.group(2) or "")).lstrip("0")
                assert len(digits) <= 6, f"{name}: {cell} has more than six significant figures"
                assert cell.lower() not in ("nan", "inf", "-inf")


def test_no_leg_identities_or_currency_names(frames, raw_rows):
    for name, rows in raw_rows.items():
        header = rows[0]
        for cell in header + [c for row in rows[1:] for c in row]:
            assert not any(code in cell.upper() for code in CODES), f"{name}: {cell!r} names a currency"
        labels = {row[header.index(c)] for row in rows[1:] for c in ("regime", "status") if c in header}
        assert labels <= STATUSES | {"zero_rate", "hiking"}, f"{name}: unexpected labels {labels}"
        df = frames[name]
        if "regime" in df:
            assert (df.regime == np.where(df.date < REGIME_BREAK, "zero_rate", "hiking")).all()


def test_status_and_missing_values(frames):
    p, e, la = frames["primary.csv"], frames["extended.csv"], frames["long_atm.csv"]
    ret = [*RET, "sigma_fx"]
    assert (p.status == "ok").sum() == 159 and list(p.status[p.status != "ok"]) == ["missing_spot"]
    assert p.loc[p.status != "ok", ret].isna().all().all() and p[E1 + ["FD"]].notna().all().all()
    assert p.loc[p.status == "ok"].drop(columns=["regime", "status"]).notna().all().all()
    assert (e.status == "ok").sum() == 72 and set(e.status) == {"ok", "too_few_currencies"}
    assert e.loc[e.status != "ok"].drop(columns=["date", "regime", "status"]).isna().all().all()
    assert e.loc[e.status == "ok"].drop(columns=["regime", "status"]).notna().all().all()
    assert (la.status == "ok").all() and set(la.n_legs) == {2, 3} and (la.n_legs == 2).sum() == 36
    assert la.notna().all().all()


def test_exact_identities_hold_to_rounding(frames):
    p = frames["primary.csv"]
    ok = p[p.status == "ok"]
    rel = lambda a, b, scale: np.max(np.abs(a - b) / scale)
    assert rel(ok.c3, -ok.C_skew_market_10d, ok.C_skew_market_10d.abs()) < 1e-5
    terms = ok[["U", "c1", "c2", "c3"]]
    assert np.max(np.abs(ok.H10 - terms.sum(axis=1)) / terms.abs().max(axis=1)) < 2e-5
    for c in ("market_10d",):
        assert rel(p[f"phi_{c}"], p[f"C_skew_{c}"] / p.FD, p[f"phi_{c}"].abs()) < 1e-5
    e = frames["extended.csv"]
    ok = e[e.status == "ok"]
    # extended.csv is rounded to four significant figures (README), so its identities hold to about 1e-3
    assert rel(ok.c3, -ok.C_skew, ok.C_skew.abs()) < 1e-5 and rel(ok.phi, ok.C_skew / ok.FD, ok.phi.abs()) < 1e-3
    terms = ok[["U", "c1", "c2", "c3"]]
    assert np.max(np.abs(ok.H10 - terms.sum(axis=1)) / terms.abs().max(axis=1)) < 2e-3


def test_readme_states_theta0_and_the_licence():
    text = (PUB / "README.md").read_text()
    m = re.search(r"theta0_mean_10d = ([0-9.eE+-]+)", text)
    assert m and within(float(m.group(1)), "0.10")  # paper, Section 5.3
    assert "1 October 2026" in text and "\u2014" not in text
    # the README says what is not published and why, including the columns removed after the derived-data check
    unpublished = text[text.index("## What is not published, and why"):text.index("## Reproduction")]
    assert all(c in unpublished for c in REMOVED)


def test_e1_regime_means_of_phi(frames):
    # Paper, Table 1 (Section 5.1): 0.809 (0.185), 0.454 (0.066), difference -0.355, HAC s.e. 0.204.
    p = frames["primary.csv"]
    s = pd.DataFrame({"date": p.date, "status": "ok", "phi": p.phi_market_10d})
    r = regime_summary(s, "phi", B=19)
    assert (r["n_zero_rate"], r["n_hiking"]) == (104, 56)
    for k, v in (("mean_zero_rate", "0.809"), ("se_zero_rate", "0.185"), ("mean_hiking", "0.454"),
                 ("se_hiking", "0.066"), ("difference", "-0.355"), ("se_difference_hac", "0.204")):
        assert within(r[k], v), (k, r[k])


def test_e3_means_and_theta_ub(frames):
    # Paper, Table 4 and Section 5.3: HML^U 17.7 (12.1), HML^H 8.3 (11.9) bp; theta_UB 0.53, 0.49, 0.45.
    p = frames["primary.csv"]
    ok = p[p.status == "ok"]
    u, h = mean_and_se(ok.U.to_numpy()), mean_and_se(ok.H10.to_numpy())
    assert within(1e4 * u["mean"], "17.7") and within(1e4 * u["se"], "12.1")
    assert within(1e4 * h["mean"], "8.3") and within(1e4 * h["se"], "11.9")
    assert within(u["mean"] / u["se"], "1.47")
    for col, v in (("H10", "0.53"), ("H25", "0.49"), ("Hatm", "0.45")):
        assert within(1 - ok[col].mean() / ok.U.mean(), v)


def test_e5_fx_volatility_factor(frames):
    # Paper, Table 5 (Section 5.4): the AR(1) coefficient of the FX-volatility level is 0.76. As in
    # estimate_e5.main, the AR(1) is fitted to sigma_fx over the months with a realised return of both samples.
    p, e = frames["primary.csv"], frames["extended.csv"]
    s = pd.concat([f.loc[f.status == "ok", ["date", "sigma_fx"]] for f in (e, p)]).sort_values("date")
    assert len(s) == 72 + 159 and s.date.is_unique
    x = s.sigma_fx.to_numpy()
    c, rho = np.linalg.lstsq(np.column_stack([np.ones(len(x) - 1), x[:-1]]), x[1:], rcond=None)[0]
    assert within(rho, "0.76") and (x > 0).all()


def test_long_sample_means_and_bounded_set(frames):
    # Paper, Section 5.3: 379 months, HML^U 37.2 (13.2), t 2.83; ATM-hedged 16.3 (7.7); set [0.27, 0.94].
    la = frames["long_atm.csv"]
    U, H = la.U.to_numpy(), la.Hatm.to_numpy()
    u, h, d = mean_and_se(U), mean_and_se(H), mean_and_se(H - U)
    assert len(la) == 379
    assert within(1e4 * u["mean"], "37.2") and within(1e4 * u["se"], "13.2") and within(u["mean"] / u["se"], "2.83")
    assert within(1e4 * h["mean"], "16.3") and within(1e4 * h["se"], "7.7")
    assert within(1e4 * d["mean"], "-20.9") and within(1e4 * d["se"], "7.9")
    kind, lo, hi, _ = estimate_long_atm.theta_set(H, U)
    assert kind == "interval" and within(lo, "0.27") and within(hi, "0.94")
