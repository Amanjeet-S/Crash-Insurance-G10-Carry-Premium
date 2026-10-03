"""Export the portfolio-level month series of the crash-insurance project to data/public/.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/export_public_series.py

Reads my private results of the LSEG retrieval of 23 September 2026
(data/private/results/2026-09-23/) and writes, into
data/public/fx_carry_portfolio_series/:

- primary.csv: the primary sample (composite quotes, three long and three short
  legs), formation month-ends May 2013 to August 2026, from e1_series.csv,
  stage4_months.csv, moments_1m.csv (through
  estimate_moments.portfolio_predictors) and e5_series.csv;
- extended.csv: the extended sample (Fenics quotes, two and two legs),
  January 2007 to April 2013, from e1_series.csv, stage4_months.csv and
  e5_series.csv;
- long_atm.csv: the long at-the-money sample (composite quotes), January 1995
  to July 2026, from long_atm_months.csv (the return columns U and Hatm only);
- README.md: what each file and column is, with units, definitions and the code
  that produces them;
- manifest.json: SHA-256 hash, size, number of rows, columns and first and last
  date of each CSV file.

Every published value is a portfolio-level quantity: an average over the legs
of a carry portfolio, or a ratio or sum of such averages. The publication set
is ALLOWED below, which I fixed after checking a first version against the
derived-data test:
that version also held rr_or, the oriented 10-delta risk reversal, which is a
plain signed average of six raw risk-reversal quotes in quote units, and the
payoff and premium of the long sample's ATM hedge, whose payoff depends on a
single leg in the months in which only one leg's option ends in the money.
Both are now in FORBIDDEN. The script refuses (by assertion) to write any
column outside ALLOWED or inside FORBIDDEN, any status other than those in
STATUSES, or any string that names a currency. It rounds every number to six
significant figures, sorts the rows by date and writes missing values as empty
fields. Before rounding it checks the exact identities stated in the README and
that the published sigma_fx of both samples determines the E5 factor of
estimate_e5.py.
The licence holder confirmed on 1 October 2026 that manipulated or transformed
data may be published and that raw data may not (research log, 1 October 2026).
The README's section on why the series cannot be traced back to the quotes
reports the result of scripts/check_public_reversibility.py.

Running the script twice on the same private results writes byte-identical
files. It reads no raw quote, calls no data provider and changes nothing under
data/private/.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from estimate_e1 import PRIMARY_START, REGIME_BREAK  # noqa: E402
from estimate_moments import portfolio_predictors  # noqa: E402

from qef.fx.conventions import G10  # noqa: E402

RETRIEVAL = "2026-09-23"
RETRIEVAL_TEXT = "23 September 2026"
OUT = ROOT / "data" / "public" / "fx_carry_portfolio_series"
SIG = 6  # significant figures of the published numbers
SIG_FILE = {"primary.csv": 6, "extended.csv": 4, "long_atm.csv": 6}  # Fenics quotes sit on a coarse grid (README)
COMBOS = (("market", 0.10), ("market", 0.25), ("smile", 0.10), ("smile", 0.25))


def combo_name(reading: str, delta: float) -> str:
    return f"{reading}_{round(100 * delta)}d"


PUBLISHED_COMBOS = COMBOS[:1]  # market reading, 10-delta (E1); the other three are withheld (README)
E1_COLUMNS = [c for r, d in PUBLISHED_COMBOS for c in (f"C_skew_{combo_name(r, d)}", f"phi_{combo_name(r, d)}")]
WITHHELD_E1 = [c for r, d in COMBOS[1:] for c in (f"C_skew_{combo_name(r, d)}", f"phi_{combo_name(r, d)}")]
RETURNS = ["U", "H10", "H25", "Hatm", "c1", "c2", "c3"]  # the Stage 4 return columns, empty without a realised return
ALLOWED = {
    "primary.csv": ["date", "regime", "status", *E1_COLUMNS, "FD", *RETURNS, "var_lo", "var_hi", "oskew_lo",
                    "oskew_hi", "sigma_fx"],
    "extended.csv": ["date", "regime", "status", "C_skew", "FD", "phi", *RETURNS, "sigma_fx"],
    "long_atm.csv": ["date", "status", "n_legs", "U", "Hatm"],
}
# Names that must never appear as a column: leg identities, per-currency and audit fields, parameters, and the
# columns removed after the derived-data check (rr_or: signed average of raw quotes; payoff and premium of the
# long sample's ATM hedge: the payoff can depend on a single leg).
FORBIDDEN = {"longs", "shorts", "currency", "currencies", "n_ccy", "n_not_converged", "theta0", "theta0_pa",
             "spot_sub", "quote_sub", "no_rate", "stale_atm", "vol_ok", "alpha", "rho", "nu", "r_vix",
             "vx_contract_expiry", "d_sigma_fx", "end", "sample", "reading", "delta", "rr_or", "payoff", "premium",
             *WITHHELD_E1}
STATUSES = {"ok", "missing_spot", "too_few_currencies"}
REGIMES = {"zero_rate", "hiking"}
WINDOWS = {  # first and last formation month-end of each file
    "primary.csv": ("2013-05-31", "2026-08-31"),
    "extended.csv": ("2007-01-31", "2013-04-30"),
    "long_atm.csv": ("1995-01-31", "2026-07-31"),
}
IDENTITY_TOL = 1e-12  # c3 = -C_skew (market, 10-delta) and H10 = U + c1 + c2 + c3, before rounding
# Months of long_atm.csv in which exactly one leg's ATM option ends in the money: the reason payoff and premium
# are not published. long_atm_months.csv holds no per-leg payoff; I recomputed the count leg by leg from the
# private inputs of estimate_long_atm.py on 1 October 2026 (379 months; legs and payoffs reproduced to 1e-16).
ONE_LEG_ITM_MONTHS = 42


def regime(dates: pd.Series) -> pd.Series:
    return pd.Series(np.where(dates < REGIME_BREAK, "zero_rate", "hiking"), index=dates.index)


def round_sig(x: float, sig: int = SIG) -> float:
    return float(f"{x:.{sig}g}") if np.isfinite(x) else np.nan


def build_primary(d: Path) -> tuple[pd.DataFrame, dict]:
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    e1 = load("e1_series.csv")
    e1 = e1[e1["sample"] == "primary"]
    assert set(e1.status) == {"ok"}, "every primary E1 month-end is expected to have status ok"
    out, fd, legs = None, {}, {}
    for reading, delta in COMBOS:
        s = e1[(e1.reading == reading) & np.isclose(e1.delta, delta)].sort_values("date").set_index("date")
        name = combo_name(reading, delta)
        fd[name], legs[name] = s.FD, s.longs + "|" + s.shorts
        part = s[["C_skew", "phi"]].rename(columns={"C_skew": f"C_skew_{name}", "phi": f"phi_{name}"})
        out = part if out is None else out.join(part, how="outer")
    base = combo_name(*COMBOS[0])
    for name in fd:  # the legs, hence FD, depend only on the ranking, not on the reading or the hedge delta
        assert fd[name].equals(fd[base]) and legs[name].equals(legs[base]), f"FD or legs differ in {name}"
    out["FD"] = fd[base]
    out = out.reset_index()

    s4 = load("stage4_months.csv")
    s4 = s4[s4["sample"] == "primary"].set_index("date")
    assert set(out.date) == set(s4.index), "Stage 4 and E1 month-ends differ"
    out = out.join(s4[["status", *RETURNS]], on="date")
    ok = out.status == "ok"
    assert np.all(np.abs(out.loc[ok, "c3"] + out.loc[ok, f"C_skew_{base}"]) < IDENTITY_TOL)
    assert np.all(np.abs(out.loc[ok, "H10"] - out.loc[ok, ["U", "c1", "c2", "c3"]].sum(axis=1)) < IDENTITY_TOL)
    assert out.loc[~ok, RETURNS].isna().all().all(), "a month without a realised return carries return columns"
    out = out.drop(columns=WITHHELD_E1)  # computed for the checks above, never written

    mom = load("moments_1m.csv")
    e1_10 = e1[(e1.reading == "market") & np.isclose(e1.delta, 0.10)].sort_values("date")
    pp = portfolio_predictors(e1_10, mom)
    good = (pp.status == "ok") & (pp.status_i == "ok")
    pred = pp.loc[good, ["date", "var_lo_i", "var_hi_i", "oskew_lo_i", "oskew_hi_i"]]
    pred = pred.rename(columns=lambda c: c.removesuffix("_i"))
    out = out.merge(pred, on="date", how="left")

    out = add_sigma_fx(out, load("e5_series.csv"), "primary")
    out["regime"] = regime(out.date)

    pm = s4[s4.status == "ok"]
    info = {"theta0_mean_10d": round_sig(float(pm.theta0.mean())), "theta0_months": int(len(pm))}
    return out, info


def add_sigma_fx(out: pd.DataFrame, e5: pd.DataFrame, sample: str) -> pd.DataFrame:
    """Merge the FX-volatility level of estimate_e5.py, present exactly in the months with a realised return."""
    out = out.merge(e5.loc[e5["sample"] == sample, ["date", "sigma_fx"]], on="date", how="left")
    ok = out.status == "ok"
    assert out.loc[ok, "sigma_fx"].notna().all() and out.loc[~ok, "sigma_fx"].isna().all(), f"{sample}: sigma_fx"
    return out


def build_extended(d: Path) -> pd.DataFrame:
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    e1 = load("e1_series.csv")
    e1 = e1[(e1["sample"] == "extended") & (e1.reading == "market") & np.isclose(e1.delta, 0.10)]
    out = e1[["date", "status", "C_skew", "FD", "phi"]].rename(columns={"status": "status_e1"})
    s4 = load("stage4_months.csv")
    s4 = s4[s4["sample"] == "extended"]
    out = out.merge(s4[["date", "status", *RETURNS]], on="date", how="left")
    assert set(s4.date) <= set(out.date), "a Stage 4 extended month-end has no E1 row"
    out["status"] = np.where(out.status_e1 != "ok", out.status_e1, out.status)
    ok = out.status == "ok"
    assert np.all(np.abs(out.loc[ok, "c3"] + out.loc[ok, "C_skew"]) < IDENTITY_TOL)
    assert np.all(np.abs(out.loc[ok, "H10"] - out.loc[ok, ["U", "c1", "c2", "c3"]].sum(axis=1)) < IDENTITY_TOL)
    out = add_sigma_fx(out, load("e5_series.csv"), "extended")
    out["regime"] = regime(out.date)
    return out.drop(columns="status_e1")


def check_e5_factor(d: Path, primary: pd.DataFrame, extended: pd.DataFrame) -> None:
    """The published sigma_fx of both samples, in date order, gives the E5 factor of estimate_e5.main.

    estimate_e5.main fits one AR(1) to sigma_fx over the months with a realised
    return of the extended and primary samples together; its residual is the
    FX-volatility innovation. I check, before rounding, that the months and the
    residual follow from the columns published here.
    """
    e5 = pd.read_csv(d / "e5_series.csv", parse_dates=["date"]).sort_values("date").reset_index(drop=True)
    pub = pd.concat([f.loc[f.status == "ok", ["date", "sigma_fx"]] for f in (extended, primary)])
    pub = pub.sort_values("date").reset_index(drop=True)
    assert pub.date.equals(e5.date), "the E5 months are not the months with status ok of the two files"
    x = pub.sigma_fx.to_numpy()
    c, rho = np.linalg.lstsq(np.column_stack([np.ones(len(x) - 1), x[:-1]]), x[1:], rcond=None)[0]
    resid = x[1:] - c - rho * x[:-1]
    gap = float(np.max(np.abs(resid - e5.d_sigma_fx.to_numpy()[1:])))
    assert np.isnan(e5.d_sigma_fx.iloc[0]) and gap < IDENTITY_TOL, f"E5 factor differs by {gap}"


def build_long_atm(d: Path) -> pd.DataFrame:
    la = pd.read_csv(d / "long_atm_months.csv", parse_dates=["date"])
    la = la[la.date <= pd.Timestamp(WINDOWS["long_atm.csv"][1])].copy()
    assert set(la.status) <= STATUSES
    ok = la.status == "ok"
    # a check of the private series only: payoff and premium are not published (FORBIDDEN)
    gap = (la.Hatm - (la.U + la.payoff - la.premium))[ok].abs().max()
    assert gap < IDENTITY_TOL, f"Hatm = U + payoff - premium fails by {gap}"
    la["n_legs"] = la.n_legs.astype(int)
    return la.drop(columns=["payoff", "premium"])


def finalise(df: pd.DataFrame, name: str) -> pd.DataFrame:
    """Select the allowed columns, check them and round every number to the file's significant figures."""
    cols = ALLOWED[name]
    missing = [c for c in cols if c not in df.columns]
    assert not missing, f"{name}: columns missing from the source: {missing}"
    out = df[cols].sort_values("date").reset_index(drop=True)
    assert list(out.columns) == cols and set(out.columns) <= set(ALLOWED[name]), f"{name}: column outside the allowed list"
    assert not set(out.columns) & FORBIDDEN, f"{name}: forbidden column"
    assert out.date.is_unique and out.date.is_monotonic_increasing
    first, last = WINDOWS[name]
    assert out.date.iloc[0] == pd.Timestamp(first) and out.date.iloc[-1] == pd.Timestamp(last), f"{name}: window"
    assert set(out.status) <= STATUSES, f"{name}: unexpected status {set(out.status) - STATUSES}"
    if "regime" in out:
        assert set(out.regime) <= REGIMES
    for col in out.columns.drop("date"):
        if out[col].dtype == object:  # strings: no currency code, and only the documented labels
            vals = set(out[col].dropna())
            assert not any(c in v for v in vals for c in [*G10, "USD"]), f"{name}: {col} names a currency"
        elif pd.api.types.is_integer_dtype(out[col]):
            continue
        else:
            out[col] = out[col].astype(float).map(lambda v: round_sig(v, SIG_FILE[name]))
    return out


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False, float_format=f"%.{SIG_FILE[path.name]}g", na_rep="", date_format="%Y-%m-%d", lineterminator="\n")


def manifest(frames: dict, out: Path) -> dict:
    files = {}
    for name, df in frames.items():
        p = out / name
        files[name] = {"sha256": sha256(p), "bytes": p.stat().st_size, "rows": int(len(df)), "columns": list(df.columns),
                       "first_date": f"{df.date.iloc[0]:%Y-%m-%d}", "last_date": f"{df.date.iloc[-1]:%Y-%m-%d}"}
    return {"dataset": "fx_carry_portfolio_series",
            "author": "Amanjeet Singh",
            "licence": "CC-BY-4.0 for my original contribution; no right to the underlying LSEG observations (LICENSING.md)",
            "repository": "https://github.com/Amanjeet-S/Crash-Insurance-G10-Carry-Premium",
            "generator": "scripts/export_public_series.py",
            "source": f"private results of the LSEG retrieval of {RETRIEVAL_TEXT} (data/private/, not published)",
            "rounding": {name: f"{SIG_FILE[name]} significant figures" for name in SIG_FILE},
            "missing_value": "empty field",
            "readme": "README.md (not hashed here)",
            "files": files}


def readme(frames: dict, info: dict) -> str:
    p, e, la = frames["primary.csv"], frames["extended.csv"], frames["long_atm.csv"]
    n = lambda df, s="ok": int((df.status == s).sum())
    zr = lambda df: int((df.regime == "zero_rate").sum())
    rows = lambda df: len(df)
    span = lambda df: f"{df.date.iloc[0]:%Y-%m-%d} to {df.date.iloc[-1]:%Y-%m-%d}"
    th = info["theta0_mean_10d"]
    text = f"""# Portfolio-level month series: crash insurance and the G10 carry premium

Amanjeet Singh. This folder holds the month-level portfolio series behind the portfolio-level results of my paper
"Crash insurance and the G10 carry premium" (`paper/main.tex`).
With them, anyone can rerun the paper's main portfolio-level tests without an LSEG licence: E1 under the
market reading with 10-delta hedges, E2 with phi, the secondary moment predictors, E3 and its split by
regime, E5, the hedge-cost ratio and its confidence sets, the long at-the-money sample and the regime
comparison. `scripts/reproduce_from_public.py` recomputes those results from these files alone and
compares each with the number printed in the paper; the results that need the quotes are listed below.

Every value is a portfolio-level quantity: an average over the legs of a carry portfolio (six legs in
the primary sample, four in the extended sample, four or six in the long at-the-money sample), an
average over the nine currencies (`sigma_fx`), or a ratio or sum of such averages. No value refers to a
single currency, and no file says which currencies were held.

## Licence and provenance

The series are derived from quotes that I obtained from LSEG Workspace under a student licence provided
by my university. They are transformed data. Each value is computed from the quotes of four to nine
currencies at once and averaged over the legs of a portfolio whose composition is not published. Most
values are nonlinear functions of calibrated option smiles, forward rates and realised spot rates (option
premia and payoffs, hedged returns, hedge-cost terms and moment bounds); `FD` is an average of log forward
discounts with signs set by the unpublished legs, and `sigma_fx` an average, over the nine currencies
and the business days of a return window, of absolute daily log spot changes. My university, the
licence holder, confirmed on 1 October 2026 that manipulated or transformed data may be published and
that raw data may not (research log, 1 October 2026). I publish these series under that confirmation.
They contain no raw quote, no per-currency value, no leg identity and no calibrated parameter. The raw
LSEG data may not be redistributed, and anyone who wants to rebuild these series from the quotes needs
their own LSEG Workspace licence.

I release the series under [Creative Commons Attribution 4.0 International](../../../LICENSES/CC-BY-4.0.txt)
as my original research results; the grant covers my contribution to them and confers no right to the
underlying LSEG observations (`LICENSING.md`). For attribution, cite Amanjeet Singh, "Crash insurance
and the G10 carry premium", and this repository.

I wrote the files with `scripts/export_public_series.py` from my private results of the retrieval of
{RETRIEVAL_TEXT}. `manifest.json` gives the SHA-256 hash, the size, the number
of rows, the columns and the first and last date of each CSV file.

## Files

| File | Sample | Quotes | Legs per side | Rows | Formation month-ends | Rows with status ok |
| --- | --- | --- | --- | --- | --- | --- |
| `primary.csv` | primary | composite | 3 | {rows(p)} | {span(p)} | {n(p)} |
| `extended.csv` | extended | Fenics | 2 | {rows(e)} | {span(e)} | {n(e)} |
| `long_atm.csv` | long at-the-money | composite (spot, forward points, ATM volatility) | 2 or 3 | {rows(la)} | {span(la)} | {n(la)} |

All three files are comma-separated with a header row. Numbers are rounded to six significant figures in
`primary.csv` and `long_atm.csv` and to four in `extended.csv` (section on traceability below);
a missing value is an empty field. Rows are sorted by date.

## Samples, regimes and dates

- `date`: the formation month-end t, the last New York business day of the calendar month (Federal
  Reserve holiday calendar), written YYYY-MM-DD. Quotes are sampled at t under the five-business-day
  substitution rule of the research design. A return dated t runs from t to the expiry of the one-month
  option traded at t, about one month later, and is evaluated at the end-of-day spot on that expiry date
  (`estimate_stage4.month_rows` and `estimate_stage4.spot_on`; `qef.data.smile_inputs.option_dates`).
- Primary sample: composite quotes, three long and three short legs, {span(p)}: {rows(p)} month-ends,
  {zr(p)} in the zero-rate regime and {rows(p) - zr(p)} in the hiking regime. The return window of the last
  month-end ends after the last spot observation of the retrieval, so {n(p)} month-ends have realised returns.
- Extended sample: Fenics quotes, two long and two short legs, {span(e)}: {rows(e)} month-ends, of which
  {n(e)} have the four currencies with the calibration set that the sample needs. Fenics quotes before 2010
  update infrequently, butterflies in particular, so this sample is secondary evidence.
- Long at-the-money sample: composite spot, one-month forward points and one-month ATM volatility only,
  {span(la)}: {rows(la)} month-ends with realised returns ({int((la.n_legs == 2).sum())} with two legs per
  side, {int((la.n_legs == 3).sum())} with three). The portfolio holds three and three legs when at least six
  currencies are available and two and two when four or five are (`scripts/estimate_long_atm.py`, rules 1
  to 3). The month-end of August 2026, whose return is not yet realised, is not included.
- `regime`: `zero_rate` for formation month-ends before 1 January 2022 and `hiking` from that date
  (`estimate_e1.REGIME_BREAK`). The labels follow the research design and denote calendar periods, not
  rate paths. Every month-end of the extended sample falls before the break, so its label is `zero_rate`
  throughout; the paper's regime comparisons use the primary sample only.

## Units

Premia, returns, hedge-cost terms, `C_skew` and `FD` are decimal fractions per USD of forward notional,
not basis points: 0.0017 is 17 bp. Premia are carried to the forward's delivery date, so they are in USD
per USD of forward notional at delivery (paper, Section 4.2). `phi` is dimensionless. The variance bounds `var_lo` and `var_hi` are in squared
log-return units (the variance of the one-month log return, not annualised); the oriented skewness bounds
`oskew_lo` and `oskew_hi` are dimensionless. `sigma_fx` is a mean absolute daily log change.

## Columns of primary.csv

Each leg buys a one-month protective option that pays when the leg loses: a put on the currency for a
long leg and a call on the currency for a short leg (for the dollar-base pairs, the corresponding call or
put on USD), with notional equal to the forward notional. Strikes are found from deltas in each pair's
own convention on the calibrated SABR smile (beta = 1), and premia are Garman-Kohlhagen prices in forward
form at the stated volatility (`qef.fx.crash`, `qef.fx.gk`, `qef.fx.sabr`).

| Column | Definition | Code |
| --- | --- | --- |
| `date`, `regime` | As above. | |
| `status` | `ok`: every column is present. `missing_spot`: the spot at the option expiry is not observed in the retrieval (the month-end of August 2026, whose return window ends after the retrieval), so `U`, `H10`, `H25`, `Hatm`, `c1`, `c2`, `c3` and `sigma_fx` are empty; the E1 columns, `FD` and the moment predictors are present. Every primary month-end has E1 status ok. | `estimate_stage4.month_rows` |
| `C_skew_market_10d` | The ex-ante skew price of protection, (1/3) times the sum over the six legs of V_smile minus V_flat, where V_smile is the premium of the leg's protective option at its 10-delta strike on the calibrated smile (market-strangle reading of the 25-delta butterfly, the primary reading of E1), priced at the smile volatility, and V_flat is the premium at the same strike priced at the ATM volatility. | `estimate_e1.e1_series`, `qef.fx.crash.leg_skew_cost` |
| `phi_market_10d` | phi = C_skew / FD, the skew price per unit of forward discount: the E1 estimand. | `estimate_e1.e1_series` |
| `FD` | Forward-discount spread, (1/3) times (the sum of fd over the long legs minus the sum over the short legs), with fd = log(X/F), X the spot and F the one-month outright forward, both as USD per unit of the currency. The legs are the three highest and the three lowest fd at t. | `qef.fx.crash.forward_discount`, `rank_legs` |
| `U` | Unhedged carry return HML^U: (1/3) times the sum over the legs of the signed forward return, rx = X_T/F - 1 for a long leg and minus that for a short leg, with X_T the spot on the option expiry date. | `estimate_stage4.month_rows`, `qef.fx.crash.forward_return` |
| `H10`, `H25`, `Hatm` | Hedged carry return HML^H: U plus (1/3) times the sum over the legs of the option payoff less its premium. `H10` and `H25` use the smile 10-delta and 25-delta strikes and the smile premium; `Hatm` uses the delta-neutral-straddle strike and the premium at the ATM volatility. | `estimate_stage4.month_rows` |
| `c1` | E3 term (i), payoff: (1/3) times the sum over the legs of the 10-delta payoff less the Garman-Kohlhagen premium at the realised volatility, the annualised standard deviation of daily log spot changes over the 21 New York business days ending at t. | `estimate_stage4.month_rows`, `realised_vol` |
| `c2` | E3 term (ii), volatility level: minus (1/3) times the sum over the legs of V_flat less that realised-volatility premium, at the 10-delta strike. | `estimate_stage4.month_rows` |
| `c3` | E3 term (iii), skew: minus (1/3) times the sum over the legs of V_smile less V_flat at the 10-delta strike. | `estimate_stage4.month_rows` |
| `var_lo`, `var_hi` | Lower and upper endpoints of the risk-neutral variance of y = ln(X_T/F) under the USD forward measure, from contracts spanned by the calibrated smile between its 10-delta strikes and bounded beyond them by result R4 under tail setting (i) (boundary elasticities), averaged endpoint by endpoint over the six legs. | `estimate_moments.portfolio_predictors`, `qef.fx.moments.implied_moment_intervals` |
| `oskew_lo`, `oskew_hi` | Endpoints of the oriented skewness under the same setting: the mean over the legs of the skewness interval of a long leg and of minus the skewness interval of a short leg (a short leg's interval is [-skew_hi, -skew_lo]), averaged endpoint by endpoint. | `estimate_moments.portfolio_predictors` |
| `sigma_fx` | Global FX-volatility level over the return window: the mean over the window's business days (after t, up to the one-month expiry of EURUSD) of the cross-sectional mean absolute daily log spot change of the nine currencies against USD (those with a spot on the day). It uses daily spot only. The E5 factor is the residual of an AR(1) fitted by least squares to this series over the months with status `ok` of `extended.csv` and `primary.csv` together, in date order. | `estimate_e5.fx_abs_returns`, `estimate_e5.main` |

## Columns of extended.csv

The same definitions as `primary.csv`, with two long and two short legs (weights 1/2 in place of 1/3),
Fenics smiles and the market-strangle reading at 10 delta only.

| Column | Definition | Code |
| --- | --- | --- |
| `date`, `regime` | As above; `regime` is `zero_rate` in every row. | |
| `status` | `ok`: every column is present. `too_few_currencies`: fewer than four currencies have the calibration set at the month-end ({n(e, 'too_few_currencies')} month-ends), so every other column is empty. | `estimate_e1.e1_series`, `estimate_stage4.month_rows` |
| `C_skew`, `phi`, `FD` | As `C_skew_market_10d`, `phi_market_10d` and `FD` in `primary.csv`. | `estimate_e1.e1_series` |
| `U`, `H10`, `H25`, `Hatm`, `c1`, `c2`, `c3` | As in `primary.csv`. | `estimate_stage4.month_rows` |
| `sigma_fx` | As in `primary.csv`, over the return window of the extended month-end. | `estimate_e5.fx_abs_returns`, `estimate_e5.main` |

## Columns of long_atm.csv

| Column | Definition | Code |
| --- | --- | --- |
| `date`, `status` | As above; `status` is `ok` in every row. | `estimate_long_atm.month_rows` |
| `n_legs` | Legs per side, 2 or 3 (weights 1/n). | `estimate_long_atm.n_legs` |
| `U` | Unhedged carry return HML^U, as in `primary.csv`. | `estimate_long_atm.month_rows` |
| `Hatm` | ATM-hedged return: U plus (1/n) times the sum over the legs of the payoff less the premium of the protective option at the delta-neutral-straddle strike, the premium being the Garman-Kohlhagen premium at the quoted ATM volatility, carried to delivery. | `estimate_long_atm.month_rows` |

## Exact identities

Some columns are exact functions of others and add no information; I publish them for convenience.
Before rounding, `phi` = `C_skew` / `FD` in both files; `c3` = -`C_skew_market_10d` in
`primary.csv` and `c3` = -`C_skew` in `extended.csv` (to 1e-16); and `H10` = `U` + `c1` + `c2` + `c3` in
both files (to 3e-16). After rounding these hold to the rounding error of each file's significant figures.

## The diffusive null of the 10-delta hedge

The reference theta0 of result R6 for the 10-delta hedge is, for each month, the average over the six
legs of the absolute forward delta of the hedging options (`estimate_stage4.month_rows`). The paper uses
only its mean over the {info['theta0_months']} return months of the primary sample, so only the mean is
published:

    theta0_mean_10d = {th:.{SIG}g}

## What is not published, and why

- The skew price and phi of the market-strangle reading at 25 delta and of the smile-strangle reading at
  10 and 25 delta (`C_skew_market_25d`, `phi_market_25d`, `C_skew_smile_10d`, `phi_smile_10d`,
  `C_skew_smile_25d`, `phi_smile_25d`: E4 and the 25-delta comparison of E1). Each is another average over the same legs of a
  function of the same quotes, and together with the other columns they would let an adversary who knew
  the legs and the exact forwards narrow single risk-reversal quotes; without them that adversary cannot
  (section on traceability below).
- The oriented 10-delta risk reversal `rr_or`, the secondary E2 predictor of the paper: it is a plain
  signed average of six raw 10-delta risk-reversal quotes, in quote units, and so the least transformed
  of the portfolio-level series.
- The payoff and the premium of the at-the-money hedge in `long_atm.csv`, which `Hatm` combines: in
  {ONE_LEG_ITM_MONTHS} of the {rows(la)} months exactly one leg's at-the-money option ends in the money, so
  the published payoff would depend on that leg alone. `Hatm` - `U`, the payoff less the premium, always
  mixes all legs, because every leg pays a premium.
- Leg identities (which currencies are long and short at each month-end), currency counts other than
  `n_legs`, and any per-currency quantity: with them a reader could recover or closely approximate
  individual quotes.
- Calibrated smile parameters and the smile inputs: a calibrated smile reproduces the quotes it was fitted
  to, so its parameters are not a substantive transformation.
- theta0 for each month: the paper uses only its mean, given above.
- Counts of substituted quotes and other audit detail: they enter no test.
- The VIX roll-down factor `r_vix`: it comes from Cboe VX futures settlements, whose redistribution terms
  are not established. A reader downloads them with `scripts/acquire_cboe_vx.py`, and
  `scripts/reproduce_from_public.py` recomputes the factor.
- The series of the robustness grid and of the three-month tenor: the base and ATM rows of the
  robustness grid, and the returns and theta_UB of its 25-delta row, can be recomputed from
  `primary.csv`; the other variants need the quotes.
- Verdelhan's currency portfolios: they are downloaded with `scripts/acquire_verdelhan.py`.

## Reproduction

From the root of a clone, after the installation in `docs/development.md`:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/reproduce_from_public.py \\
        --out reproduction_report.md \\
        [--cboe-dir data/private/cboe/<retrieval date>/raw] \\
        [--verdelhan-file data/private/verdelhan/<retrieval date>/CurrencyPortfolios.xls]

By default the report goes to the git-ignored `data/private/results/public_reproduction.md`. The
script checks the hashes in `manifest.json`, recomputes the paper's portfolio-level results with
the project's own estimation functions, compares each number with the paper at its printed precision
(and, where my private results exist, with my private summaries) and exits with a non-zero status if any
number that these files determine does not match. E5 needs the Cboe files, and the comparison with
published portfolios needs Verdelhan's file and the `xlrd` package (pinned in `requirements.lock`); each
is used only if present. `tests/test_public_series.py`
checks the files and a few key numbers in the test suite.

The E5 regressions of both samples are reproduced: their FX-volatility innovation is the residual of an
AR(1) fitted to `sigma_fx` over the months with status `ok` of `extended.csv` and `primary.csv` together,
and the VIX roll-down is recomputed from the Cboe files.

Some results cannot be recomputed from these files: everything that needs a single currency's smile or
quotes (the identification of option-implied moments, the per-currency moment intervals and tail exponents,
the out-of-sample test with 5-delta quotes, the leg-by-leg comparison of contributors, the SOFR and
named-broker checks, the calibration diagnostics and the theta0 of the at-the-money hedge), E4 and the
25-delta rows of E1 (their skew prices are withheld), the robustness variants other than the base,
25-delta and ATM rows (of the 25-delta row only the returns and theta_UB), the three-month tenor, the
stale-butterfly check of the extended sample and the account of 5 August 2024. Every result that uses the
oriented risk reversal is also among them, because `rr_or` is not published: its rows of Table 3 in sample
and out of sample, its persistence and its correlations with `phi` and with the return innovations
(Section 5.2), and its change between the regimes (Section 5.1). So are the mean payoff and mean premium
of the at-the-money hedge in 2008 (Section 5.3), because the long sample's payoff and premium are not
published; the long sample's `U` and `Hatm`, and everything the paper builds from them, are reproduced.
The model validation of the moment code is synthetic and already public
(`scripts/validate_moments_models.py`).
`scripts/reproduce_from_public.py` lists each case and its reason.

To rebuild these files from the quotes, a reader with their own LSEG licence runs `scripts/reproduce.sh`
and then `scripts/export_public_series.py`.

## Why the series cannot be traced back to the quotes

LSEG's redistribution guidance treats an output as derived data if it is unrecognisable, non-reversible
and cannot be traced back to the original content without exceptional effort. I checked these files
against that test with `scripts/check_public_reversibility.py`, which reads my private data and reports
pooled statistics only. Its adversary is stronger than any reader: it knows everything about a month
except the option quotes, including the unpublished leg identities and the exact LSEG spot, forwards and
rates.

- Unrecognisable. No column is a quote, the price of a single option or a value of a vendor series.
  Every value averages at least four currencies (nine for `sigma_fx`) whose identities are not published,
  and most are nonlinear functions of calibrated smiles, forwards and realised spot rates.
- Non-reversible. In every month the published values that depend on the option quotes are fewer than
  the quotes they depend on: 9 for 18 in the primary sample (5 in its last month, whose return is not
  yet realised), 5 for 12 in the extended sample and 1 for 4 or 6 in the long sample. A continuous map from more unknowns to fewer values cannot be one-to-one, so
  infinitely many quote vectors reproduce every published value; the sensitivity matrix of the published
  values has full rank in every month, leaving 9 (13 in that last month), 7 and 3 to 5 directions
  unconstrained.
- Constructed alternatives. For three months of each sample, spread over its span, I constructed quote
  vectors that reproduce every published value of the month at the file's precision while differing
  from the true quotes materially: by at least 0.5 volatility points in every month, by 2 in two of the
  three primary months and by 4 in every extended and long month. Every alternative calibrates, passes
  the project's arbitrage checks and stays within the range of quotes seen in the data.
- Not traceable on the quote grid. Real quotes lie on a grid (0.0005 volatility points for the composite
  quotes, 0.00625 to 0.025 for Fenics). Within a box of side one volatility point centred on the true
  quotes, the quote vectors
  on that grid that reproduce every rounded value of a month number at least about 10^3 in the least
  protected month of the long sample, 10^5.9 in the extended sample and 10^12.9 in the primary sample
  (medians 10^8.8, 10^10.1 and 10^18.4). For the extended sample, whose Fenics quotes sit on the coarsest
  grid, this required rounding `extended.csv` to four significant figures instead of six.
- Best attainable estimate. Given the leg identities, the exact spot and forwards, a prior for each quote
  from public information (the realised volatility for at-the-money volatility, zero for the risk
  reversal) and even the exact sensitivities of the published values at the true quotes, the best linear
  estimate leaves a median uncertainty of 1.3 volatility points for at-the-money volatility and 0.7 for
  the 25-delta risk reversal in the primary sample (2.4 and 1.6 in the extended sample; 2.4 for
  at-the-money volatility in the long sample). No risk reversal is pinned down to within 0.10
  volatility points, and about 1% of at-the-money and risk-reversal quotes to within 0.25. The 25-delta
  butterfly varies so little that public information alone predicts it to about 0.09 volatility points;
  the series reduce that only to 0.08.

Before release I withdrew the columns that failed this review: `rr_or`, the long sample's payoff and
premium, and the skew prices and phi of the market reading at 25 delta and of the smile reading. With
those three skew-price series, the same attacker would have pinned about a fifth of primary risk
reversals to within 0.25 volatility points.
"""
    assert "\u2014" not in text, "no em dashes"
    return text


def main():
    p = argparse.ArgumentParser(description="Export the public portfolio-level series of the crash-insurance project.")
    p.add_argument("--retrieval-date", default=RETRIEVAL)
    p.add_argument("--out", type=Path, default=OUT)
    args = p.parse_args()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    primary, info = build_primary(d)
    extended = build_extended(d)
    check_e5_factor(d, primary, extended)
    frames = {"primary.csv": finalise(primary, "primary.csv"),
              "extended.csv": finalise(extended, "extended.csv"),
              "long_atm.csv": finalise(build_long_atm(d), "long_atm.csv")}
    assert frames["primary.csv"].date.iloc[0] == PRIMARY_START
    args.out.mkdir(parents=True, exist_ok=True)
    for name, df in frames.items():
        write_csv(df, args.out / name)
    (args.out / "manifest.json").write_text(json.dumps(manifest(frames, args.out), indent=2) + "\n")
    (args.out / "README.md").write_text(readme(frames, info))
    for name, df in frames.items():
        print(f"{name}: {len(df)} rows, {len(df.columns)} columns, {df.date.iloc[0]:%Y-%m-%d} to {df.date.iloc[-1]:%Y-%m-%d}")
    print(f"theta0_mean_10d = {info['theta0_mean_10d']:.{SIG}g}; written to {args.out}")


if __name__ == "__main__":
    main()
