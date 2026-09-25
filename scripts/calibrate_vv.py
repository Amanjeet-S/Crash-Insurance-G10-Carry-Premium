"""Calibrate month-end vanna–volga smiles and compare their 10-delta predictions with SABR.

Robustness variant 9 of the research design (section 8), with the choices
recorded in the research log (24 September 2026). For every currency and
month-end of the primary sample (composite one-month quotes with the complete
calibration set, from 31 May 2013), a vanna–volga smile is built under each
butterfly reading with each method of qef.fx.vannavolga; the primary method is
the price, eq. (7) of Castagna and Mercurio (2007). Under the market reading
each method is calibrated with its own smile (Reiswich, 2010, Table 3.2).

For each smile the script records the calibration status, the smile strangle,
the pillars and the largest quote residual; the 10-delta put and call strikes
on that smile (smile_delta_strike, which tries qef.fx.gk.strike_from_delta_smile
first and records whether it alone found the strike); whether each method,
with the same pillars, is defined at those strikes; static-arbitrage flags
(result R2) and the share of an ATM-centred strike grid on which the smile is
undefined; and the 10-delta risk reversal and butterfly the smile implies,
which were not used in the fit, as scripts/calibrate_smiles.py does for SABR.
Months where a method is undefined at a 10-delta strike are counted, not
dropped.

Inputs: data/private/results/<retrieval>/smile_inputs.csv, written by
scripts/calibrate_smiles.py. Outputs are LSEG-derived and are written to the
same folder: smile_panel_vv.csv (one row per currency, month-end, reading and
method) and smile_diagnostic_vv.md (restricted aggregate summary, compared
with the SABR panel smile_panel.csv).
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from qef.fx.arbitrage import check_smile
from qef.fx.conventions import G10
from qef.fx.gk import CALL, PUT
from qef.fx.smile import SmileQuotes
from qef.fx.vannavolga import METHODS, PRIMARY_METHOD, calibrate_vv, implied_quotes, smile_delta_strike

ROOT = Path(__file__).resolve().parents[1]
READINGS = ("market", "smile")
WIDTHS = (("inner", 4.0), ("wide", 10.0))  # ± ATM standard deviations of the arbitrage grid, as for SABR


def _smile_record(cal, row, conv, reading) -> dict:
    """Diagnostics of one calibrated smile (everything after the calibration itself)."""
    sm = cal.smile
    rec = {}
    sd = row.atm * np.sqrt(row.tau)
    for label, width in WIDTHS:
        with np.errstate(all="ignore"):
            rep = check_smile(sm.vol, row.F, row.tau, -width * sd, width * sd, n=2001)
            grid = row.F * np.exp(np.linspace(-width * sd, width * sd, 2001))
            rec[f"undefined_share_{label}"] = float(np.mean(np.isnan(sm.vol(grid))))
        rec[f"arb_ok_{label}"] = rep.convex_ok and rep.slope_bounds_ok and rep.density_ok
        rec[f"min_g_{label}"] = rep.min_g
    iq = implied_quotes(sm, conv, 0.10, row.df_base, reading)
    kp, kc = iq["K_put"], iq["K_call"]
    rec.update(K_put10=kp, K_call10=kc)
    for side, phi, target, k in (("put", PUT, -0.10, kp), ("call", CALL, 0.10, kc)):
        rec[f"{side}10_found_gk"] = bool(np.isfinite(smile_delta_strike(sm, target, phi, conv, row.df_base, fallback=False)))
        for m in METHODS:
            rec[f"defined_{side}10_{m}"] = bool(np.isfinite(k) and np.isfinite(sm.vol(k, m)))
    rec.update(rr10_pred=iq["rr"], bf10_pred=iq["bf"], rr10_quote=row.rr10, bf10_quote=row.bf10,
               rr10_err=iq["rr"] - row.rr10, bf10_err=iq["bf"] - row.bf10)
    return rec


def _calibrate_currency(rows: pd.DataFrame) -> list[dict]:
    out = []
    for _, row in rows.sort_values("date").iterrows():
        conv = G10[row.currency].delta
        q = SmileQuotes(row.F, row.tau, row.atm, row.rr25, row.bf25, 0.25, row.df_base, conv)
        for reading in READINGS:
            for method in METHODS:
                rec = {"currency": row.currency, "date": row.date, "reading": reading, "method": method,
                       "bf25_quote": row.bf25}
                t0 = time.perf_counter()
                try:
                    cal = calibrate_vv(q, reading, method)
                except Exception as exc:  # recorded, not dropped
                    out.append({**rec, "status": f"error: {type(exc).__name__}"})
                    continue
                rec.update(status=cal.status, bf_ss=cal.bf_ss, residual=float(np.max(np.abs(cal.residuals))),
                           res_atm=cal.residuals[0], res_rr=cal.residuals[1], res_bf=cal.residuals[2],
                           K_put_ms=cal.strikes.get("K_put_ms", np.nan), K_call_ms=cal.strikes.get("K_call_ms", np.nan))
                if cal.smile is not None:
                    sm = cal.smile
                    rec.update(K1=sm.strikes[0], K2=sm.strikes[1], K3=sm.strikes[2],
                               sigma1=sm.vols[0], sigma2=sm.vols[1], sigma3=sm.vols[2])
                    try:
                        rec.update(_smile_record(cal, row, conv, reading))
                    except Exception as exc:
                        rec["diag_status"] = f"error: {type(exc).__name__}"
                rec["seconds"] = time.perf_counter() - t0
                out.append(rec)
    print(f"{rows.currency.iloc[0]}: {len(out)} smiles", file=sys.stderr, flush=True)
    return out


def _table(frame) -> str:
    return "```\n" + frame.to_string() + "\n```"


def _both_found(frame):
    return np.isfinite(frame.K_put10) & np.isfinite(frame.K_call10)


def summarise(panel: pd.DataFrame, sabr: pd.DataFrame, start: str) -> str:
    n_cm = panel[["currency", "date"]].drop_duplicates().shape[0]
    lines = ["# Vanna–volga smiles: calibration and 10-delta diagnostic (restricted)", "",
             f"Composite one-month quotes with the complete calibration set, month-ends from {start}: {n_cm} "
             f"currency-months, {panel.date.nunique()} month-ends. Primary method: {PRIMARY_METHOD}. "
             "Volatility differences are in vol points.", ""]
    st = panel.groupby(["reading", "method"])["status"].value_counts().unstack(fill_value=0)
    lines += ["## Calibration status", "", _table(st), ""]

    cal = panel[panel.K1.notna()].copy()
    flags = [c for c in cal if c.startswith(("defined_", "arb_ok_")) or c.endswith("_found_gk")]
    # Object dtype when a flag is missing (diagnostics raised) or read back from CSV; a missing flag is a
    # failure, whereas astype(bool) would count NaN as True.
    cal[flags] = cal[flags].eq(True)
    fit = cal.groupby(["reading", "method"]).agg(
        max_residual=("residual", "max"), secs_median=("seconds", "median"),
        arb_free_inner=("arb_ok_inner", "mean"), arb_free_wide=("arb_ok_wide", "mean"),
        undefined_inner=("undefined_share_inner", "mean"), undefined_wide=("undefined_share_wide", "mean"))
    lines += ["## Fit and static arbitrage",
              "", "Smiles with pillars. max_residual in volatility units; arb_free: share passing R2 on ±4 and ±10 "
              "ATM standard deviations; undefined: mean share of that grid where the smile is undefined.", "",
              _table(fit), ""]

    mk = cal[cal.reading == "market"].assign(gap=lambda d: 100 * (d.bf_ss - d.bf25_quote))
    g = mk.groupby("method")["gap"].agg(mean="mean", median="median", max_abs=lambda x: x.abs().max())
    lines += ["## Market reading: smile strangle minus quoted butterfly", "", _table(g.round(4)), ""]

    cal["both10"] = _both_found(cal)
    cal["gk_only"] = cal.put10_found_gk.astype(bool) & cal.call10_found_gk.astype(bool)
    strikes = cal.groupby(["reading", "method"]).agg(
        n=("both10", "size"), put10_found=("K_put10", lambda x: np.isfinite(x).mean()),
        call10_found=("K_call10", lambda x: np.isfinite(x).mean()), both_found=("both10", "mean"),
        currency_months_missing=("both10", lambda x: int((~x).sum())), found_by_gk_alone=("gk_only", "mean"))
    months = cal[~cal.both10].groupby(["reading", "method"]).date.nunique().rename("month_ends_missing")
    strikes = strikes.join(months).fillna({"month_ends_missing": 0})
    lines += ["## 10-delta strikes on each smile",
              "", "Shares of smiles with pillars. A strike is missing where the method is undefined at it or no strike "
              "attains the delta on the defined range around K2. found_by_gk_alone: both strikes found by "
              "strike_from_delta_smile without the fallback scan.", "", _table(strikes.round(4)), ""]

    prim = cal[cal.method == PRIMARY_METHOD]
    cross = prim.groupby("reading")[[f"defined_{s}10_{m}" for s in ("put", "call") for m in METHODS]].mean()
    lines += [f"## Methods defined at the 10-delta strikes of the {PRIMARY_METHOD} smile (same pillars)", "",
              _table(cross.T.round(4)), ""]
    by_ccy = prim.assign(missing=~prim.both10).pivot_table(index="currency", columns="reading", values="missing", aggfunc="sum")
    lines += [f"## Currency-months with a missing 10-delta strike, {PRIMARY_METHOD} smile", "", _table(by_ccy), ""]

    # Held-out comparison on calibrated smiles only, for vanna–volga as for SABR (status ok).
    sab = sabr[(sabr.status == "ok") & (sabr.date >= pd.Timestamp(start))]
    sab = sab.merge(panel[["currency", "date"]].drop_duplicates(), on=["currency", "date"])
    vv = cal[cal.status == "ok"]
    both = pd.concat([vv[["currency", "date", "reading", "method", "rr10_err", "bf10_err"]],
                      sab[["currency", "date", "reading", "rr10_err", "bf10_err"]].assign(method="sabr")])
    both[["rr10_abs", "bf10_abs"]] = 100 * both[["rr10_err", "bf10_err"]].abs()
    both["ok"] = np.isfinite(both.rr10_abs) & np.isfinite(both.bf10_abs)
    models = list(METHODS) + ["sabr"]
    own = both[both.ok].groupby(["reading", "method"]).agg(n=("ok", "size"), rr10_mae=("rr10_abs", "mean"),
                                                           bf10_mae=("bf10_abs", "mean"))
    wide = both[both.ok].pivot_table(index=["currency", "date", "reading"], columns="method", values=["rr10_abs", "bf10_abs"])
    common = wide.dropna()
    com = pd.DataFrame({(q, m): common[(q, m)].groupby(level="reading").mean() for q in ("rr10_abs", "bf10_abs") for m in models})
    com["n"] = common.groupby(level="reading").size()
    lines += ["## Held-out 10-delta quotes: mean absolute error",
              "", "Smiles with status ok. Each model on the currency-months where it predicts both quotes:", "",
              _table(own.round(4)), "",
              "On the currency-months where all four models predict both quotes:", "", _table(com.round(4).T), ""]
    pair = wide[[(q, m) for q in ("rr10_abs", "bf10_abs") for m in (PRIMARY_METHOD, "sabr")]].dropna()
    rows = {}
    for reading, g in pair.groupby(level="reading"):
        rows[reading] = {"n": len(g)}
        for q in ("rr10_abs", "bf10_abs"):
            rows[reading][f"{q}_mae_vv"] = g[(q, PRIMARY_METHOD)].mean()
            rows[reading][f"{q}_mae_sabr"] = g[(q, "sabr")].mean()
            rows[reading][f"{q}_share_vv_better"] = (g[(q, PRIMARY_METHOD)] < g[(q, "sabr")]).mean()
    lines += [f"Paired, {PRIMARY_METHOD} against SABR on their common currency-months:", "",
              _table(pd.DataFrame(rows).T.round(4)), ""]
    ccy = pair.groupby(level=["currency", "reading"]).mean()
    ccy.columns = [f"{q}_{m}" for q, m in ccy.columns]
    lines += [f"By currency, {PRIMARY_METHOD} against SABR (same sample):", "", _table(ccy.unstack("reading").round(3)), ""]
    return "\n".join(lines)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--start", default="2013-05-31")
    p.add_argument("--workers", type=int, default=9)
    p.add_argument("--summary-only", action="store_true", help="rebuild the summary from an existing smile_panel_vv.csv")
    args = p.parse_args()
    out_dir = ROOT / "data" / "private" / "results" / args.retrieval_date
    if args.summary_only:
        panel = pd.read_csv(out_dir / "smile_panel_vv.csv", parse_dates=["date"])
    else:
        inputs = pd.read_csv(out_dir / "smile_inputs.csv", parse_dates=["date"])
        inputs = inputs[inputs.complete_calib & (inputs.date >= pd.Timestamp(args.start))]
        groups = [g for _, g in inputs.groupby("currency")]
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            records = [r for part in ex.map(_calibrate_currency, groups) for r in part]
        panel = pd.DataFrame(records)
        panel.to_csv(out_dir / "smile_panel_vv.csv", index=False)
    sabr = pd.read_csv(out_dir / "smile_panel.csv", parse_dates=["date"])
    (out_dir / "smile_diagnostic_vv.md").write_text(summarise(panel, sabr, args.start))
    print(f"{len(panel)} smiles written to {out_dir}")


if __name__ == "__main__":
    main()
