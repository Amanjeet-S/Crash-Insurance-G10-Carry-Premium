"""Calibrate month-end SABR smiles and run the butterfly-convention diagnostic.

Stage 1 check 7 and the Stage 2 panel calibration (research design, sections
3 and 11). For every currency and month-end of the primary sample, a SABR
smile (beta = 1) is calibrated to the ATM volatility, 25-delta risk reversal
and 25-delta butterfly under each butterfly reading. The calibrated smile then
predicts the 10-delta risk reversal and butterfly, which were not used in the
fit. The reading that predicts the held-out quotes better is supporting
evidence only. Each smile is also checked for static arbitrage (R2).

Outputs are LSEG-derived and are written to data/private/results/.
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from qef.data.panel import ny_month_ends
from qef.data.smile_inputs import build_month_end_inputs
from qef.fx.arbitrage import check_smile
from qef.fx.conventions import G10
from qef.fx.smile import SmileQuotes, calibrate_sabr, quotes_from_smile

ROOT = Path(__file__).resolve().parents[1]
READINGS = ("market", "smile")


def _calibrate_currency(rows: pd.DataFrame) -> list[dict]:
    out = []
    warm = {r: None for r in READINGS}
    for _, row in rows.sort_values("date").iterrows():
        base = {"currency": row.currency, "date": row.date}
        if not row.complete_calib:
            out += [{**base, "reading": r, "status": "missing_input"} for r in READINGS]
            continue
        conv = G10[row.currency].delta
        q = SmileQuotes(row.F, row.tau, row.atm, row.rr25, row.bf25, 0.25, row.df_base, conv)
        for reading in READINGS:
            rec = {**base, "reading": reading}
            t0 = time.perf_counter()
            try:
                res = calibrate_sabr(q, reading, x0=warm[reading])
            except Exception as exc:  # recorded, not dropped
                out.append({**rec, "status": f"error: {type(exc).__name__}"})
                continue
            sm = res.smile
            rec.update(status="ok" if res.success else "not_converged", alpha=sm.alpha, rho=sm.rho, nu=sm.nu,
                       max_abs_residual=float(np.max(np.abs(res.residuals))), jacobian_cond=res.jacobian_cond,
                       seconds=time.perf_counter() - t0, **res.strikes)
            if res.success:
                warm[reading] = np.array([np.log(sm.alpha), np.arctanh(sm.rho), np.log(sm.nu)])
            sd = row.atm * np.sqrt(row.tau)
            for label, width in (("inner", 4.0), ("wide", 10.0)):
                rep = check_smile(sm.vol, row.F, row.tau, -width * sd, width * sd, n=2001)
                rec[f"arb_ok_{label}"] = rep.convex_ok and rep.slope_bounds_ok and rep.density_ok
                rec[f"min_g_{label}"] = rep.min_g
            if not row.has_10d:
                out.append(rec)
                continue
            try:
                _, rr10, bf10 = quotes_from_smile(sm, row.tau, conv, 0.10, row.df_base, reading)
                rec.update(rr10_pred=rr10, bf10_pred=bf10, rr10_quote=row.rr10, bf10_quote=row.bf10,
                           rr10_err=rr10 - row.rr10, bf10_err=bf10 - row.bf10)
            except Exception as exc:
                rec["pred_status"] = f"error: {type(exc).__name__}"
            print(f"{row.currency} {row.date:%Y-%m} {reading} {rec['status']} {rec.get('seconds', float('nan')):.2f}s", file=sys.stderr, flush=True)
            out.append(rec)
    return out


def _table(frame) -> str:
    return "```\n" + frame.to_string() + "\n```"


def summarise(panel: pd.DataFrame) -> str:
    lines = ["# Smile calibration and butterfly-convention diagnostic (restricted)", ""]
    n = panel.groupby("reading")["status"].value_counts().unstack(fill_value=0)
    lines += ["## Calibration status", "", _table(n), ""]
    ok = panel[panel.status == "ok"]
    agg = ok.groupby("reading").agg(max_resid=("max_abs_residual", "max"), cond_median=("jacobian_cond", "median"),
                                    cond_p99=("jacobian_cond", lambda x: x.quantile(0.99)),
                                    arb_inner=("arb_ok_inner", "mean"), arb_wide=("arb_ok_wide", "mean"),
                                    secs=("seconds", "median"))
    lines += ["## Fit quality and static arbitrage (share arbitrage-free)", "", _table(agg), ""]
    err = ok.assign(rr10_abs=100 * ok.rr10_err.abs(), bf10_abs=100 * ok.bf10_err.abs())
    tab = err.pivot_table(index="currency", columns="reading", values=["rr10_abs", "bf10_abs"], aggfunc="median")
    lines += ["## 10-delta prediction, median absolute error (vol points)", "", _table(tab.round(3)), ""]
    wide = err.pivot_table(index=["currency", "date"], columns="reading", values=["rr10_abs", "bf10_abs"])
    comp = []
    for q in ("rr10_abs", "bf10_abs"):
        d = (wide[(q, "market")] - wide[(q, "smile")]).dropna()
        g = d.groupby(level="currency")
        comp.append(pd.DataFrame({f"{q}: share market better": g.apply(lambda x: (x < 0).mean()),
                                  f"{q}: median diff (mkt - smile)": g.median()}))
        comp.append(pd.DataFrame({f"{q}: share market better": [(d < 0).mean()],
                                  f"{q}: median diff (mkt - smile)": [d.median()]}, index=["ALL"]))
    c = pd.concat([comp[0], comp[1]]).join(pd.concat([comp[2], comp[3]]))
    lines += ["## Paired comparison across month-ends", "", _table(c.round(4)), ""]
    return "\n".join(lines)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--start", default="2013-05-31")
    p.add_argument("--end", default="2026-08-31")
    p.add_argument("--workers", type=int, default=9)
    p.add_argument("--contributor", default="", help="volatility contributor suffix, e.g. FN for Fenics; empty for composite")
    args = p.parse_args()
    raw = ROOT / "data" / "private" / "lseg" / args.retrieval_date / "raw"
    out_dir = ROOT / "data" / "private" / "results" / args.retrieval_date
    out_dir.mkdir(parents=True, exist_ok=True)
    month_ends = ny_month_ends(args.start, args.end)
    inputs = build_month_end_inputs(raw, month_ends, contributor=args.contributor)
    tag = f"_{args.contributor.lower()}" if args.contributor else ""
    inputs.to_csv(out_dir / f"smile_inputs{tag}.csv", index=False)
    groups = [g for _, g in inputs.groupby("currency")]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        records = [r for part in ex.map(_calibrate_currency, groups) for r in part]
    panel = pd.DataFrame(records)
    panel.to_csv(out_dir / f"smile_panel{tag}.csv", index=False)
    (out_dir / f"smile_diagnostic{tag}.md").write_text(summarise(panel))
    print(f"{len(panel)} calibrations written to {out_dir}")


if __name__ == "__main__":
    main()
