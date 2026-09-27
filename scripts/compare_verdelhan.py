"""Sign and magnitude check of HML_FX against the portfolios of Lustig, Roussanov and Verdelhan.

    .venv/bin/python scripts/acquire_verdelhan.py --retrieval-date 2026-09-27
    .venv/bin/python scripts/compare_verdelhan.py

Research design, section 10: "HML_FX is compared in sign and magnitude with
Verdelhan's public portfolios." The public file (CurrencyPortfolios.xls, notes
sheet) holds monthly excess returns in levels, gross of transaction costs, of
currency portfolios sorted on the forward discount as in Lustig, Roussanov and
Verdelhan (2011); its "Developed currencies" sheet has five portfolios and
HML = P5 − P1, to May 2021.

Rules fixed before computing the comparison (27 September 2026):
- the comparison series is the gross developed-currency HML (the G10 panel is a
  subset of developed currencies; the all-currency sheet adds emerging markets);
- HML^U of the primary sample (stage4_months.csv, status ok) dated at month-end
  t, which runs from t to the option expiry about one month later, is paired
  with their return dated at the next month-end;
- the overlap is every such pair with their date no later than May 2021;
- reported: both means in bp per month with Newey–West standard errors, the
  Newey–West standard error of the mean difference, the correlation, the share
  of months with the same sign, and the slope of HML^U on their HML.
The check is descriptive, as the design specifies: the two portfolios differ in
currencies (nine against USD here, a wider developed set there), in the number
of portfolios and in timing, so they should agree in sign and in order of
magnitude, not exactly. Output: verdelhan_comparison.md in the private results
folder; it contains pooled statistics only.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from qef.stats.hac import mean_and_se, ols_hac

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "private" / "verdelhan" / "2026-09-27" / "CurrencyPortfolios.xls"


def main():
    d = ROOT / "data" / "private" / "results" / "2026-09-23"
    lrv = pd.read_excel(SOURCE, sheet_name="Developed currencies")
    lrv = lrv[pd.to_datetime(lrv["Dates"], errors="coerce").notna()][["Dates", "HML = P5 - P1"]]
    lrv = lrv.rename(columns={"Dates": "lrv_date", "HML = P5 - P1": "lrv_hml"}).astype({"lrv_hml": float})
    lrv["key"] = pd.to_datetime(lrv.lrv_date).dt.to_period("M")
    s4 = pd.read_csv(d / "stage4_months.csv", parse_dates=["date"])
    s4 = s4[(s4["sample"] == "primary") & (s4.status == "ok")][["date", "U"]].copy()
    s4["key"] = s4.date.dt.to_period("M") + 1  # the return from month-end t is paired with their next month-end
    m = s4.merge(lrv, on="key", how="inner").sort_values("date")
    u, h = m.U.to_numpy(), m.lrv_hml.to_numpy()
    mu, mh, md = mean_and_se(u), mean_and_se(h), mean_and_se(u - h)
    slope = ols_hac(u, np.column_stack([np.ones(len(h)), h]))
    manifest = json.loads((SOURCE.parent / "manifest.json").read_text())
    lines = ["# HML_FX against the developed-currency HML of Lustig, Roussanov and Verdelhan (restricted summary)", "",
             f"Source file sha256 {manifest['sha256'][:12]}, last modified {manifest['last_modified']}. "
             f"Overlap: {len(m)} months, their dates {m.lrv_date.min():%Y-%m} to {m.lrv_date.max():%Y-%m}.", "",
             "| Quantity | Value |", "| --- | --- |",
             f"| Mean HML^U (bp per month, Newey–West s.e.) | {1e4 * mu['mean']:.1f} ({1e4 * mu['se']:.1f}) |",
             f"| Mean developed HML (bp per month, Newey–West s.e.) | {1e4 * mh['mean']:.1f} ({1e4 * mh['se']:.1f}) |",
             f"| Mean difference (bp per month, Newey–West s.e.) | {1e4 * md['mean']:.1f} ({1e4 * md['se']:.1f}) |",
             f"| Correlation of monthly returns | {np.corrcoef(u, h)[0, 1]:.3f} |",
             f"| Share of months with the same sign | {(np.sign(u) == np.sign(h)).mean():.1%} |",
             f"| Slope of HML^U on their HML (Newey–West s.e.) | {slope['beta'][1]:.3f} ({slope['se'][1]:.3f}) |",
             f"| Standard deviation, HML^U and theirs (bp per month) | {1e4 * u.std(ddof=1):.1f}, {1e4 * h.std(ddof=1):.1f} |", ""]
    (d / "verdelhan_comparison.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
