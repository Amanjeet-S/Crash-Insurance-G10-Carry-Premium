"""E5: systemic-risk exposure of hedged and unhedged carry (research design, section 6).

Factors, measured over each carry return window (month-end t to the one-month
option expiry E_t, the dates of the Stage 4 returns):

- VIX roll-down return: minus the percentage change of the second-month VIX
  future's settlement price over the window, holding the same contract. This
  is a month-end analogue of the roll-down strategy of Caballero and Doyle
  (2012, NBER WP 18644), who short the VIX future whose expiry matches the
  one-month forward's maturity and hold it to expiry. Contract prices are
  Cboe daily settlements.
- Global FX-volatility innovation: the mean over the window's days of the
  cross-sectional average absolute daily log spot change of the nine
  currencies against USD (Menkhoff, Sarno, Schmeling and Schrimpf, 2012,
  eq. 4, applied to this currency set and window), and the residual of an
  AR(1) fitted to that series.

Regressions of HML^U and HML^H (10Δ) on the factors use Newey–West errors.
The exposure-adjusted hedge-cost ratio is 1 − α_H/α_U, with a test-inversion
set from regressions of HML^H − (1 − θ)HML^U on the factors. As a data check,
the second-month Cboe settlement at month-ends is compared with the LSEG
second-month continuation (VXc2). Outputs go to data/private/results/.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from qef.data.panel import read_raw
from qef.data.vix import load_vx, price_on
from qef.data.smile_inputs import option_dates
from qef.fx.conventions import G10
from qef.stats.hac import ols_hac

ROOT = Path(__file__).resolve().parents[1]


def vix_rolldown(contracts: dict, t, end) -> tuple[float, pd.Timestamp]:
    live = [exp for exp in contracts if exp > t]
    if len(live) < 2:
        return np.nan, pd.NaT
    exp = live[1]  # second-month contract at t
    if exp < end:
        return np.nan, exp
    s = contracts[exp]
    p0, p1 = price_on(s, t), price_on(s, end)
    return -(p1 / p0 - 1.0), exp


def fx_abs_returns(raw: Path) -> pd.Series:
    """Daily cross-sectional mean of |Δ log spot| over the nine currencies (available ones)."""
    cols = {}
    for c in G10:
        f = read_raw(raw / "spot" / f"{c}=.csv")
        f = f[~f.index.duplicated(keep="last")].sort_index()
        mid = ((f["BID"] + f["ASK"]) / 2).dropna()
        mid = mid[mid.index.dayofweek < 5]
        cols[c] = np.log(mid).diff().abs()
    return pd.DataFrame(cols).mean(axis=1, skipna=True).dropna()


def ratio_set(H, U, F, grid=np.linspace(-10, 10, 4001)):
    X = np.column_stack([np.ones(len(H)), F])
    ok = []
    for th in grid:
        r = ols_hac(H - (1 - th) * U, X)
        ok.append(abs(r["beta"][0]) <= 1.96 * r["se"][0])
    ok = np.array(ok)
    if not ok.any():
        return "empty", np.nan, np.nan
    kind = "unbounded" if ok[0] or ok[-1] else ("interval" if np.all(np.diff(np.where(ok)[0]) == 1) else "union of intervals")
    return kind, float(grid[ok].min()), float(grid[ok].max())


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--cboe-date", default="2026-09-24")
    args = p.parse_args()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    raw = ROOT / "data" / "private" / "lseg" / args.retrieval_date / "raw"
    contracts, rescaled = load_vx(ROOT / "data" / "private" / "cboe" / args.cboe_date / "raw")
    months = pd.read_csv(d / "stage4_months.csv", parse_dates=["date"])
    months = months[months.status == "ok"].sort_values("date").reset_index(drop=True)
    absr = fx_abs_returns(raw)

    rows = []
    for _, m in months.iterrows():
        end = option_dates(m.date, "EUR")[2]
        r_vix, exp = vix_rolldown(contracts, m.date, end)
        window = absr.loc[(absr.index > m.date) & (absr.index <= end)]
        rows.append({"date": m.date, "end": end, "vx_contract_expiry": exp, "r_vix": r_vix,
                     "sigma_fx": float(window.mean()) if len(window) else np.nan})
    f = pd.DataFrame(rows)
    x = f.sigma_fx.to_numpy()
    c, rho = np.linalg.lstsq(np.column_stack([np.ones(len(x) - 1), x[:-1]]), x[1:], rcond=None)[0]
    f["d_sigma_fx"] = np.r_[np.nan, x[1:] - c - rho * x[:-1]]
    df = months.merge(f, on="date")
    df.to_csv(d / "e5_series.csv", index=False)

    out = {}
    for sample in ("primary", "extended"):
        s = df[(df["sample"] == sample)].dropna(subset=["r_vix", "d_sigma_fx"])
        Fm = s[["r_vix", "d_sigma_fx"]].to_numpy()
        X = np.column_stack([np.ones(len(s)), Fm])
        res = {"n": len(s)}
        for name in ("U", "H10"):
            r = ols_hac(s[name].to_numpy(), X)
            for k, lab in enumerate(("alpha_bp", "beta_vix", "beta_fxvol")):
                scale = 1e4 if k == 0 else 1.0
                res[f"{name}_{lab}"] = scale * r["beta"][k]
                res[f"{name}_{lab}_t"] = r["beta"][k] / r["se"][k]
            fit = X @ r["beta"]
            res[f"{name}_R2"] = 1 - np.var(s[name] - fit) / np.var(s[name])
        res["theta_UB_raw"] = 1 - s.H10.mean() / s.U.mean()
        res["theta_UB_adjusted"] = 1 - res["H10_alpha_bp"] / res["U_alpha_bp"]
        kind, lo, hi = ratio_set(s.H10.to_numpy(), s.U.to_numpy(), Fm)
        res.update(theta_adj_set=kind, theta_adj_lo=lo, theta_adj_hi=hi)
        res["corr_U_rvix"] = float(np.corrcoef(s.U, s.r_vix)[0, 1])
        res["corr_H10_rvix"] = float(np.corrcoef(s.H10, s.r_vix)[0, 1])
        out[sample] = res
    out["factor_ar1_rho_sigma_fx"] = {"primary": rho, "extended": rho}

    # Data check: Cboe second-month settlement against LSEG VXc2 at month-ends.
    vxc2 = read_raw(raw / "vix" / "VXc2.csv")
    vxc2 = vxc2[~vxc2.index.duplicated(keep="last")]
    col = "SETTLE" if "SETTLE" in vxc2 else ("TRDPRC_1" if "TRDPRC_1" in vxc2 else vxc2.columns[0])
    chk = []
    for t in months.date:
        live = [e for e in contracts if e > t]
        if len(live) >= 2 and t in vxc2.index:
            chk.append((price_on(contracts[live[1]], t), float(vxc2.loc[t, col])))
    chk = np.array([c for c in chk if np.all(np.isfinite(c))])
    check = {"n": len(chk), "lseg_field": col, "mean_abs_diff": float(np.mean(np.abs(chk[:, 0] - chk[:, 1]))),
             "share_equal_within_0.05": float(np.mean(np.abs(chk[:, 0] - chk[:, 1]) <= 0.05)),
             "corr": float(np.corrcoef(chk[:, 0], chk[:, 1])[0, 1])}

    fmt = lambda v: f"{v:.6g}" if isinstance(v, (float, np.floating)) else str(v)
    lines = ["# E5: systemic-risk exposure (restricted)", "", "```",
             pd.DataFrame({k: v for k, v in out.items() if k in ("primary", "extended")}).to_string(float_format=fmt), "```", "",
             f"AR(1) coefficient of the FX-volatility series: {rho:.4f}", "",
             f"Archive rescaling detected (file, first date at index scale): {rescaled}", "",
             "## Data check: Cboe second-month settlement vs LSEG VXc2 at month-ends", "", "```",
             pd.Series(check).to_string(float_format=fmt), "```", ""]
    (d / "e5_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
