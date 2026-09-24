"""Stage 4: carry portfolio returns, E2 (predictability) and E3 (hedge-cost decomposition).

Research design, sections 5 to 7. Legs are those of the E1 series (market
reading). For each month-end t and leg:

- the forward return and the option payoffs are evaluated at the spot quote
  on the option expiry date (value on delivery), from the daily spot mid;
- the protective option is struck at the smile 10Δ strike (primary), the 25Δ
  strike or the ATM (delta-neutral straddle) strike, and bought at the smile
  premium carried to delivery;
- realised volatility σ̂P is the annualised standard deviation (√252) of daily
  log spot changes over the 21 New York business days ending at t.

Month-level quantities, in USD per USD of forward notional at delivery:
HML^U = (1/n) Σ signed forward returns; HML^H = HML^U + (1/n) Σ (payoff − V_smile);
E3 terms (1/n) Σ (payoff − V^GK(σ̂P)), −(1/n) Σ (V_flat − V^GK(σ̂P)) and
−(1/n) Σ (V_smile − V_flat), which add up to HML^H − HML^U.

E2 regresses next-period HML^U on φ_t (primary predictor) and on the oriented
10Δ risk reversal, in sample with Newey–West errors and a residual
bootstrap under the null of no predictability (Stambaugh, 1999, eq. 18, gives
the first-order bias used as a check), and out of sample with an expanding
window and the MSPE-adjusted statistic of Clark and West (2007). E3 reports θ_UB = 1 − mean(HML^H)/mean(HML^U)
with a test-inversion confidence set and the three hedge-cost terms.
Outputs are written to data/private/results/.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from qef.data.panel import NewYorkCalendar, read_raw
from qef.data.smile_inputs import option_dates
from qef.fx.conventions import G10
from qef.fx.crash import forward_return, leg_skew_cost, option_payoff, protective_option
from qef.fx.gk import atm_dns_strike, d_plus_minus, forward_premium
from qef.fx.sabr import sabr_vol
from qef.stats.hac import mean_and_se, ols_hac

ROOT = Path(__file__).resolve().parents[1]
FIRST_OOS = pd.Timestamp("2016-12-30")  # forecast made at end-December 2016 for January 2017
BDAY = pd.offsets.CustomBusinessDay(calendar=NewYorkCalendar())


def daily_spot_mid(raw: Path, ccy: str) -> pd.Series:
    f = read_raw(raw / "spot" / f"{ccy}=.csv")
    f = f[~f.index.duplicated(keep="last")].sort_index()
    ok = f["BID"].notna() & f["ASK"].notna()
    return ((f.loc[ok, "BID"] + f.loc[ok, "ASK"]) / 2).astype(float)


def spot_on(series: pd.Series, date: pd.Timestamp, window: int = 5):
    """Spot on `date`, or the last quote within `window` New York business days before it.

    A date after the last available observation has not been realised yet and
    returns NaN, so incomplete return windows at the end of the sample are excluded.
    """
    if date > series.index[-1]:
        return np.nan, -1
    if date in series.index:
        return float(series.loc[date]), 0
    prior = series.loc[:date]
    if prior.empty or len(pd.date_range(prior.index[-1], date, freq=BDAY)) - 1 > window:
        return np.nan, -1
    return float(prior.iloc[-1]), 1


def realised_vol(series: pd.Series, date: pd.Timestamp, n: int = 21) -> float:
    s = series.loc[:date]
    s = s[s.index.dayofweek < 5]
    r = np.diff(np.log(s.to_numpy()[-(n + 1):]))
    return float(r.std(ddof=1) * np.sqrt(252.0)) if r.size == n else np.nan


def month_rows(e1, inputs, panel, spots, n_legs):
    inp = inputs.set_index(["currency", "date"])
    fit = panel[panel.reading == "market"].set_index(["currency", "date"])
    rows = []
    for _, m in e1.iterrows():
        t = m.date
        legs = [(c, True) for c in m.longs.split()] + [(c, False) for c in m.shorts.split()]
        acc = dict(U=0.0, H10=0.0, H25=0.0, Hatm=0.0, c1=0.0, c2=0.0, c3=0.0, theta0=0.0, rr_or=0.0, spot_sub=0, vol_ok=True)
        for c, is_long in legs:
            r, f = inp.loc[(c, t)], fit.loc[(c, t)]
            conv = G10[c].delta
            vol = lambda K, r=r, f=f: sabr_vol(K, r.F, r.tau, f.alpha, f.rho, f.nu, 1.0)
            expiry = option_dates(t, c)[2]
            S_end, sub = spot_on(spots[c], expiry)
            if not np.isfinite(S_end):
                acc = None
                break
            acc["spot_sub"] += sub
            phi = protective_option(c, is_long)
            rx = forward_return(c, is_long, S_end, r.F)
            K10, vs10, vf10 = leg_skew_cost(vol, r.F, r.tau, r.df_base, c, r.atm, is_long, 0.10)
            K25, vs25, _ = leg_skew_cost(vol, r.F, r.tau, r.df_base, c, r.atm, is_long, 0.25)
            k_atm = atm_dns_strike(r.F, r.atm, r.tau, conv)
            v_atm = float(forward_premium(r.F, k_atm, r.atm, r.tau, phi)) / r.F
            sig_p = realised_vol(spots[c], t)
            if not np.isfinite(sig_p):
                acc["vol_ok"] = False
                sig_p = r.atm
            vp10 = float(forward_premium(r.F, K10, sig_p, r.tau, phi)) / r.F
            pay10 = option_payoff(c, is_long, K10, r.F, S_end)
            w = 1.0 / n_legs
            acc["U"] += w * rx
            acc["H10"] += w * (rx + pay10 - vs10)
            acc["H25"] += w * (rx + option_payoff(c, is_long, K25, r.F, S_end) - vs25)
            acc["Hatm"] += w * (rx + option_payoff(c, is_long, k_atm, r.F, S_end) - v_atm)
            acc["c1"] += w * (pay10 - vp10)
            acc["c2"] += -w * (vf10 - vp10)
            acc["c3"] += -w * (vs10 - vf10)
            d1, _ = d_plus_minus(r.F, K10, float(vol(K10)), r.tau)
            # R6 diffusive null per leg: |forward delta| of the hedge; averaged over the 2n legs
            acc["theta0"] += float(norm.cdf(phi * float(d1))) / (2 * n_legs)
            # oriented 10Δ risk reversal: protective-option vol minus opposite-option vol (quoted)
            acc["rr_or"] += w * (r.rr10 if phi == 1 else -r.rr10)
        if acc is None:
            rows.append({"date": t, "status": "missing_spot"})
            continue
        rows.append({"date": t, "status": "ok", "phi": m.phi, "C_skew": m.C_skew, "FD": m.FD, **acc})
    return pd.DataFrame(rows)


def stambaugh_bootstrap(y, x, B, seed=20260924):
    """Residual bootstrap of the predictive slope under H0: no predictability.

    Stambaugh (1999) shows that the slope estimate is biased when the predictor
    is persistent and its innovations are correlated with returns. The null
    distribution is simulated with y*_k = ȳ + û_k and x*_{k+1} = ĉ + ρ̂ x*_k + v̂_{k+1},
    where û are residuals from the unrestricted regression of y on x (so they
    carry no predictability) and v̂ are AR(1) residuals of x. The pairs
    (û_k, v̂_{k+1}) are drawn independently with replacement, which keeps their
    contemporaneous correlation and removes any link between returns and past
    predictor shocks. The analytic first-order bias −(σ_uv/σ_v²)(1 + 3ρ)/T is
    reported alongside as a check.
    """
    y, x = np.asarray(y, float), np.asarray(x, float)
    T = len(y)
    X = np.column_stack([np.ones(T), x])
    a_hat, b_hat = np.linalg.lstsq(X, y, rcond=None)[0]
    u_full = y - a_hat - b_hat * x
    c, rho = np.linalg.lstsq(np.column_stack([np.ones(T - 1), x[:-1]]), x[1:], rcond=None)[0]
    v = x[1:] - c - rho * x[:-1]
    u = u_full[:-1]
    s_uv, s_vv = np.cov(u, v)[0, 1], v.var(ddof=1)
    rng = np.random.default_rng(seed)
    n = T - 1
    draws = np.empty(B)
    for i in range(B):
        idx = rng.integers(n, size=n)
        us, vs = u[idx], v[idx]
        xs = np.empty(n)
        xs[0] = x[0]
        for k in range(n - 1):
            xs[k + 1] = c + rho * xs[k] + vs[k]
        ys = y.mean() + us
        draws[i] = np.linalg.lstsq(np.column_stack([np.ones(n), xs]), ys, rcond=None)[0][1]
    return {"b_hat": b_hat, "bias_bootstrap": draws.mean(), "bias_analytic": -(s_uv / s_vv) * (1 + 3 * rho) / T,
            "b_bias_corrected": b_hat - draws.mean(), "p_one_sided_b_gt_0": float((draws >= b_hat).mean()),
            "rho_predictor": rho, "corr_u_v": float(np.corrcoef(u, v)[0, 1])}


def clark_west(dates, y, x):
    """Expanding-window forecasts from FIRST_OOS; Clark and West (2007) MSPE-adjusted statistic."""
    ok = np.isfinite(x) & np.isfinite(y)
    dates, y, x = dates[ok], y[ok], x[ok]
    f, e_m, e_b = [], [], []
    for i in np.where(dates >= FIRST_OOS.to_datetime64())[0]:
        if i < 24:
            continue
        X = np.column_stack([np.ones(i), x[:i]])
        a, b = np.linalg.lstsq(X, y[:i], rcond=None)[0]
        yhat, ybar = a + b * x[i], y[:i].mean()
        e_m.append(y[i] - yhat)
        e_b.append(y[i] - ybar)
        f.append((y[i] - ybar) ** 2 - ((y[i] - yhat) ** 2 - (ybar - yhat) ** 2))
    f, e_m, e_b = map(np.asarray, (f, e_m, e_b))
    # Clark and West (2006, Section 2): regress f on a constant; for one-step forecasts the
    # usual least-squares standard error is used (primary); Newey-West reported alongside.
    t_ls = f.mean() / (f.std(ddof=1) / np.sqrt(len(f)))
    r = mean_and_se(f)
    return {"n_forecasts": len(f), "oos_r2": 1 - (e_m @ e_m) / (e_b @ e_b), "cw_mean": f.mean(),
            "cw_t": t_ls, "cw_p_one_sided": float(1 - norm.cdf(t_ls)),
            "cw_t_nw": r["mean"] / r["se"], "cw_p_one_sided_nw": float(1 - norm.cdf(r["mean"] / r["se"]))}


def theta_confidence_set(H, U, grid=np.linspace(-10, 10, 20001)):
    accept = []
    for th in grid:
        r = mean_and_se(H - (1 - th) * U)
        accept.append(abs(r["mean"]) <= 1.96 * r["se"])
    accept = np.array(accept)
    if not accept.any():
        return "empty", np.nan, np.nan
    inside = grid[accept]
    bounded = not (accept[0] or accept[-1])
    contiguous = np.all(np.diff(np.where(accept)[0]) == 1)
    kind = ("interval" if contiguous else "union of intervals") if bounded else "unbounded"
    return kind, float(inside.min()), float(inside.max())


def e3_summary(df):
    out = {}
    U, H = df.U.to_numpy(), df.H10.to_numpy()
    for name in ("U", "H10", "H25", "Hatm", "c1", "c2", "c3"):
        r = mean_and_se(df[name].to_numpy())
        out[f"mean_{name}_bp"], out[f"se_{name}_bp"] = 1e4 * r["mean"], 1e4 * r["se"]
    for hedge in ("H10", "H25", "Hatm"):
        out[f"theta_UB_{hedge}"] = 1 - df[hedge].mean() / df.U.mean()
    kind, lo, hi = theta_confidence_set(H, U)
    out.update(theta_UB_H10_set=kind, theta_UB_H10_lo=lo, theta_UB_H10_hi=hi)
    out["theta0_reference"] = df.theta0.mean()
    out["n_months"] = len(df)
    out["spot_substitutions"] = int(df.spot_sub.sum())
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--bootstrap", type=int, default=9999)
    args = p.parse_args()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    raw = ROOT / "data" / "private" / "lseg" / args.retrieval_date / "raw"
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    e1 = load("e1_series.csv")
    spots = {c: daily_spot_mid(raw, c) for c in G10}
    sel = lambda sample: e1[(e1["sample"] == sample) & (e1.reading == "market") & (e1.delta == 0.10) & (e1.status == "ok")]
    prim = month_rows(sel("primary"), load("smile_inputs.csv"), load("smile_panel.csv"), spots, 3)
    ext = month_rows(sel("extended"), load("smile_inputs_fn.csv"), load("smile_panel_fn.csv"), spots, 2)
    prim["sample"], ext["sample"] = "primary", "extended"
    months = pd.concat([ext, prim], ignore_index=True).sort_values("date")
    months.to_csv(d / "stage4_months.csv", index=False)
    pm = prim[prim.status == "ok"].sort_values("date")
    allm = months[months.status == "ok"].sort_values("date")

    e2 = {}
    for name in ("phi", "rr_or"):
        X = np.column_stack([np.ones(len(pm)), pm[name]])
        reg = ols_hac(pm.U.to_numpy(), X)
        e2[name] = {"b": reg["beta"][1], "se_hac": reg["se"][1], "t_hac": reg["beta"][1] / reg["se"][1], "hac_lag": reg["lag"],
                    **stambaugh_bootstrap(pm.U.to_numpy(), pm[name].to_numpy(), args.bootstrap),
                    **clark_west(allm.date.to_numpy(), allm.U.to_numpy(), allm[name].to_numpy())}
    e2["corr_phi_rr_or"] = {"value": float(np.corrcoef(pm.phi, pm.rr_or)[0, 1])}
    e3 = {"primary": e3_summary(pm), "extended": e3_summary(ext[ext.status == "ok"])}

    fmt = lambda v: f"{v:.6g}" if isinstance(v, (float, np.floating)) else str(v)
    lines = ["# Stage 4: returns, E2 and E3 (restricted)", "",
             f"Months: primary {len(pm)} ({pm.date.min():%Y-%m} to {pm.date.max():%Y-%m}); extended {int((ext.status == 'ok').sum())}", "",
             "## E2: predictive regressions of HML^U on the predictor (primary sample)", "", "```",
             pd.DataFrame(e2).T.to_string(float_format=fmt), "```", "",
             "## E3: hedge-cost decomposition (bp per month) and theta_UB", "", "```",
             pd.DataFrame(e3).to_string(float_format=fmt), "```", ""]
    (d / "stage4_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
