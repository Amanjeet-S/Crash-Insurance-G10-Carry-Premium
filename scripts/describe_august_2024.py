"""Descriptive account of smiles, φ and carry returns around 5 August 2024.

The research design (section 6, with E5) names a descriptive account of
smiles, φ_t and returns around 5 August 2024. It was not computed with the
main results; the research log of 2 October 2026 records why it is computed
now. The rules below were fixed before the script was run.

- Window: New York business days from 15 July to 30 August 2024. Quotes are
  composite, selected on each day under the same five-business-day rule as
  at month-ends (``build_month_end_inputs`` with the days as sample dates).
- Smiles: on each day, the one-month at-the-money volatility and the absolute
  25Δ risk reversal, averaged over the nine currencies.
- φ_t: on each day, a market-reading SABR smile (β = 1) is calibrated for
  every currency with the calibration set, warm-started from the previous
  day, and checked for static arbitrage on the ±4 and ±10 ATM standard
  deviation grids as at month-ends, with non-convergence and arbitrage flags
  counted, not dropped (design, section 3; check added on 3 October 2026), and the three-long, three-short portfolio is formed on that day's
  forward discounts as at a month-end; C_skew_t, FD_t and φ_t are those of
  E1 for one-month 10Δ protection bought on that day.
- Returns: the position opened at the 31 July 2024 month-end, the one held
  over 5 August in the primary sample. Its unhedged value on each day is the
  leg-averaged forward return with the day's spot in place of the spot at
  expiry, and its options are marked at their intrinsic value on that day;
  neither is a traded price, and the time value of the options is ignored.
  The realised month (unhedged and 10Δ-hedged returns to expiry) is taken
  from the Stage 4 series.

The 31 July 2024 calibration must reproduce the stored E1 legs and skew
price, or the script stops. Every value reported is an average over at least
six currencies. Outputs are written to data/private/results/2026-09-23/:
august_2024_days.csv and august_2024_summary.md.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from estimate_stage4 import daily_spot_mid, spot_on  # noqa: E402

from qef.data.panel import ny_business_days  # noqa: E402
from qef.data.smile_inputs import build_month_end_inputs  # noqa: E402
from qef.fx.arbitrage import check_smile  # noqa: E402
from qef.fx.conventions import G10  # noqa: E402
from qef.fx.crash import forward_discount, forward_return, leg_skew_cost, option_payoff, rank_legs  # noqa: E402
from qef.fx.sabr import sabr_vol  # noqa: E402
from qef.fx.smile import SmileQuotes, calibrate_sabr  # noqa: E402

START, END = pd.Timestamp("2024-07-15"), pd.Timestamp("2024-08-30")
POSITION_DATE = pd.Timestamp("2024-07-31")
EVENT = pd.Timestamp("2024-08-05")
N_LEGS = 3


def calibrate_days(inputs: pd.DataFrame):
    fits, warm, flags = {}, {}, {"missing_input": 0, "not_converged": 0, "arbitrage_inner": 0, "arbitrage_wide": 0}
    for _, row in inputs.sort_values(["currency", "date"]).iterrows():
        if not row.complete_calib:
            flags["missing_input"] += 1
            continue
        q = SmileQuotes(row.F, row.tau, row.atm, row.rr25, row.bf25, 0.25, row.df_base, G10[row.currency].delta)
        res = calibrate_sabr(q, "market", x0=warm.get(row.currency))
        if not res.success:
            flags["not_converged"] += 1
            continue
        sm = res.smile
        sd = row.atm * np.sqrt(row.tau)
        for label, width in (("inner", 4.0), ("wide", 10.0)):
            rep = check_smile(sm.vol, row.F, row.tau, -width * sd, width * sd, n=2001)
            if not (rep.convex_ok and rep.slope_bounds_ok and rep.density_ok):
                flags[f"arbitrage_{label}"] += 1
        warm[row.currency] = np.array([np.log(sm.alpha), np.arctanh(sm.rho), np.log(sm.nu)])
        fits[(row.currency, row.date)] = sm
    return fits, flags


def day_portfolio(day_rows: pd.DataFrame, fits: dict, t):
    avail = {r.currency: r for _, r in day_rows.iterrows() if (r.currency, t) in fits}
    if len(avail) < 2 * N_LEGS:
        return None
    fd = {c: forward_discount(r.S, r.F, G10[c].usd_base) for c, r in avail.items()}
    longs, shorts = rank_legs(fd, N_LEGS)
    c_skew, legs = 0.0, []
    for c, is_long in [(c, True) for c in longs] + [(c, False) for c in shorts]:
        r, sm = avail[c], fits[(c, t)]
        vol = lambda K, r=r, sm=sm: sabr_vol(K, r.F, r.tau, sm.alpha, sm.rho, sm.nu, 1.0)
        K, vs, vf = leg_skew_cost(vol, r.F, r.tau, r.df_base, c, r.atm, is_long, 0.10)
        c_skew += (vs - vf) / N_LEGS
        legs.append((c, is_long, K, r.F, vs))
    FD = (sum(fd[c] for c in longs) - sum(fd[c] for c in shorts)) / N_LEGS
    return {"longs": " ".join(longs), "shorts": " ".join(shorts), "C_skew": c_skew, "FD": FD,
            "phi": c_skew / FD if FD > 0 else np.nan, "legs": legs}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    args = p.parse_args()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    raw = ROOT / "data" / "private" / "lseg" / args.retrieval_date / "raw"
    days = ny_business_days(START, END)
    inputs = build_month_end_inputs(raw, list(days))
    fits, flags = calibrate_days(inputs)
    spots = {c: daily_spot_mid(raw, c) for c in G10}

    e1 = pd.read_csv(d / "e1_series.csv", parse_dates=["date"])
    ref = e1[(e1["sample"] == "primary") & (e1.reading == "market") & (e1.delta == 0.10) & (e1.date == POSITION_DATE)].iloc[0]
    pos = day_portfolio(inputs[inputs.date == POSITION_DATE], fits, POSITION_DATE)
    if pos["longs"] != ref.longs or pos["shorts"] != ref.shorts or abs(pos["C_skew"] - ref.C_skew) > 1e-10:
        raise SystemExit("the 31 July 2024 calibration does not reproduce the stored E1 month")

    rows = []
    for t in days:
        day = inputs[inputs.date == t]
        calib = day[day.complete_calib]
        row = {"date": t, "n_currencies": len(calib), "atm_mean": calib.atm.mean(), "abs_rr25_mean": calib.rr25.abs().mean()}
        port = day_portfolio(day, fits, t)
        if port is not None:
            row.update(longs=port["longs"], shorts=port["shorts"], C_skew=port["C_skew"], FD=port["FD"], phi=port["phi"],
                       same_legs_as_position=(port["longs"], port["shorts"]) == (pos["longs"], pos["shorts"]))
        if t >= POSITION_DATE:
            U = H = 0.0
            for c, is_long, K, F0, prem in pos["legs"]:
                S, _ = spot_on(spots[c], t)
                rx = forward_return(c, is_long, S, F0)
                U += rx / N_LEGS
                H += (rx + option_payoff(c, is_long, K, F0, S) - prem) / N_LEGS
            row.update(position_U_marked=U, position_H_marked=H)
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(d / "august_2024_days.csv", index=False)

    s4 = pd.read_csv(d / "stage4_months.csv", parse_dates=["date"])
    month = s4[(s4["sample"] == "primary") & (s4.date == POSITION_DATE)].iloc[0]
    o = out.set_index("date")
    pre = o.loc[:POSITION_DATE]
    post = o.loc[POSITION_DATE:]
    trough = post.position_U_marked.idxmin()
    peak_atm = o.atm_mean.idxmax()
    fmt = lambda v: f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
    bp = lambda x: f"{1e4 * x:.1f}"
    vp = lambda x: f"{100 * x:.2f}"
    lines = ["# Smiles, phi and carry returns around 5 August 2024 (restricted)", "",
             f"Window {START:%Y-%m-%d} to {END:%Y-%m-%d}: {len(out)} New York business days; calibrated smiles {len(fits)}.",
             f"Flags (counted, not dropped): missing input {flags['missing_input']}, not converged {flags['not_converged']}, "
             f"failing the static-arbitrage checks within 4 ATM s.d. {flags['arbitrage_inner']} and within 10 ATM s.d. {flags['arbitrage_wide']}.",
             f"Position of {POSITION_DATE:%Y-%m-%d}: long {pos['longs']}, short {pos['shorts']} (reproduces E1).", "",
             "## Smiles (averages over the nine currencies, volatility points)", "",
             f"- ATM: {vp(pre.atm_mean.mean())} on average from 15 to 31 July; {vp(o.loc[POSITION_DATE, 'atm_mean'])} on 31 July; "
             f"{vp(o.loc[EVENT, 'atm_mean'])} on 5 August; peak {vp(o.atm_mean.max())} on {peak_atm:%d %B}; {vp(o.atm_mean.iloc[-1])} on 30 August.",
             f"- Absolute 25-delta risk reversal: {vp(o.loc[POSITION_DATE, 'abs_rr25_mean'])} on 31 July; {vp(o.loc[EVENT, 'abs_rr25_mean'])} on 5 August; "
             f"peak {vp(o.abs_rr25_mean.max())} on {o.abs_rr25_mean.idxmax():%d %B}; {vp(o.abs_rr25_mean.iloc[-1])} on 30 August.", "",
             "## phi_t of one-month 10-delta protection on the portfolio formed each day", "",
             f"- C_skew (bp): {bp(o.loc[POSITION_DATE, 'C_skew'])} on 31 July; {bp(o.loc[EVENT, 'C_skew'])} on 5 August; "
             f"maximum {bp(o.C_skew.max())} on {o.C_skew.idxmax():%d %B}; {bp(o.C_skew.iloc[-1])} on 30 August.",
             f"- FD (bp): {bp(o.loc[POSITION_DATE, 'FD'])} on 31 July; {bp(o.loc[EVENT, 'FD'])} on 5 August; {bp(o.FD.iloc[-1])} on 30 August.",
             f"- phi: {o.loc[POSITION_DATE, 'phi']:.3f} on 31 July; {o.loc[EVENT, 'phi']:.3f} on 5 August; maximum {o.phi.max():.3f} on "
             f"{o.phi.idxmax():%d %B}; {o.phi.iloc[-1]:.3f} on 30 August; range over 15 to 31 July {pre.phi.min():.3f} to {pre.phi.max():.3f}.",
             f"- Days on which the daily portfolio has the legs of the 31 July position: {int(o.same_legs_as_position.sum())} of {int(o.same_legs_as_position.notna().sum())}.", "",
             "## The 31 July 2024 position (bp of notional, marked at the day's spot; options at intrinsic value)", "",
             f"- Unhedged: {bp(o.loc[EVENT, 'position_U_marked'])} on 5 August; lowest {bp(post.position_U_marked.min())} on {trough:%d %B}; "
             f"{bp(o.position_U_marked.iloc[-1])} on 30 August.",
             f"- With 10-delta options: {bp(o.loc[EVENT, 'position_H_marked'])} on 5 August; {bp(o.loc[trough, 'position_H_marked'])} on {trough:%d %B}; "
             f"{bp(o.position_H_marked.iloc[-1])} on 30 August.",
             f"- Realised month (Stage 4, to expiry): unhedged {bp(month.U)}, hedged {bp(month.H10)}, payoff term (i) {bp(month.c1)}, "
             f"skew term (iii) {bp(month.c3)}.", "",
             "## Daily values", "", "```",
             out.drop(columns=["longs", "shorts"]).to_string(index=False, float_format=fmt), "```", ""]
    (d / "august_2024_summary.md").write_text("\n".join(lines))
    print("\n".join(lines[:22]))


if __name__ == "__main__":
    main()
