"""Long ATM sample: ATM-hedged carry from the start of the ATM and forward series.

Research design, section 4: "The long ATM sample is ATM-hedged carry
(Burnside et al., 2011) from the start of the ATM and forward series." The
paper lists it as not reported. Citations below are inherited from the design,
references.md and the modules named; none was consulted again for this script.

Rules already fixed by the design, reused from the project's code:

- New York month-ends and the five-business-day substitution rule
  (qef.data.panel, through qef.data.smile_inputs.build_month_end_inputs).
- Spot S, the outright forward F = S + points × pip factor and volatility time
  τ = (expiry − trade)/365 from build_month_end_inputs; option expiry and
  delivery from qef.data.smile_inputs.option_dates.
- Ranking on the forward discount fd = log(X/F_X) and legs long the n highest
  and short the n lowest (qef.fx.crash.forward_discount and rank_legs;
  Lustig, Roussanov and Verdelhan, 2011, not consulted).
- The hedge is the 'Hatm' leg of scripts/estimate_stage4.py: each leg buys the
  protective option (qef.fx.crash.protective_option) at the delta-neutral
  straddle strike in the pair's convention, K = F e^{+σ²τ/2} for pips delta
  and F e^{−σ²τ/2} for premium-adjusted delta (qef.fx.gk.atm_dns_strike;
  Reiswich and Wystup, 2012, not consulted), with σ the quoted ATM volatility,
  bought at the Garman–Kohlhagen premium carried to delivery, V = φ[F Φ(φd+)
  − K Φ(φd−)]/F in USD per USD of forward notional (qef.fx.crash; the
  construction of Burnside et al., 2011, Section 1.2, as recorded in
  references.md, not consulted).
- Returns at the end-of-day spot on the option expiry date
  (estimate_stage4.spot_on: substitution within five New York business days;
  a window ending after the last spot observation is unrealised and
  excluded). Per month, with w = 1/n,
  HML^U = w Σ_legs rx and HML^H = w Σ_legs (rx + payoff − V).
- Newey–West standard errors with the automatic bandwidth
  (qef.stats.hac.mean_and_se) and the test-inversion set for
  θ_UB = 1 − mean(HML^H)/mean(HML^U),
  CS = {θ : |mean(HML^H − (1 − θ) HML^U)| ≤ 1.96 se_HAC} on the grid [−10, 10]
  (estimate_stage4.theta_confidence_set; Fieller, 1954, not consulted); a set
  that reaches the edge of the grid is reported as unbounded.
  Amendment, 27 September 2026, after an independent review: whether the set
  is bounded is now decided from the limit, not from the grid. The automatic
  Newey–West bandwidth is scale invariant, so as |θ| grows the statistic
  tends to |t_U| = |mean(HML^U)|/se_HAC(HML^U). When |t_U| < 1.96 the set is
  unbounded on both sides: its ends are reported as −inf and +inf, and its
  pieces are listed only within the scanned range [−10, 10]. Every θ far
  enough beyond the scan on either side belongs to the set, and the set
  between the scan edge and infinity is not described. When |t_U| > 1.96 the
  set is bounded, and the script stops (the grid is too narrow) if the set
  on the grid reaches an edge or if any θ of a log-spaced check with
  10 ≤ |θ| ≤ 10⁶ is accepted. At |t_U| = 1.96 exactly the limit decides
  nothing, and the script stops. Under the original rule the set for
  January 2009 to December 2021, which is unbounded, was printed with a
  finite lower end that was an artefact of the grid. The acceptance rule and
  the grid are unchanged, and no tolerance or estimate was tuned to results.

Rules fixed before computing the results. The design leaves these open. I
fixed them from the availability records of the private audit (coverage and
two-sidedness), before I computed any return, premium or hedged return of
the months before the primary start. The 159 primary-window months reproduce
the Stage 4 HML^U and Hatm series (rule 7), which I had already computed and
reported. The Stage 4 extended-sample Hatm (Fenics quotes, two and two, 72
months before the primary start) was also known. An earlier version of this
paragraph said that I fixed the rules before I looked at any return, premium
or hedged return of this sample, which overstated what I had not seen; I
corrected it on 27 September 2026, after an independent review. Output files
of the same names had been written on 27 September 2026 by an earlier run
of a draft of this script; I overwrote them without opening them. Relative to that draft, rule 5 gained the windows "2009-01
onwards" and "before the primary start", and rule 8 was added.

1. Quotes. Composite spot, composite one-month forward points and the
   composite one-month ATM volatility only. Fenics quotes are not spliced in,
   so each month uses one contributor. A currency enters the ranking at
   month-end t when all three are observed or substituted under the month-end
   rule and the ATM mid is positive. Risk reversals, butterflies and rates are
   not required.
2. Legs. With N currencies available at t, n = 3 when N ≥ 6 and n = 2 when N
   is 4 or 5; a month with N < 4 is excluded and counted. Three long and
   three short is the primary construction (design, section 5) and is used
   whenever it is feasible without a currency on both sides. Two and two is
   the extended-sample rule, used only when three and three is infeasible.
   With N = 6 (or 4) every available currency is held, split at the median
   forward discount. The rule mirrors the requirement of the other samples
   that at least twice the number of legs be available.
3. First month. The sample starts at the first New York month-end at which
   N ≥ 4 and runs to 31 August 2026, the last month-end before the retrieval;
   later months with N < 4 stay in the count as excluded. The retrieval
   requested history from 1 January 1995 (reports/data_audit.md), so the
   start of this sample is bounded by the retrieval, not necessarily by the
   provider's history. Checked on 1 October 2026: a coverage retrieval of
   spot, one-month forward points and one-month ATM volatility from
   1 January 1970 to 31 January 1995 (data/private/lseg/2026-10-01/) returned
   no composite ATM quote before 6 January 1995 for any currency, while spot
   quotes go back to 1971 and forward points to 1982. The first month-end of the sample is
   therefore also the first month-end of the provider's ATM series, as the
   design requires; the summary states this from that retrieval's manifest.
4. Discount factors. None enters, so months before the one-month rate series
   start are kept, and the rate, D_q and D_b columns of
   build_month_end_inputs are not used. D_b enters the delta convention only
   as a common factor: the delta-neutral-straddle strike solves
   Δ_call + Δ_put = 0, which in spot delta is D_b[Φ(d+) − Φ(−d+)] = 0 for
   pips delta and D_b (K/F)[Φ(d−) − Φ(−d−)] = 0 for premium-adjusted delta,
   so d+ = 0 or d− = 0 and the strike does not depend on D_b (forward delta
   gives the same strike). The premium carried to delivery is the
   undiscounted premium divided by F (qef.fx.crash). For EURUSD-type pairs
   the premium D_q V_fwd/F is paid in USD and carried at the USD rate behind
   D_q. For USD-base pairs the premium in USD is D_q V_fwd/S = D_b V_fwd/F,
   carried at the USD rate implied by the forward, D_b = F D_q/S. Either way
   the rate cancels, as in the Hatm leg of estimate_stage4.py. The forward
   return and the payoff need only F and the spot at expiry. The number of
   used legs whose quote-currency rate is missing is reported.
5. Sub-periods, by formation month-end t: the start to December 2007; 2008,
   the formation month-ends of January to December 2008, whose returns are
   realised from February 2008 to January 2009; January 2009 onwards;
   January 2009 to December 2021; January 2022 onwards, the hiking regime of
   the design (the break of estimate_e1.py); the full sample; the full
   sample without 2008; before the primary start (31 May 2013); and the
   primary window, from the primary start. Means are also reported by leg
   count. Twelve months of 2008 give an imprecise Newey–West standard error;
   that window is reported to show what it contributes, not as a test.
6. Reference θ₀. The diffusive null of R6 for the ATM hedge: the mean over
   months of the average over the 2n legs of the hedge's forward delta
   Φ(φ d+) at the ATM strike and volatility, the estimate_stage4.py formula
   applied to the ATM hedge instead of the 10Δ hedge.
   Amendment, 27 September 2026, after an independent review: the rule is
   unchanged, but for USD-base legs it is an approximation. Their returns are
   in USD per USD, so the hedge ratio in USD returns is the premium-adjusted
   forward delta of the quoted pair, (K/F)Φ(φ d−), and the pips form Φ(φ d+)
   is its leading-order approximation. At the premium-adjusted
   delta-neutral-straddle strike d− = 0, so the exact ratio is e^{−σ²τ/2}/2
   against Φ(φσ√τ). As a diagnostic that changes no result, the months file
   also holds θ₀ with (K/F)Φ(φ d−) for USD-base legs (theta0_pa), and the
   summary reports its largest difference from θ₀ over the windows. No
   tolerance or estimate was tuned to results.
7. Checks. The rebuilt inputs must agree with smile_inputs.csv on S, F, τ and
   the ATM volatility (relative difference below 1e-14). In the overlap with
   the primary sample, every month whose legs coincide with those of the E1
   series must reproduce HML^U and Hatm of stage4_months.csv within 1e-12;
   otherwise the script stops. The correlations of HML^U and Hatm with the
   primary series over the overlap are reported, and those with the
   extended-sample Hatm (Fenics, two and two) as a comparison under
   different rules.
8. Counts and a diagnostic. As the design requires (section 4), substituted
   month-end quotes and excluded currency-months are counted. A
   currency-month is an exclusion when the currency was available at an
   earlier month-end of the sample and is not available at t; month-ends
   before a currency's first available one are not exclusions. As a
   diagnostic that changes no result, I count the leg-months whose selected
   ATM quote lies in a run of five or more unchanged business days, the
   audit's staleness rule for butterflies (robustness.stale_flags) applied to
   the ATM quote.

Outputs are LSEG-derived and are written to data/private/results/2026-09-23/:
long_atm_months.csv (month-level, private) and long_atm_summary.md (pooled
statistics). On 27 September 2026, after the same review, I also added to
the summary, without changing any estimate: a note for windows whose mean
HML^U is not positive, where θ_UB = −(mean(HML^H) − mean(HML^U))/mean(HML^U)
reverses sign, so that a positive value means the hedge raised the mean
return and is not a cost share; a line saying that the start is set by the
retrieval (rule 3); labels that name formation month-ends; and the rule 7
tolerance in basis points beside the differences. Reproduction, after
estimate_e1.py and estimate_stage4.py:

    .venv/bin/python scripts/estimate_long_atm.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
COVERAGE_RETRIEVAL = "2026-10-01"  # rule 3: pre-1995 coverage of the long-sample series
sys.path.insert(0, str(ROOT / "scripts"))
from estimate_e1 import PRIMARY_START, REGIME_BREAK  # noqa: E402
from estimate_stage4 import daily_spot_mid, spot_on, theta_confidence_set  # noqa: E402

from qef.data.panel import on_business_days, ny_month_ends, read_raw, sample_month_ends, stale_mask, two_sided  # noqa: E402
from qef.data.smile_inputs import build_month_end_inputs, option_dates  # noqa: E402
from qef.fx.conventions import G10  # noqa: E402
from qef.fx.crash import forward_discount, forward_return, option_payoff, protective_option, rank_legs  # noqa: E402
from qef.fx.gk import atm_dns_strike, d_plus_minus, forward_premium  # noqa: E402
from qef.stats.hac import mean_and_se  # noqa: E402

FIELDS = ("spot", "fwd", "atm")  # the quotes a currency needs (rule 1)
CHECK_COLUMNS = ("S", "F", "tau", "atm")
STALE_RUN = 5  # business days, design section 4
THETA_GRID = theta_confidence_set.__defaults__[0]  # the grid of estimate_stage4.theta_confidence_set
Z = 1.96  # the critical value of estimate_stage4.theta_confidence_set
BEYOND_SCAN = np.geomspace(10.0, 1e6, 601)  # |θ| checked beyond the grid for a bounded set (amendment)
RULE7_TOL = 1e-12  # rule 7, in decimal return units


def n_legs(n_ccy: int) -> int:
    """Legs per side for N available currencies (rule 2); 0 excludes the month."""
    return 3 if n_ccy >= 6 else 2 if n_ccy >= 4 else 0


def availability(inputs: pd.DataFrame) -> np.ndarray:
    """Currency-months with spot, forward points and a positive ATM mid (rule 1)."""
    ok = np.ones(len(inputs), dtype=bool)
    for f in FIELDS:
        ok &= (inputs[f"status_{f}"] != "missing").to_numpy()
    ok &= np.isfinite(inputs[["S", "F", "atm"]].to_numpy()).all(axis=1) & (inputs.atm > 0).to_numpy()
    return ok


def stale_atm(raw: Path, month_ends) -> set:
    """(currency, month-end) whose selected ATM quote is stale; robustness.stale_flags for the ATM (rule 8)."""
    out = set()
    for c in G10:
        f = on_business_days(read_raw(raw / "vol_atm" / f"{c}1MO=.csv"))
        mask = stale_mask(f[two_sided(f)], STALE_RUN)
        me = sample_month_ends(f, month_ends)
        me = me[me.status != "missing"]
        out |= {(c, m) for m, src in me.source_date.items() if bool(mask.get(src, False))}
    return out


def month_rows(inputs: pd.DataFrame, spots: dict, stale: set) -> pd.DataFrame:
    """One row per month-end from the first with at least four currencies (rule 3)."""
    inp = inputs.set_index(["date", "currency"]).sort_index()
    rows, started = [], False
    for t in inp.index.get_level_values("date").unique():
        avail = inp.loc[t]
        avail = avail[avail.available]
        n = n_legs(len(avail))
        if not started and n == 0:
            continue
        started = True
        base = {"date": t, "n_ccy": len(avail), "n_legs": n}
        if n == 0:
            rows.append({**base, "status": "too_few_currencies"})
            continue
        fd = {c: forward_discount(r.S, r.F, G10[c].usd_base) for c, r in avail.iterrows()}
        longs, shorts = rank_legs(fd, n)
        w = 1.0 / n
        acc = dict(U=0.0, Hatm=0.0, payoff=0.0, premium=0.0, theta0=0.0, theta0_pa=0.0,
                   spot_sub=0, quote_sub=0, no_rate=0, stale_atm=0)
        status = "ok"
        for c, is_long in [(c, True) for c in longs] + [(c, False) for c in shorts]:
            r = avail.loc[c]
            conv = G10[c].delta
            expiry = option_dates(t, c)[2]
            S_end, sub = spot_on(spots[c], expiry)
            if not np.isfinite(S_end):
                status = "return_unrealised" if expiry > spots[c].index[-1] else "missing_spot"
                break
            phi = protective_option(c, is_long)
            rx = forward_return(c, is_long, S_end, r.F)
            k_atm = atm_dns_strike(r.F, r.atm, r.tau, conv)
            v_atm = float(forward_premium(r.F, k_atm, r.atm, r.tau, phi)) / r.F
            pay = option_payoff(c, is_long, k_atm, r.F, S_end)
            # the two accumulations below are those of estimate_stage4.month_rows
            acc["U"] += w * rx
            acc["Hatm"] += w * (rx + pay - v_atm)
            acc["payoff"] += w * pay
            acc["premium"] += w * v_atm
            d1, d2 = d_plus_minus(r.F, k_atm, r.atm, r.tau)
            acc["theta0"] += float(norm.cdf(phi * float(d1))) / (2 * n)  # rule 6
            # rule 6 diagnostic: premium-adjusted forward delta for USD-base legs
            pa = (k_atm / r.F) * norm.cdf(phi * float(d2)) if G10[c].usd_base else norm.cdf(phi * float(d1))
            acc["theta0_pa"] += float(pa) / (2 * n)
            acc["spot_sub"] += sub
            acc["quote_sub"] += int(sum(r[f"status_{f}"] == "substituted" for f in FIELDS))
            acc["no_rate"] += int(r["status_rate"] == "missing")
            acc["stale_atm"] += int((c, t) in stale)
        row = {**base, "status": status, "longs": " ".join(longs), "shorts": " ".join(shorts),
               "FD": (sum(fd[c] for c in longs) - sum(fd[c] for c in shorts)) / n}
        if status == "ok":
            row.update(acc)
        rows.append(row)
    return pd.DataFrame(rows)


WINDOWS = {
    "full": lambda d: d.notna(),
    "start to 2007-12": lambda d: d < "2008-01-01",
    "2008": lambda d: (d >= "2008-01-01") & (d < "2009-01-01"),
    "2009-01 onwards": lambda d: d >= "2009-01-01",
    "2009-01 to 2021-12": lambda d: (d >= "2009-01-01") & (d < REGIME_BREAK),
    "2022-01 onwards": lambda d: d >= REGIME_BREAK,
    "full without 2008": lambda d: (d < "2008-01-01") | (d >= "2009-01-01"),
    "before the primary start": lambda d: d < PRIMARY_START,
    "primary window": lambda d: d >= PRIMARY_START,
}


def window_stats(df: pd.DataFrame) -> dict:
    U, H = df.U.to_numpy(), df.Hatm.to_numpy()
    out = {"months": len(df), "first": f"{df.date.min():%Y-%m}", "last": f"{df.date.max():%Y-%m}"}
    for name, x in (("U", U), ("Hatm", H), ("Hatm-U", H - U), ("payoff", df.payoff.to_numpy()),
                    ("premium", df.premium.to_numpy())):
        r = mean_and_se(x)
        out[f"{name} mean"], out[f"{name} se"] = 1e4 * r["mean"], 1e4 * r["se"]
        out[f"{name} t"], out[f"{name} lag"] = r["mean"] / r["se"], r["lag"]
    out["theta_UB"] = 1 - H.mean() / U.mean()
    kind, lo, hi, pieces = theta_set(H, U)
    text = " u ".join(f"[{a:.3f}, {b:.3f}]" for a, b in pieces)
    if kind == "unbounded":
        text = (text or "nothing") + " (within the scan only; far enough out, both tails belong to the set)"
    out.update({"CS": kind, "CS lo": lo, "CS hi": hi, "CS pieces": text or "empty",
                "theta0": df.theta0.mean(), "theta0_pa": df.theta0_pa.mean()})
    return out


def accepted(H, U, th) -> bool:
    """The acceptance rule of estimate_stage4.theta_confidence_set at one θ."""
    r = mean_and_se(H - (1 - th) * U)
    return abs(r["mean"]) <= Z * r["se"]


def cs_pieces(H, U, grid=THETA_GRID) -> list:
    """Contiguous pieces of the test-inversion set for θ_UB on the grid.

    The acceptance rule is that of estimate_stage4.theta_confidence_set, which
    reports only the outer bounds. Because the automatic Newey–West lag is
    chosen afresh for each θ, the set need not be an interval.
    """
    idx = np.flatnonzero([accepted(H, U, th) for th in grid])
    if idx.size == 0:
        return []
    cut = np.flatnonzero(np.diff(idx) > 1)
    starts, ends = np.r_[idx[0], idx[cut + 1]], np.r_[idx[cut], idx[-1]]
    return [(float(grid[a]), float(grid[b])) for a, b in zip(starts, ends)]


def theta_set(H, U, grid=THETA_GRID):
    """Kind, ends and pieces of the test-inversion set for θ_UB (amended rule).

    Boundedness is decided from the limit: as |θ| grows the statistic tends to
    |t_U| = |mean(U)|/se_HAC(U), because the automatic bandwidth is scale
    invariant. With |t_U| < 1.96 the set is unbounded on both sides, the ends
    are −inf and +inf, and the pieces are those within the grid only. With
    |t_U| > 1.96 the set on the grid must stay inside the grid, and no θ of
    BEYOND_SCAN on either side may be accepted; otherwise the grid is too
    narrow and the script stops, as it does at |t_U| = 1.96 exactly.
    """
    rU = mean_and_se(U)
    t_U = abs(rU["mean"]) / rU["se"]
    grid_kind, grid_lo, grid_hi = theta_confidence_set(H, U, grid)
    pieces = cs_pieces(H, U, grid)
    same = (not pieces and grid_kind == "empty") or (
        bool(pieces) and pieces[0][0] == grid_lo and pieces[-1][1] == grid_hi
        and (grid_kind == "unbounded" or (len(pieces) > 1) == (grid_kind == "union of intervals")))
    if not same:
        raise SystemExit("theta set: the pieces disagree with estimate_stage4.theta_confidence_set")
    if t_U == Z:
        raise SystemExit("theta set: |t_U| equals 1.96, where the limit does not decide boundedness")
    if t_U < Z:
        return "unbounded", -np.inf, np.inf, pieces
    if grid_kind == "unbounded":
        raise SystemExit(f"theta set: |t_U| = {t_U:.3f} > 1.96 but the set reaches the edge of the grid; "
                         "the grid is too narrow")
    if any(accepted(H, U, s * th) for th in BEYOND_SCAN for s in (-1.0, 1.0)):
        raise SystemExit(f"theta set: |t_U| = {t_U:.3f} > 1.96 but a theta beyond the grid is accepted; "
                         "the grid is too narrow")
    return grid_kind, grid_lo, grid_hi, pieces


def input_check(inputs: pd.DataFrame, d: Path) -> str:
    """Rule 7: the rebuilt inputs agree with those behind the calibration, E1 and Stage 4."""
    ref = pd.read_csv(d / "smile_inputs.csv", parse_dates=["date"]).set_index(["date", "currency"])
    mine = inputs.set_index(["date", "currency"]).loc[ref.index]
    for col in CHECK_COLUMNS:
        a, b = mine[col].to_numpy(float), ref[col].to_numpy(float)
        both = np.isfinite(a) & np.isfinite(b)
        if not (np.array_equal(np.isfinite(a), np.isfinite(b)) and np.all(np.abs(a[both] - b[both]) <= 1e-14 * np.abs(b[both]))):
            raise SystemExit(f"input check: {col} differs from smile_inputs.csv")
    return f"Rebuilt inputs agree with smile_inputs.csv on {', '.join(CHECK_COLUMNS)} over {len(ref)} currency-months."


def overlap_check(months: pd.DataFrame, d: Path) -> tuple[dict, str]:
    """Rule 7: reproduction of the primary Hatm where the legs coincide, and correlations."""
    load = lambda name: pd.read_csv(d / name, parse_dates=["date"])
    s4, e1 = load("stage4_months.csv"), load("e1_series.csv")
    e1 = e1[(e1["sample"] == "primary") & (e1.reading == "market") & (e1.delta == 0.10) & (e1.status == "ok")].set_index("date")
    mine = months[months.status == "ok"].set_index("date")
    out, note = {}, ""
    for sample in ("primary", "extended"):
        ref = s4[(s4["sample"] == sample) & (s4.status == "ok")].set_index("date")
        common = ref.index.intersection(mine.index)
        dH = (mine.loc[common, "Hatm"] - ref.loc[common, "Hatm"]).abs()
        dU = (mine.loc[common, "U"] - ref.loc[common, "U"]).abs()
        res = {"n_reference_months": len(ref), "n_overlap": len(common),
               "corr_Hatm": float(np.corrcoef(mine.loc[common, "Hatm"], ref.loc[common, "Hatm"])[0, 1]),
               "corr_U": float(np.corrcoef(mine.loc[common, "U"], ref.loc[common, "U"])[0, 1]),
               "max_abs_diff_Hatm_bp": 1e4 * float(dH.max()), "max_abs_diff_U_bp": 1e4 * float(dU.max())}
        if sample == "primary":
            same = np.array([t in e1.index and mine.loc[t, "longs"] == e1.loc[t, "longs"]
                             and mine.loc[t, "shorts"] == e1.loc[t, "shorts"] for t in common])
            res["n_same_legs"] = int(same.sum())
            if not (same.any() and dH[same].max() <= RULE7_TOL and dU[same].max() <= RULE7_TOL):
                raise SystemExit("overlap check: the long ATM sample does not reproduce the primary Hatm where legs coincide")
            if (~same).any():
                res["max_abs_diff_Hatm_bp_different_legs"] = 1e4 * float(dH[~same].max())
            res["rule 7 tolerance_bp"] = 1e4 * RULE7_TOL
            exp = lambda v: f"{v:.0e}".replace("e-0", "e-")
            note = (f"Rule 7 tolerance: {exp(RULE7_TOL)} in decimal return units, that is {exp(1e4 * RULE7_TOL)} bp, "
                    f"on HML^U and Hatm in the {int(same.sum())} months whose legs coincide; largest differences "
                    f"there {1e4 * float(dH[same].max()):.3g} bp (Hatm) and {1e4 * float(dU[same].max()):.3g} bp "
                    "(HML^U); passed.")
        out[sample] = res
    return out, note


def counts(months: pd.DataFrame, inputs: pd.DataFrame) -> dict:
    ok = months[months.status == "ok"]
    out = {"month-ends in sample": len(months),
           "first month-end": f"{months.date.min():%Y-%m}", "last month-end": f"{months.date.max():%Y-%m}",
           "first formation month-end with a realised return": f"{ok.date.min():%Y-%m}",
           "last formation month-end with a realised return": f"{ok.date.max():%Y-%m}"}
    out.update({f"status {k}": int(v) for k, v in months.status.value_counts().items()})
    out.update({f"realised months with n = {k}": int(v) for k, v in ok.n_legs.value_counts().sort_index().items()})
    for n in sorted(ok.n_legs.unique()):
        s = ok[ok.n_legs == n]
        out[f"n = {n}: first, last formation month-end with a realised return"] = (
            f"{s.date.min():%Y-%m} to {s.date.max():%Y-%m}")
    out.update({f"month-ends with N = {k} currencies": int(v)
                for k, v in months.n_ccy.value_counts().sort_index().items()})
    # rule 8: exclusions and substitutions over all currency-months of the sample
    inp = inputs[inputs.date >= months.date.min()]
    a = inp.pivot(index="date", columns="currency", values="available").astype(int)
    out["currency-months available"] = int(a.to_numpy().sum())
    out["currency-months excluded after the currency's first available month-end"] = int(
        ((a.cummax() == 1) & (a == 0)).to_numpy().sum())
    av = inp[inp.available]
    for f in FIELDS:
        out[f"substituted {f} quotes, available currency-months"] = int((av[f"status_{f}"] == "substituted").sum())
    out.update({"leg-months": int(2 * ok.n_legs.sum()),
                "substituted month-end quotes on legs (spot, forward, ATM)": int(ok.quote_sub.sum()),
                "spot substitutions at expiry": int(ok.spot_sub.sum()),
                "leg-months without a quote-currency rate (not needed, rule 4)": int(ok.no_rate.sum()),
                "leg-months with a stale ATM quote (diagnostic, rule 8)": int(ok.stale_atm.sum()),
                "  of which before the primary start": int(ok.loc[ok.date < PRIMARY_START, "stale_atm"].sum())})
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--start", default="1995-01-01", help="first date searched; the retrieval starts here")
    p.add_argument("--end", default="2026-08-31")
    args = p.parse_args()
    d = ROOT / "data" / "private" / "results" / args.retrieval_date
    raw = ROOT / "data" / "private" / "lseg" / args.retrieval_date / "raw"
    month_ends = ny_month_ends(args.start, args.end)
    inputs = build_month_end_inputs(raw, month_ends, quote_status=True)
    inputs["available"] = availability(inputs)
    check_inputs = input_check(inputs, d)
    spots = {c: daily_spot_mid(raw, c) for c in G10}
    months = month_rows(inputs, spots, stale_atm(raw, month_ends))
    months.to_csv(d / "long_atm_months.csv", index=False)
    overlap, overlap_note = overlap_check(months, d)

    ok = months[months.status == "ok"]
    table = {name: window_stats(ok[mask(ok.date)]) for name, mask in WINDOWS.items()}
    for n in sorted(ok.n_legs.unique()):
        table[f"n = {n} months"] = window_stats(ok[ok.n_legs == n])
    table = pd.DataFrame(table).T
    ret_cols = ["months", "first", "last", "U mean", "U se", "U t", "Hatm mean", "Hatm se", "Hatm t",
                "Hatm-U mean", "Hatm-U se"]
    theta_cols = ["theta_UB", "CS", "CS lo", "CS hi", "theta0", "premium mean", "premium se",
                  "payoff mean", "payoff se", "U lag", "Hatm lag"]
    pieces_col = table["CS pieces"]

    def fmt(v):
        if isinstance(v, (float, np.floating)):
            return "+inf" if v == np.inf else f"{v:.4g}"
        return str(v)

    # rule 3: the start is set by the retrieval (amendment of 27 September 2026)
    first, searched = months.date.min(), pd.Timestamp(args.start)
    requested = f"history requested from {searched.day} {searched:%B %Y}"
    start_note = (f"The first month-end, {first:%Y-%m}, is the first month-end of the retrieval ({requested}), "
                  "not necessarily the start of the provider's series." if first == month_ends[0] else
                  f"The first month-end, {first:%Y-%m}, is the first with at least four currencies ({requested}).")
    # rule 3, checked on 1 October 2026 against a coverage retrieval from 1 January 1970
    coverage = ROOT / "data" / "private" / "lseg" / COVERAGE_RETRIEVAL / "manifest.json"
    if coverage.exists():
        import json
        man = json.loads(coverage.read_text())
        files = [(Path(f["file"]).parts, f["first_date"]) for f in man["files"] if f.get("first_date")]
        firsts = {blk: min((d for parts, d in files if len(parts) > 2 and parts[1] == blk), default=None)
                  for blk in ("vol_atm", "spot", "forward")}
        if firsts["vol_atm"] is not None and pd.Timestamp(firsts["vol_atm"]) > month_ends[0] - pd.offsets.MonthBegin(1):
            start_note = (f"The first month-end, {first:%Y-%m}, is the first month-end of the retrieval ({requested}) "
                          f"and of the provider's composite one-month ATM series: a coverage retrieval from "
                          f"{man['requested_start']} (retrieval of {COVERAGE_RETRIEVAL}) returned no ATM quote before "
                          f"{firsts['vol_atm']}, while spot and forwards begin on {firsts['spot']} and "
                          f"{firsts['forward']} at the earliest.")
    # theta_UB where mean HML^U <= 0 (amendment of 27 September 2026); the value itself is as pre-registered
    neg = table[table["U mean"].astype(float) <= 0]
    theta_notes = ["theta_UB = 1 - mean(Hatm)/mean(HML^U) = -(mean(Hatm) - mean(HML^U))/mean(HML^U), as pre-registered."]
    if len(neg):
        theta_notes += [
            "Where mean HML^U <= 0 the sign is reversed: a positive theta_UB means that the hedge raised the mean",
            "return, so it is not a cost share and is not comparable with the rows where mean HML^U > 0.",
            "Windows with mean HML^U <= 0: " + "; ".join(
                f"{k} (mean HML^U {r['U mean']:.1f} bp, Hatm - U {r['Hatm-U mean']:+.1f} bp, theta_UB {r['theta_UB']:.4g})"
                for k, r in neg.iterrows()) + "."]
    # rule 6 diagnostic (amendment of 27 September 2026)
    diff = (table["theta0_pa"] - table["theta0"]).astype(float)
    k = diff.abs().idxmax()
    theta_notes += [
        "Rule 6 diagnostic, which changes no result: theta0 with the premium-adjusted forward delta (K/F)Phi(phi d-)",
        "for USD-base legs, their hedge ratio in USD returns, differs from theta0 by at most "
        f"{abs(diff[k]):.4f} over the windows",
        f"({diff[k]:+.4f} in '{k}'; the differences lie between {diff.min():+.4f} and {diff.max():+.4f})."]
    lines = ["# Long ATM sample: ATM-hedged carry (restricted)", "",
             "Rules: scripts/estimate_long_atm.py (docstring). Composite spot, forward points and ATM volatility;",
             "n = 3 legs per side when at least six currencies are available, n = 2 with four or five.",
             "Returns in basis points per month, USD per USD of forward notional at delivery; Newey-West",
             "standard errors with the automatic bandwidth; windows, and their first and last months, by",
             "formation month-end.",
             "CS: 95% test-inversion set for theta_UB, scanned on the grid [-10, 10]; bounded or unbounded by the",
             "limit |t_U| against 1.96 (see the note on the pieces).", "",
             "## Sample", "", "```", pd.Series(counts(months, inputs)).to_string(), "```", "", start_note, "",
             "## HML^U and ATM-hedged HML^H by window (bp per month)", "", "```",
             table[ret_cols].to_string(float_format=fmt), "```", "",
             "## theta_UB, hedge premium and payoff by window", "", "```",
             table[theta_cols].to_string(float_format=fmt), "```", "", *theta_notes, "",
             "Pieces of each test-inversion set on the scan [-10, 10], grid step 0.001. A set is unbounded on",
             "both sides when |t_U| = |mean HML^U|/se_HAC(HML^U) (U t above) is below 1.96: the automatic",
             "bandwidth is scale invariant, so the statistic tends to |t_U| as |theta| grows. Such a set has",
             "ends -inf and +inf, and its pieces are listed only within the scan, so an end at -10.000 or",
             "10.000 is the edge of the scan, not a bound. Every theta far enough beyond the scan on either side",
             "belongs to the set; the set between the scan edge and infinity is not described. A bounded set",
             "(|t_U| > 1.96) reaches neither edge of the scan, and no theta of a log-spaced check with",
             "10 <= |theta| <= 1e6 is accepted.", "", "```",
             *(f"{k:<26}{v}" for k, v in pieces_col.items()), "```", "",
             "## Checks", "", check_inputs, "",
             "Overlap with the Stage 4 series (primary: rules coincide; extended: Fenics ATM, two and two):", "", "```",
             pd.DataFrame(overlap).to_string(float_format=lambda v: f"{v:.6g}"), "```", "", overlap_note, ""]
    (d / "long_atm_summary.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
