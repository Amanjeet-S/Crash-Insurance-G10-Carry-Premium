"""Stage 1 audit of the LSEG FX panel (data plan, audit checks 1 to 5, 8 and 9).

    .venv/bin/python scripts/audit_fx_panel.py [--retrieval-date 2026-09-23]

Reads data/private/lseg/<date>/ (raw tables, metadata.csv, requests.jsonl,
manifest.json) and writes data/private/audit/<date>/:

    coverage.csv        first and last date with a two-sided quote, per RIC
    two_sidedness.csv   month-end availability and substitutions under the five-day rule
    staleness.csv       runs of unchanged quotes; butterfly runs of five or more flagged
    conventions.csv     underlying pair, units, forward-point scaling, risk-reversal sign
    splice.csv          Fenics and TIFO minus composite on the overlap
    plausibility.csv    crossed quotes, negative volatilities, forward-point signs
    samples.csv         month-end availability behind the primary and extended samples
    audit_report.md     restricted summary with counts and dates

Month-end rules follow the research design, section 4. Daily statistics use
quotes dated on New York business days. Every output is LSEG-derived and
stays in data/private/.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from qef.data.panel import (
    CURRENCIES,
    TENORS,
    VOL_QUOTES,
    instrument_catalogue,
    mid,
    ny_business_days,
    ny_month_ends,
    on_business_days,
    parse_ric,
    read_raw,
    sample_month_ends,
    stale_mask,
    two_sided,
    unchanged_runs,
)
from qef.fx.conventions import G10

ROOT = Path(__file__).resolve().parents[1]
WINDOW = 5  # business days for month-end substitution (research design, section 4)
STALE_RUN = 5  # business days of unchanged quotes that flag a butterfly
SIGN_MIN_DIFF = 0.10  # percentage points; smaller deposit differentials are not sign-checked
SCALING_MIN_DIFF = 1.0  # percentage points; differentials used to test point scaling
FIRST_MONTH = "1995-01-01"
CODES = {v: k for k, v in VOL_QUOTES.items()}  # quote type -> RIC code
SMILE = ("atm", "rr25", "bf25", "rr10", "bf10")


def vol_ric(ccy, quote, tenor="1M", contributor=""):
    return f"{ccy}{tenor}{CODES[quote]}={contributor}"


def d(x):
    """Date as ISO string, or empty when missing."""
    return "" if x is None or pd.isna(x) else str(pd.Timestamp(x).date())


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Loading


def load(src: Path):
    manifest = json.loads((src / "manifest.json").read_text())
    mismatched = [f["file"] for f in manifest["files"]
                  if not (src / f["file"]).exists() or sha256(src / f["file"]) != f["sha256"]]
    frames = {p.stem: read_raw(p) for p in sorted((src / "raw").glob("*/*.csv"))}
    meta = pd.read_csv(src / "metadata.csv").set_index("requested_ric")
    requests = pd.read_json(src / "requests.jsonl", lines=True)
    return manifest, mismatched, frames, meta, requests


class Panel:
    """Raw tables with cached business-day views, mids and month-end samples."""

    def __init__(self, frames, month_ends):
        self.raw = frames
        self.month_ends = month_ends
        self._bd, self._me = {}, {}

    def has(self, ric):
        return ric in self.raw and two_sided(self.raw[ric]).any()

    def bd(self, ric):
        if ric not in self._bd:
            self._bd[ric] = on_business_days(self.raw[ric])
        return self._bd[ric]

    def mid(self, ric):
        f = self.bd(ric)
        return mid(f).dropna()

    def me(self, ric):
        """Month-end sample over the full month-end range (missing if no data)."""
        if ric not in self._me:
            frame = self.raw.get(ric, pd.DataFrame(columns=["BID", "ASK"]))
            self._me[ric] = sample_month_ends(frame, self.month_ends, WINDOW)
        return self._me[ric]

    def available(self, ric):
        return self.me(ric)["status"] != "missing"


# ---------------------------------------------------------------------------
# Checks 1 and 2: coverage and two-sidedness


def coverage(panel, catalogue, requests):
    errors = (requests[requests["kind"] == "history"].dropna(subset=["error"])
              .groupby("ric")["error"].first())
    rows = []
    for inst in catalogue:
        row = {"ric": inst.ric, "block": inst.block, "currency": inst.currency, "tenor": inst.tenor,
               "quote": inst.quote, "contributor": inst.contributor}
        f = panel.raw.get(inst.ric)
        if f is None:
            row.update(rows=0, note=f"no data: {errors.get(inst.ric, '')}"[:200])
            rows.append(row)
            continue
        ts = two_sided(f)
        bd = panel.bd(inst.ric)
        row.update(
            rows=len(f), n_duplicate_dates=int(f.index.duplicated().sum()),
            n_non_business_rows=len(f) - len(bd) - int(f.index.duplicated().sum()),
            n_weekend_rows=int((f.index.dayofweek >= 5).sum()),
            first_row=d(f.index.min()), last_row=d(f.index.max()),
            first_two_sided=d(f.index[ts].min()) if ts.any() else "",
            last_two_sided=d(f.index[ts].max()) if ts.any() else "",
            n_two_sided=int(ts.sum()),
        )
        if ts.any():
            days = ny_business_days(f.index[ts].min(), f.index[ts].max())
            row["share_business_days_two_sided"] = round(float(two_sided(bd).sum()) / len(days), 4)
        if inst.block == "vix" and "SETTLE" in f:
            s = f["SETTLE"].dropna()
            row.update(first_settle=d(s.index.min()), last_settle=d(s.index.max()), n_settle=len(s))
        rows.append(row)
    return pd.DataFrame(rows)


def two_sidedness(panel, catalogue, primary_start):
    rows = []
    for inst in catalogue:
        if inst.block == "vix" or not panel.has(inst.ric):
            continue
        s = panel.me(inst.ric)
        first = panel.raw[inst.ric].index[two_sided(panel.raw[inst.ric])].min()
        w = s[s.index >= first]
        row = {"ric": inst.ric, "block": inst.block, "currency": inst.currency, "tenor": inst.tenor,
               "quote": inst.quote, "contributor": inst.contributor,
               "window_start": d(w.index.min()), "window_end": d(w.index.max()),
               "n_month_ends": len(w), "n_two_sided_on_day": int((w.status == "observed").sum()),
               "share_two_sided_on_day": round(float((w.status == "observed").mean()), 4) if len(w) else np.nan,
               "n_substituted": int((w.status == "substituted").sum()),
               "n_missing": int((w.status == "missing").sum()),
               "max_lag": w["lag"].max()}
        if primary_start is not None:
            p, e = s[s.index >= primary_start], w[w.index < primary_start]
            row.update(n_month_ends_primary=len(p), n_substituted_primary=int((p.status == "substituted").sum()),
                       n_missing_primary=int((p.status == "missing").sum()),
                       n_month_ends_pre_primary=len(e),
                       n_substituted_pre_primary=int((e.status == "substituted").sum()),
                       n_missing_pre_primary=int((e.status == "missing").sum()))
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Samples (research design, section 4)


def primary_sample(panel):
    required = [vol_ric(c, q) for c in CURRENCIES for q in SMILE] + [f"{c}1M=" for c in CURRENCIES]
    avail = pd.DataFrame({r: panel.available(r) for r in required})
    strict = pd.DataFrame({r: panel.me(r)["status"] == "observed" for r in required})
    all_ok = avail.all(axis=1)
    start = all_ok.idxmax() if all_ok.any() else None
    strict_start = strict.all(axis=1).idxmax() if strict.all(axis=1).any() else None
    first_avail = {r: avail.index[avail[r]].min() if avail[r].any() else pd.NaT for r in required}
    out = {"required": required, "avail": avail, "start": start, "strict_start": strict_start,
           "first_avail": first_avail}
    if start is not None:
        prev = avail.index[avail.index < start]
        out["binding"] = [r for r in required if len(prev) and not avail.loc[prev[-1], r]]
        after = avail.loc[start:]
        out["months"] = len(after)
        out["incomplete_months"] = [d(m) for m in after.index[~after.all(axis=1)]]
        # A currency is excluded in a month if any of its required series is missing.
        out["exclusions"] = {c: int((~after[[r for r in required if r.startswith(c)]].all(axis=1)).sum())
                             for c in CURRENCIES}
    return out


def extended_sample(panel, primary_start):
    """Fenics availability per currency (full smile set and the calibration set)."""
    sets = {"full": SMILE, "calibration": ("atm", "rr25", "bf25")}
    per_ccy, counts = {}, {}
    for name, quotes in sets.items():
        avail = pd.DataFrame({c: pd.concat([panel.available(vol_ric(c, q, contributor="FN")) for q in quotes]
                                           + [panel.available(f"{c}1M=")], axis=1).all(axis=1)
                              for c in CURRENCIES})
        per_ccy[name] = avail
        counts[name] = avail.sum(axis=1)
    out = {}
    for name, avail in per_ccy.items():
        n = counts[name]
        before = n[n.index < primary_start] if primary_start is not None else n
        ok = before[before >= 4]
        start = ok.index.min() if len(ok) else None
        gaps = [d(m) for m in before.index[(before.index >= start) & (before < 4)]] if start is not None else []
        out[name] = {
            "first_by_ccy": {c: avail.index[avail[c]].min() if avail[c].any() else pd.NaT for c in CURRENCIES},
            "start_4ccy": start,
            "ccys_at_start": [c for c in CURRENCIES if start is not None and avail.loc[start, c]],
            "months_before_primary": int((before >= 4).sum()),
            "gaps": gaps,
            "count": n,
        }
    return out


def long_atm(panel):
    rows = {}
    for c in CURRENCIES:
        fwd = panel.available(f"{c}1M=")
        comp = (panel.available(vol_ric(c, "atm")) & fwd)
        fen = (panel.available(vol_ric(c, "atm", contributor="FN")) & fwd)
        either = comp | fen
        rows[c] = {k: (s.index[s].min() if s.any() else pd.NaT)
                   for k, s in (("composite", comp), ("fenics", fen), ("either", either))}
    return rows


# ---------------------------------------------------------------------------
# Check 3: staleness


def staleness(panel, catalogue, primary_start):
    rows = []
    for inst in catalogue:
        if inst.block == "vix" or not panel.has(inst.ric):
            continue
        bd = panel.bd(inst.ric)
        q = bd[two_sided(bd)]
        runs = unchanged_runs(q)
        mask = stale_mask(q, STALE_RUN)
        long = runs.loc[runs["length"].idxmax()]
        me = panel.me(inst.ric)
        me = me[me.status != "missing"]
        me_stale = mask.reindex(pd.DatetimeIndex(me["source_date"]), fill_value=False).to_numpy(dtype=bool)
        row = {"ric": inst.ric, "block": inst.block, "currency": inst.currency, "tenor": inst.tenor,
               "quote": inst.quote, "contributor": inst.contributor,
               "flag_rule_applies": inst.quote in ("bf25", "bf10"),
               "n_obs": len(q),
               "share_unchanged_from_previous": round((len(q) - len(runs)) / max(len(q) - 1, 1), 4),
               "n_runs_ge5": int((runs["length"] >= STALE_RUN).sum()),
               "share_stale_days": round(float(mask.mean()), 4),
               "longest_run": int(long["length"]), "longest_run_start": d(long["start"]),
               "longest_run_end": d(long["end"]),
               "n_month_ends": len(me), "n_stale_month_ends": int(me_stale.sum())}
        if primary_start is not None:
            inp = me.index >= primary_start
            qp, qe = mask[mask.index >= primary_start], mask[mask.index < primary_start]
            row.update(n_month_ends_primary=int(inp.sum()), n_stale_month_ends_primary=int(me_stale[inp].sum()),
                       share_stale_days_primary=round(float(qp.mean()), 4) if len(qp) else np.nan,
                       n_month_ends_pre_primary=int((~inp).sum()),
                       n_stale_month_ends_pre_primary=int(me_stale[~inp].sum()),
                       share_stale_days_pre_primary=round(float(qe.mean()), 4) if len(qe) else np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Checks 4 and 5: underlying pair, units, scaling and risk-reversal sign


def meta_pair(meta, ric):
    if ric not in meta.index or not bool(meta.loc[ric, "found"]):
        return ""
    return f"{meta.loc[ric, 'FirstCurrency']}{meta.loc[ric, 'SecondCurrency']}"


def deposit_differential(panel, ccy, usd_base):
    """Quote-currency minus base-currency one-month deposit rate, per cent (mid)."""
    r_c, r_u = panel.mid(f"{ccy}1MD="), panel.mid("USD1MD=")
    return (r_c - r_u) if usd_base else (r_u - r_c)


def forward_frame(panel, ccy):
    spot, pts = panel.mid(f"{ccy}="), panel.mid(f"{ccy}1M=")
    f = pd.concat({"spot": spot, "points": pts}, axis=1).dropna()
    days = panel.bd(f"{ccy}1M=").get("DAYS_MAT")
    f["days"] = (days.reindex(f.index) if days is not None else pd.Series(np.nan, index=f.index)).fillna(30.0)
    return f


def realised_forward(logs, n=21):
    """Annualised volatility (per cent) of daily log changes over the next n business days."""
    return logs.diff().rolling(n).std().shift(-n) * np.sqrt(252) * 100


def conventions(panel, meta):
    rows = []
    for c in CURRENCIES:
        cfg = G10[c]
        vol_rics = [vol_ric(c, q, t, s) for q in SMILE for t in TENORS for s in ("", "FN")]
        if c == "EUR":
            vol_rics += ["EUR1MO=TIFO", "EUR1MRR=TIFO", "EUR1MBF=TIFO"]
        found = [r for r in vol_rics if r in meta.index and bool(meta.loc[r, "found"])]
        vol_pairs = sorted({meta_pair(meta, r) for r in found})
        und = sorted({str(meta.loc[r, "UnderlyingQuoteRIC"]) for r in found
                      if pd.notna(meta.loc[r, "UnderlyingQuoteRIC"])})
        scaling = meta.loc[f"{c}1M=", "ScalingFactor"] if f"{c}1M=" in meta.index else np.nan
        row = {"currency": c, "pair_config": cfg.pair, "usd_base_config": cfg.usd_base,
               "spot_pair_metadata": meta_pair(meta, f"{c}="),
               "forward_pair_metadata": meta_pair(meta, f"{c}1M="),
               "vol_pairs_metadata": "|".join(vol_pairs), "n_vol_rics_found": len(found),
               "n_vol_rics_requested": len(vol_rics),
               "vol_underlying_rics_metadata": "|".join(und),
               "vol_pair_matches_config": vol_pairs == [cfg.pair],
               "forward_scaling_metadata": scaling, "pip_factor_config": cfg.pip_factor,
               "pip_factor_matches": bool(np.isclose(1.0 / scaling, cfg.pip_factor)) if pd.notna(scaling) else False,
               "delta_convention_config": ("spot, premium-adjusted" if cfg.delta.premium_adjusted else "spot, pips"),
               "delta_convention_source": "Reiswich and Wystup (2012), reporting Clark (2011); provisional"}

        # Scaling test: implied (F/S - 1)·360/days against the deposit differential.
        if panel.has(f"{c}1M=") and panel.has(f"{c}=") and panel.has(f"{c}1MD=") and panel.has("USD1MD="):
            f = forward_frame(panel, c)
            implied = f["points"] / scaling / f["spot"] * 360 / f["days"] * 100
            diff = deposit_differential(panel, c, cfg.usd_base).reindex(f.index)
            use = diff.abs() >= SCALING_MIN_DIFF
            ratio = (implied[use] / diff[use]).dropna()
            row.update(implied_to_deposit_diff_median=round(float(ratio.median()), 3) if len(ratio) else np.nan,
                       implied_to_deposit_diff_iqr=(f"{ratio.quantile(.25):.3f} to {ratio.quantile(.75):.3f}"
                                                    if len(ratio) else ""),
                       n_scaling_days=len(ratio))

        # Units: medians of the composite one-month quotes.
        if panel.has(vol_ric(c, "atm")):
            row["atm_1m_mid_median"] = round(float(panel.mid(vol_ric(c, "atm")).median()), 3)
        if panel.has(f"{c}1MD="):
            row["deposit_mid_median"] = round(float(panel.mid(f"{c}1MD=").median()), 3)

        # Risk-reversal sign: under "call minus put on the quoted pair" the RR carries the
        # sign of the quoted pair's spot-volatility correlation.
        atm, rr = panel.mid(vol_ric(c, "atm")), panel.mid(vol_ric(c, "rr25"))
        logs = np.log(panel.mid(f"{c}="))
        j = pd.concat({"ls": logs, "atm": atm, "rr": rr}, axis=1).dropna()
        if len(j) > 250:
            ch = j[["ls", "atm"]].diff().dropna()
            corr = ch["ls"].corr(ch["atm"])
            yearly = [(np.sign(g["rr"].mean()), np.sign(ch.loc[ch.index.isin(g.index), "ls"]
                                                        .corr(ch.loc[ch.index.isin(g.index), "atm"])))
                      for _, g in j.groupby(j.index.year) if len(g) > 100]
            row.update(rr25_1m_mean=round(float(j["rr"].mean()), 3),
                       rr25_1m_share_positive=round(float((j["rr"] > 0).mean()), 3),
                       corr_dlogspot_datm=round(float(corr), 3),
                       rr_sign_matches_spot_vol_corr=bool(np.sign(j["rr"].mean()) == np.sign(corr)),
                       share_years_sign_match=round(float(np.mean([a == b for a, b in yearly])), 3) if yearly else np.nan,
                       n_years=len(yearly))
        rr10, bf25, bf10 = (panel.mid(vol_ric(c, q)) for q in ("rr10", "bf25", "bf10"))
        k = pd.concat({"rr25": rr, "rr10": rr10}, axis=1).dropna()
        k = k[(k["rr25"].abs() > 0.1)]
        if len(k):
            row.update(rr10_rr25_sign_agreement=round(float((np.sign(k["rr10"]) == np.sign(k["rr25"])).mean()), 3),
                       rr10_over_rr25_median=round(float((k["rr10"] / k["rr25"]).median()), 3))
        b = pd.concat({"bf25": bf25, "bf10": bf10}, axis=1).dropna()
        b = b[b["bf25"] > 0.05]
        if len(b):
            row["bf10_over_bf25_median"] = round(float((b["bf10"] / b["bf25"]).median()), 3)

        # Empirical pair check: 1M ATM against realised volatility of the USD pair and of the EUR cross.
        if c != "EUR" and panel.has(vol_ric(c, "atm")):
            eurusd = np.log(panel.mid("EUR="))
            ls = logs
            usd_base_sign = 1.0 if cfg.usd_base else -1.0  # log EUR/ccy = log EURUSD + log USD/ccy
            cross = eurusd + usd_base_sign * ls
            r = pd.concat({"atm": atm, "rv_usd": realised_forward(ls), "rv_eur": realised_forward(cross)}, axis=1).dropna()
            if len(r) > 250:
                row.update(rmse_atm_vs_rv_usd_pair=round(float(np.sqrt(((r.atm - r.rv_usd) ** 2).mean())), 3),
                           rmse_atm_vs_rv_eur_cross=round(float(np.sqrt(((r.atm - r.rv_eur) ** 2).mean())), 3),
                           mean_rv_usd_pair=round(float(r.rv_usd.mean()), 3),
                           mean_rv_eur_cross=round(float(r.rv_eur.mean()), 3),
                           mean_atm_1m=round(float(r.atm.mean()), 3))
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Check 8: splice


def splice(panel, month_ends):
    pairs = [(vol_ric(c, q, t, "FN"), vol_ric(c, q, t), "fenics-composite")
             for c in CURRENCIES for t in TENORS for q in SMILE]
    for q in ("atm", "rr25", "bf25"):
        pairs.append((vol_ric("EUR", q, "1M", "TIFO"), vol_ric("EUR", q), "tifo-composite"))
        pairs.append((vol_ric("EUR", q, "1M", "TIFO"), vol_ric("EUR", q, "1M", "FN"), "tifo-fenics"))
    rows = []
    for a, b, label in pairs:
        inst = parse_ric(b)
        row = {"comparison": label, "ric": a, "reference_ric": b, "currency": inst.currency,
               "tenor": inst.tenor, "quote": inst.quote}
        if not (panel.has(a) and panel.has(b)):
            rows.append({**row, "n": 0})
            continue
        fa, fb = panel.bd(a), panel.bd(b)
        j = pd.concat({"m_a": panel.mid(a), "m_b": panel.mid(b),
                       "s_a": (fa["ASK"] - fa["BID"]), "s_b": (fb["ASK"] - fb["BID"])}, axis=1).dropna()
        diff = j["m_a"] - j["m_b"]
        me = diff[diff.index.isin(month_ends)]
        row.update(n=len(diff), first_overlap=d(diff.index.min()) if len(diff) else "",
                   last_overlap=d(diff.index.max()) if len(diff) else "")
        if len(diff):
            row.update(mean=round(float(diff.mean()), 4), sd=round(float(diff.std()), 4),
                       mean_abs=round(float(diff.abs().mean()), 4), median=round(float(diff.median()), 4),
                       share_equal_mid=round(float((diff.abs() < 1e-12).mean()), 4),
                       mean_spread=round(float(j["s_a"].mean()), 4),
                       mean_spread_reference=round(float(j["s_b"].mean()), 4),
                       n_month_ends=len(me), mean_month_ends=round(float(me.mean()), 4) if len(me) else np.nan,
                       sd_month_ends=round(float(me.std()), 4) if len(me) > 1 else np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Check 9: plausibility


def plausibility(panel, catalogue, month_ends):
    rows = []
    for inst in catalogue:
        if inst.ric not in panel.raw:
            continue
        bd = panel.bd(inst.ric)
        row = {"ric": inst.ric, "block": inst.block, "currency": inst.currency, "tenor": inst.tenor,
               "quote": inst.quote, "contributor": inst.contributor}
        if inst.block == "vix":
            s = bd["SETTLE"].dropna() if "SETTLE" in bd else pd.Series(dtype=float)
            row.update(n_obs=len(s), n_non_positive=int((s <= 0).sum()),
                       min_value=s.min() if len(s) else np.nan, max_value=s.max() if len(s) else np.nan)
            rows.append(row)
            continue
        q = bd[two_sided(bd)]
        bid, ask = q["BID"].astype(float), q["ASK"].astype(float)
        crossed = q.index[bid > ask]
        m = (bid + ask) / 2
        row.update(n_obs=len(q), n_crossed=len(crossed), first_crossed=d(crossed.min()) if len(crossed) else "",
                   last_crossed=d(crossed.max()) if len(crossed) else "",
                   n_zero_width=int((bid == ask).sum()),
                   min_mid=m.min() if len(m) else np.nan, max_mid=m.max() if len(m) else np.nan,
                   median_spread=(ask - bid).median() if len(q) else np.nan)
        if inst.quote == "atm":
            neg = q.index[(bid < 0) | (ask < 0) | (m <= 0)]
            row.update(n_negative_vol=len(neg), first_negative_vol=d(neg.min()) if len(neg) else "",
                       n_atm_above_100=int((m > 100).sum()))
        if inst.quote in ("bf25", "bf10"):
            row["n_negative_bf_mid"] = int((m < 0).sum())
        if inst.block == "spot":
            row["n_non_positive"] = int((m <= 0).sum())
        if inst.block == "forward" and panel.has(f"{inst.currency}1MD=") and panel.has("USD1MD="):
            diff = deposit_differential(panel, inst.currency, G10[inst.currency].usd_base).reindex(m.index)
            chk = pd.concat({"pts": m, "diff": diff}, axis=1).dropna()
            chk = chk[(chk["diff"].abs() >= SIGN_MIN_DIFF) & (chk["pts"] != 0)]
            bad = chk[np.sign(chk["pts"]) != np.sign(chk["diff"])]
            row.update(n_sign_checked=len(chk), n_sign_inconsistent=len(bad),
                       share_sign_inconsistent=round(len(bad) / len(chk), 4) if len(chk) else np.nan,
                       first_sign_inconsistent=d(bad.index.min()) if len(bad) else "",
                       last_sign_inconsistent=d(bad.index.max()) if len(bad) else "",
                       max_abs_diff_inconsistent=round(float(bad["diff"].abs().max()), 3) if len(bad) else np.nan,
                       date_max_abs_diff_inconsistent=d(bad["diff"].abs().idxmax()) if len(bad) else "",
                       n_sign_checked_month_ends=int(chk.index.isin(month_ends).sum()),
                       n_sign_inconsistent_month_ends=int(bad.index.isin(month_ends).sum()))
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Report


def md_table(df: pd.DataFrame, columns=None) -> str:
    if columns is not None:
        df = df.reindex(columns=columns)
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    def fmt(v):
        if isinstance(v, float) and v.is_integer() and abs(v) >= 1:
            return str(int(v))
        return "" if pd.isna(v) else str(v)
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(fmt(v) for v in r) + " |")
    return "\n".join(lines)


def report(ctx) -> str:
    cov, ts, st, conv, spl, pl = (ctx[k] for k in ("coverage", "two_sided", "stale", "conv", "splice", "plaus"))
    prim, ext, latm, req = ctx["primary"], ctx["extended"], ctx["long_atm"], ctx["requests"]
    hist = req[req["kind"] == "history"]
    failed = cov[cov["rows"] == 0]
    chunked = (hist[hist["rows"] > 0].groupby("ric").size() > 1).sum()
    L = [f"# Private data audit: LSEG FX panel, retrieval {ctx['date']}", "",
         "Restricted. Every value, date and count below is LSEG-derived and stays in `data/private/` "
         "until the licence terms for publication are confirmed. Generated by `scripts/audit_fx_panel.py`.", "",
         "## Inputs", "",
         f"- RICs requested: {len(cov)}; with data: {len(cov) - len(failed)}; without data: {len(failed)}.",
         f"- History requests: {len(hist)}, of which {int(len(hist) - hist['ric'].nunique())} backward extensions "
         f"or retries; RICs whose history needed more than one non-empty request (truncation): {int(chunked)}.",
         f"- Manifest: {len(ctx['manifest']['files'])} files; hash mismatches: {len(ctx['mismatched'])}.",
         f"- Rows on dates that are not New York business days: {int(cov['n_non_business_rows'].fillna(0).sum())}, "
         f"of which on weekends {int(cov['n_weekend_rows'].fillna(0).sum())} (RICs with weekend rows: "
         f"{', '.join(cov.loc[cov['n_weekend_rows'].fillna(0) > 0, 'ric'])}); the rest are New York holidays. "
         f"Duplicate dates: {int(cov['n_duplicate_dates'].fillna(0).sum())}.",
         f"- Month-end range audited: {d(ctx['month_ends'].min())} to {d(ctx['month_ends'].max())} "
         "(last complete month before retrieval).", ""]
    if len(failed):
        L += ["RICs without data:", ""] + [f"- `{r.ric}`: {r.note}" for r in failed.itertuples()] + [""]

    L += ["## Underlying pairs, units and scaling (checks 4 and 5)", ""]
    L += [md_table(conv, ["currency", "pair_config", "vol_pairs_metadata", "n_vol_rics_found",
                          "forward_pair_metadata", "forward_scaling_metadata", "pip_factor_matches",
                          "implied_to_deposit_diff_median", "n_scaling_days"]), ""]
    ns = conv[conv.currency.isin(["NOK", "SEK"])]
    for r in ns.itertuples():
        L.append(f"- {r.currency}: metadata pair {r.vol_pairs_metadata} across {r.n_vol_rics_found} of "
                 f"{r.n_vol_rics_requested} volatility RICs (underlying RIC {r.vol_underlying_rics_metadata}). "
                 f"RMSE of 1M ATM against next-21-day realised volatility: {r.rmse_atm_vs_rv_usd_pair} "
                 f"(USD pair) and {r.rmse_atm_vs_rv_eur_cross} (EUR cross).")
    L += ["", "Risk-reversal sign (1M composite 25Δ; the RR should carry the sign of the quoted pair's "
          "spot-volatility correlation under the call-minus-put convention):", ""]
    L += [md_table(conv, ["currency", "rr25_1m_mean", "corr_dlogspot_datm", "rr_sign_matches_spot_vol_corr",
                          "share_years_sign_match", "rr10_rr25_sign_agreement", "rr10_over_rr25_median",
                          "bf10_over_bf25_median", "atm_1m_mid_median"]), ""]

    L += ["## Samples (research design, section 4)", ""]
    if prim["start"] is not None:
        fa = prim["first_avail"]
        latest = sorted(fa, key=lambda r: fa[r])[-5:]
        L += [f"- Primary start: **{d(prim['start'])}**, the first month-end at which all nine currencies have "
              "two-sided composite 1M ATM, 25Δ and 10Δ RR and BF, and forward points, with the five-day "
              f"substitution rule. Without substitution: {d(prim['strict_start'])}.",
              f"- Binding series (unavailable at the preceding month-end): {', '.join(prim['binding']) or 'none'}.",
              f"- Latest first availability among the required series: "
              + "; ".join(f"`{r}` {d(fa[r])}" for r in latest) + ".",
              f"- Primary window: {prim['months']} month-ends to {d(ctx['month_ends'].max())}. Month-ends with at "
              f"least one currency excluded: {len(prim['incomplete_months'])}"
              + (f" ({', '.join(prim['incomplete_months'][:12])}{' ...' if len(prim['incomplete_months']) > 12 else ''})"
                 if prim["incomplete_months"] else "") + ".",
              "- Currency-months excluded in the primary window: "
              + ", ".join(f"{k} {v}" for k, v in prim["exclusions"].items()) + ".", ""]
    else:
        L += ["- Primary start: no month-end satisfies the rule.", ""]
    rows = []
    for ccy in CURRENCIES:
        rows.append({"currency": ccy,
                     "primary set (composite)": d(max(prim["first_avail"][r] for r in prim["required"] if r.startswith(ccy))),
                     "Fenics full set": d(ext["full"]["first_by_ccy"][ccy]),
                     "Fenics ATM, 25Δ RR, 25Δ BF": d(ext["calibration"]["first_by_ccy"][ccy]),
                     "long ATM (composite or Fenics)": d(latm[ccy]["either"])})
    L += ["First month-end available, by currency (each set includes composite forward points):", "",
          md_table(pd.DataFrame(rows)), ""]
    for name, label in (("full", "full Fenics smile set (ATM, 25Δ and 10Δ RR and BF)"),
                        ("calibration", "Fenics calibration set (ATM, 25Δ RR and BF)")):
        e = ext[name]
        L.append(f"- Extended sample, {label}: first month-end with at least four currencies "
                 f"{d(e['start_4ccy']) or 'none'} ({', '.join(e['ccys_at_start'])}); month-ends before the primary "
                 f"start with four or more currencies: {e['months_before_primary']}; gaps: "
                 f"{', '.join(e['gaps'][:10]) or 'none'}.")
    ets = ts[(ts.contributor == "fenics") & (ts.tenor == "1M")]
    if "n_missing_pre_primary" in ets:
        miss = ets[ets.n_missing_pre_primary > 0].sort_values("n_missing_pre_primary", ascending=False)
        L.append("- Fenics 1M series with missing month-ends before the primary start (after the five-day rule): "
                 + "; ".join(f"`{r.ric}` {r.n_missing_pre_primary}" for r in miss.itertuples()) + ".")
    L.append("")

    L += ["## Two-sidedness and substitutions (check 2)", ""]
    req_ts = ts[ts.ric.isin(prim["required"])]
    if "n_substituted_primary" in req_ts:
        L.append(f"- Required primary series, primary window: {int(req_ts.n_substituted_primary.sum())} substitutions "
                 f"and {int(req_ts.n_missing_primary.sum())} missing month-ends across {len(req_ts)} series.")
        worst = req_ts.sort_values("n_substituted_primary", ascending=False).head(5)
        L.append("- Most substitutions: " + "; ".join(f"`{r.ric}` {r.n_substituted_primary}"
                                                     for r in worst.itertuples()) + ".")
    fn = ts[(ts.contributor == "fenics") & (ts.tenor == "1M")]
    if "n_substituted_pre_primary" in fn:
        L.append(f"- Fenics 1M series before the primary start: {int(fn.n_substituted_pre_primary.sum())} "
                 f"substitutions and {int(fn.n_missing_pre_primary.sum())} missing month-ends across {len(fn)} series "
                 f"({int(fn.n_month_ends_pre_primary.sum())} series-month-ends).")
    L.append("")

    L += ["## Staleness (check 3)", ""]
    bf = st[st.flag_rule_applies & (st.tenor == "1M") & (st.contributor == "composite")]
    if len(bf):
        L += [md_table(bf, ["ric", "share_stale_days", "share_stale_days_primary", "n_runs_ge5", "longest_run",
                            "longest_run_start", "n_stale_month_ends_primary", "n_month_ends_primary"]), ""]
    bff = st[st.flag_rule_applies & (st.tenor == "1M") & (st.contributor == "fenics")]
    if len(bff):
        L += ["Fenics 1M butterflies (the extended sample uses the window before the primary start):", "",
              md_table(bff, ["ric", "share_stale_days_pre_primary", "n_stale_month_ends_pre_primary",
                             "n_month_ends_pre_primary", "share_stale_days_primary", "longest_run",
                             "longest_run_start", "longest_run_end"]), ""]
    other = st[~st.flag_rule_applies & (st.tenor == "1M") & st.quote.isin(["atm", "rr25", "rr10"])]
    if "share_stale_days_pre_primary" in other:
        g = other.groupby("contributor")[["share_stale_days_pre_primary", "share_stale_days_primary"]].median()
        g = g.loc[[k for k in ("composite", "fenics") if k in g.index]]
        L.append("- Median share of stale days for 1M ATM and risk reversals (runs of five or more): "
                 + "; ".join(f"{k} {r.share_stale_days_pre_primary:.3f} before and {r.share_stale_days_primary:.3f} "
                             f"from the primary start" for k, r in g.iterrows()) + ".")
    L.append("")

    L += ["## Splice (check 8)", "", "Fenics minus composite, 1M, daily mids on the overlap (volatility points):", ""]
    s1 = spl[(spl.comparison == "fenics-composite") & (spl.tenor == "1M") & (spl.n > 0)]
    if len(s1):
        cell = s1.apply(lambda r: f"{r['mean']:+.3f} ({r['sd']:.3f}, {int(r['n'])})", axis=1)
        piv = s1.assign(cell=cell).pivot(index="currency", columns="quote", values="cell")
        piv = piv.reindex(index=list(CURRENCIES), columns=list(SMILE))
        L += ["Cells: mean (sd, n).", "", md_table(piv.reset_index()), ""]
    ti = spl[spl.comparison.str.startswith("tifo") & (spl.n > 0)]
    for r in ti.itertuples():
        L.append(f"- `{r.ric}` minus `{r.reference_ric}`: mean {r.mean:+.3f}, sd {r.sd:.3f}, n {r.n}, "
                 f"{r.first_overlap} to {r.last_overlap}.")
    L.append("")

    L += ["## Plausibility (check 9)", ""]
    L.append(f"- Crossed quotes (bid above ask): {int(pl.n_crossed.fillna(0).sum())} observations in "
             f"{int((pl.n_crossed.fillna(0) > 0).sum())} RICs.")
    L.append(f"- Zero-width quotes (bid equal to ask): {int(pl.n_zero_width.fillna(0).sum())} observations in "
             f"{int((pl.n_zero_width.fillna(0) > 0).sum())} RICs.")
    if "n_negative_vol" in pl:
        L.append(f"- Negative or zero ATM volatilities: {int(pl.n_negative_vol.fillna(0).sum())}.")
    if "n_negative_bf_mid" in pl:
        nb = pl[pl.n_negative_bf_mid.fillna(0) > 0]
        L.append(f"- Negative butterfly mids: {int(nb.n_negative_bf_mid.sum())} in {len(nb)} RICs ("
                 + ", ".join(f"`{r.ric}` {int(r.n_negative_bf_mid)}" for r in nb.itertuples()) + ").")
    fw = pl[pl.block == "forward"]
    if "n_sign_checked" in fw:
        L += ["", f"Forward points against the sign of the 1M deposit differential (days with |differential| "
              f"≥ {SIGN_MIN_DIFF} percentage points):", "",
              md_table(fw, ["ric", "n_sign_checked", "n_sign_inconsistent", "share_sign_inconsistent",
                            "max_abs_diff_inconsistent", "date_max_abs_diff_inconsistent",
                            "n_sign_inconsistent_month_ends"]), ""]
    worst = pl.sort_values("n_crossed", ascending=False).head(5)
    worst = worst[worst.n_crossed > 0]
    if len(worst):
        L.append("- Most crossed quotes: " + "; ".join(f"`{r.ric}` {int(r.n_crossed)} ({r.first_crossed} to "
                                                      f"{r.last_crossed})" for r in worst.itertuples()) + ".")
    L += ["", "Flagged observations are kept; the clean layer carries the flags.", "",
          "## Not covered here", "",
          "- Butterfly convention (check 7): provider documentation and the 10Δ diagnostic under both readings.",
          "- Delta and premium conventions (check 6): provider documentation; the configuration is provisional.",
          "- Quote revisions: need a second retrieval to compare with this one.", ""]
    return "\n".join(L)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--retrieval-date", default="2026-09-23")
    p.add_argument("--lseg-root", default=str(ROOT / "data" / "private" / "lseg"))
    p.add_argument("--out-root", default=str(ROOT / "data" / "private" / "audit"))
    args = p.parse_args()

    src = Path(args.lseg_root) / args.retrieval_date
    out = Path(args.out_root) / args.retrieval_date
    out.mkdir(parents=True, exist_ok=True)
    manifest, mismatched, frames, meta, requests = load(src)
    catalogue = instrument_catalogue()
    # Month-ends up to the last complete month before the retrieval date.
    month_ends = ny_month_ends(FIRST_MONTH, pd.Timestamp(args.retrieval_date) - pd.Timedelta(days=1))
    month_ends = month_ends[month_ends < pd.Timestamp(args.retrieval_date).to_period("M").start_time]
    panel = Panel(frames, month_ends)

    prim = primary_sample(panel)
    ctx = {"date": args.retrieval_date, "manifest": manifest, "mismatched": mismatched, "requests": requests,
           "month_ends": month_ends, "primary": prim}
    ctx["coverage"] = coverage(panel, catalogue, requests)
    ctx["two_sided"] = two_sidedness(panel, catalogue, prim["start"])
    ctx["stale"] = staleness(panel, catalogue, prim["start"])
    ctx["conv"] = conventions(panel, meta)
    ctx["splice"] = splice(panel, month_ends)
    ctx["plaus"] = plausibility(panel, catalogue, month_ends)
    ctx["extended"] = extended_sample(panel, prim["start"])
    ctx["long_atm"] = long_atm(panel)

    names = {"coverage": "coverage.csv", "two_sided": "two_sidedness.csv", "stale": "staleness.csv",
             "conv": "conventions.csv", "splice": "splice.csv", "plaus": "plausibility.csv"}
    for k, name in names.items():
        ctx[k].to_csv(out / name, index=False)
    samples = pd.concat({"primary_all_available": prim["avail"].all(axis=1),
                         "fenics_full_n_ccy": ctx["extended"]["full"]["count"],
                         "fenics_calibration_n_ccy": ctx["extended"]["calibration"]["count"]}, axis=1)
    samples.index.name = "month_end"
    samples.to_csv(out / "samples.csv")
    (out / "audit_report.md").write_text(report(ctx))
    print(f"audit written to {out}")
    print(f"primary start {d(prim['start'])} (strict {d(prim['strict_start'])}); binding {prim.get('binding')}")


if __name__ == "__main__":
    main()
