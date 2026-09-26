"""Stage 1 quote-revision check: compare two LSEG retrievals of the same series.

    .venv/bin/python scripts/check_quote_revisions.py [--new 2026-09-26] [--old 2026-09-23 2026-09-24]

The data plan (storage and provenance) separates the retrieval date from the
observation date because a current retrieval need not reproduce what an
earlier one returned. This script measures that directly. Every instrument in
the catalogue is taken from the first ``--old`` retrieval that holds it and
compared with the ``--new`` retrieval over the dates both requested (up to
``--end``). Both retrievals store the returned tables unmodified at full
precision, so an unrevised series gives identical values.

For each instrument the script records whether the files are byte-identical;
otherwise the dates only one retrieval returns, and, on common dates, the
cells that differ, split into bid and ask (the fields the estimation uses) and
all other fields, and into changed values, values that appeared and values
that disappeared. The month-end rule of the research design (section 4) is
then applied to both retrievals, and a month-end counts as revised if its
status, source date, bid or ask differs. Month-ends are counted over the whole
history and over the primary sample, which starts at the 31 May 2013
month-end.

Writes data/private/audit/<new>/revisions.csv (one row per instrument) and
revisions_report.md (pooled counts). Both are LSEG-derived and stay in
data/private/; the console shows pooled counts only.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from qef.data.panel import instrument_catalogue, ny_month_ends, raw_path, read_raw, sample_month_ends

ROOT = Path(__file__).resolve().parents[1]
WINDOW = 5  # business days for month-end substitution (research design, section 4)
PRIMARY_START = pd.Timestamp("2013-05-31")  # private audit, primary-sample rule
QUOTE_FIELDS = ("BID", "ASK")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load(path: Path, end: pd.Timestamp) -> pd.DataFrame:
    frame = read_raw(path)
    frame = frame[~frame.index.duplicated(keep="last")].sort_index()
    return frame[frame.index <= end]


def cell_changes(a: pd.DataFrame, b: pd.DataFrame, cols) -> dict:
    """Counts of changed, appeared and disappeared values on the common dates.

    Values are compared as numbers where both parse, otherwise as strings, so
    that date-valued fields are covered too.
    """
    out = {"cells": 0, "changed": 0, "appeared": 0, "disappeared": 0, "first": pd.NaT, "last": pd.NaT}
    dates = a.index.intersection(b.index)
    touched = pd.Series(False, index=dates)
    for c in cols:
        x, y = a.loc[dates, c], b.loc[dates, c]
        xn, yn = pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")
        numeric = (xn.notna() | x.isna()) & (yn.notna() | y.isna())
        x_miss, y_miss = x.isna(), y.isna()
        both = ~x_miss & ~y_miss
        differ = np.where(numeric, xn.to_numpy() != yn.to_numpy(), x.astype(str).to_numpy() != y.astype(str).to_numpy())
        changed = both & pd.Series(differ, index=dates)
        appeared, disappeared = x_miss & ~y_miss, ~x_miss & y_miss
        out["cells"] += int((~x_miss | ~y_miss).sum())
        out["changed"] += int(changed.sum())
        out["appeared"] += int(appeared.sum())
        out["disappeared"] += int(disappeared.sum())
        touched |= changed | appeared | disappeared
    if touched.any():
        out["first"], out["last"] = touched[touched].index.min(), touched[touched].index.max()
    return out


def month_end_changes(a: pd.DataFrame, b: pd.DataFrame, end: pd.Timestamp) -> dict:
    """Month-ends whose selected quote (status, source date, bid or ask) differs between retrievals."""
    if not all(f in a and f in b for f in QUOTE_FIELDS):
        return {"month_ends": 0, "month_ends_revised": 0, "primary_month_ends": 0, "primary_month_ends_revised": 0}
    start = min(a.index.min(), b.index.min())
    mes = ny_month_ends(start, end)
    sa = sample_month_ends(a, mes, WINDOW)
    sb = sample_month_ends(b, mes, WINDOW)
    same_status = sa["status"].to_numpy() == sb["status"].to_numpy()
    same_date = (sa["source_date"].to_numpy() == sb["source_date"].to_numpy()) | (sa["source_date"].isna() & sb["source_date"].isna()).to_numpy()
    same_val = np.ones(len(mes), dtype=bool)
    for f in ("bid", "ask"):
        x, y = sa[f].to_numpy(dtype=float), sb[f].to_numpy(dtype=float)
        same_val &= (x == y) | (np.isnan(x) & np.isnan(y))
    revised = ~(same_status & same_date & same_val)
    primary = mes >= PRIMARY_START
    return {"month_ends": len(mes), "month_ends_revised": int(revised.sum()),
            "primary_month_ends": int(primary.sum()), "primary_month_ends_revised": int((revised & primary).sum())}


def compare(inst, old_roots, new_root: Path, end: pd.Timestamp) -> dict:
    row = {"ric": inst.ric, "block": inst.block, "currency": inst.currency, "tenor": inst.tenor}
    old = next((raw_path(r, inst) for r in old_roots if raw_path(r, inst).exists()), None)
    new = raw_path(new_root, inst)
    row["old_retrieval"] = old.parents[2].name if old else ""
    if old is None or not new.exists():
        row["status"] = "missing_old" if old is None and new.exists() else "missing_new" if old else "missing_both"
        return row
    row["identical_bytes"] = sha256(old) == sha256(new)
    a, b = load(old, end), load(new, end)
    row.update({"rows_old": len(a), "rows_new": len(b),
                "dates_only_old": len(a.index.difference(b.index)), "dates_only_new": len(b.index.difference(a.index)),
                "columns_only_old": len(a.columns.difference(b.columns)),
                "columns_only_new": len(b.columns.difference(a.columns))})
    common = [c for c in a.columns if c in b.columns]
    quote = [c for c in common if c in QUOTE_FIELDS]
    other = [c for c in common if c not in QUOTE_FIELDS]
    for name, cols in (("quote", quote), ("other", other)):
        ch = cell_changes(a, b, cols)
        row.update({f"{name}_cells": ch["cells"], f"{name}_changed": ch["changed"],
                    f"{name}_appeared": ch["appeared"], f"{name}_disappeared": ch["disappeared"],
                    f"{name}_first_revised": ch["first"], f"{name}_last_revised": ch["last"]})
    row.update(month_end_changes(a, b, end))
    row["status"] = "compared"
    return row


def report(res: pd.DataFrame, args) -> str:
    comp = res[res.status == "compared"]
    L = [f"# Quote revisions: retrieval of {args.new} against {', '.join(args.old)}", "",
         "Restricted: LSEG-derived counts and dates. Pooled counts may be published; values may not.", "",
         f"Dates compared: up to {args.end}. Instruments in the catalogue: {len(res)}; compared: {len(comp)}; "
         f"missing from a retrieval: {int((res.status != 'compared').sum())}.", ""]
    for status, n in res.status.value_counts().items():
        L.append(f"- {status}: {n}")
    L += ["", "| Block | Compared | Byte-identical | Any change | Quote cells | Quote cells revised | Other cells revised "
          "| Dates added | Dates removed | Month-ends revised (all) | Month-ends revised (primary) |",
          "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    comp = comp.assign(any_change=(comp.quote_changed + comp.quote_appeared + comp.quote_disappeared
                                   + comp.other_changed + comp.other_appeared + comp.other_disappeared
                                   + comp.dates_only_old + comp.dates_only_new + comp.columns_only_old
                                   + comp.columns_only_new) > 0)
    groups = [(b, g) for b, g in comp.groupby("block")] + [("all", comp)]
    for b, g in groups:
        q_rev = int((g.quote_changed + g.quote_appeared + g.quote_disappeared).sum())
        o_rev = int((g.other_changed + g.other_appeared + g.other_disappeared).sum())
        L.append(f"| {b} | {len(g)} | {int(g.identical_bytes.sum())} | {int(g.any_change.sum())} | {int(g.quote_cells.sum())} "
                 f"| {q_rev} | {o_rev} | {int(g.dates_only_new.sum())} | {int(g.dates_only_old.sum())} "
                 f"| {int(g.month_ends_revised.sum())} of {int(g.month_ends.sum())} "
                 f"| {int(g.primary_month_ends_revised.sum())} of {int(g.primary_month_ends.sum())} |")
    changed = comp[comp.any_change]
    if len(changed):
        L += ["", "Instruments with any change (first and last revised date of bid or ask, then of other fields):", ""]
        for r in changed.itertuples():
            L.append(f"- `{r.ric}`: quote cells revised {r.quote_changed + r.quote_appeared + r.quote_disappeared} "
                     f"({r.quote_first_revised} to {r.quote_last_revised}); other cells revised "
                     f"{r.other_changed + r.other_appeared + r.other_disappeared}; dates added {r.dates_only_new}, "
                     f"removed {r.dates_only_old}; month-ends revised {r.month_ends_revised} "
                     f"(primary {r.primary_month_ends_revised}).")
    return "\n".join(L) + "\n"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--new", default="2026-09-26")
    p.add_argument("--old", nargs="+", default=["2026-09-23", "2026-09-24"])
    p.add_argument("--end", default="2026-09-22", help="last date requested by both retrievals")
    p.add_argument("--lseg-root", default=str(ROOT / "data" / "private" / "lseg"))
    p.add_argument("--out-root", default=str(ROOT / "data" / "private" / "audit"))
    args = p.parse_args()

    lseg = Path(args.lseg_root)
    old_roots = [lseg / d for d in args.old]
    end = pd.Timestamp(args.end)
    res = pd.DataFrame([compare(i, old_roots, lseg / args.new, end) for i in instrument_catalogue()])
    out = Path(args.out_root) / args.new
    out.mkdir(parents=True, exist_ok=True)
    res.to_csv(out / "revisions.csv", index=False)
    text = report(res, args)
    (out / "revisions_report.md").write_text(text)
    comp = res[res.status == "compared"]
    q = int((comp.quote_changed + comp.quote_appeared + comp.quote_disappeared).sum())
    print(f"{len(comp)} instruments compared, {int(comp.identical_bytes.sum())} byte-identical; "
          f"bid and ask cells revised: {q} of {int(comp.quote_cells.sum())}; "
          f"month-ends revised: {int(comp.month_ends_revised.sum())} of {int(comp.month_ends.sum())} "
          f"(primary sample {int(comp.primary_month_ends_revised.sum())} of {int(comp.primary_month_ends.sum())})")
    print(f"written to {out}")


if __name__ == "__main__":
    main()
