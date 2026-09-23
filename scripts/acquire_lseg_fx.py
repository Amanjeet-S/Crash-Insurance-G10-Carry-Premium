"""Acquire the daily LSEG FX panel for project 01 into the private raw layer.

Run in the LSEG environment (requirements-lseg.lock) with LSEG Workspace open:

    .venv-lseg/bin/python scripts/acquire_lseg_fx.py

Output, under data/private/lseg/<retrieval date>/:

    raw/<block>/<RIC>.csv   returned tables, unmodified, full precision, date index
    requests.jsonl          one line per request
    metadata.csv            search metadata per RIC: title, currencies, underlying, scaling
    manifest.json           SHA-256, rows, first and last date of every file

A history request returns all default fields. If the first returned date lies
more than a week after the requested start, the earlier range is requested
again, until a request returns nothing or reaches the start. Existing raw files
are kept unless --overwrite is given, so an interrupted run can be resumed.
The App Key is read from ~/.lseg/app_key and is never printed or written.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
import time
import warnings
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qef.data.panel import instrument_catalogue, raw_path  # noqa: E402

import lseg.data as ld  # noqa: E402
from lseg.data import discovery  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

INTERVAL = "daily"
KEY_FILE = Path("~/.lseg/app_key").expanduser()
# Search fields that identify the underlying pair, contributor and point scaling.
SEARCH_FIELDS = [
    "RIC", "DocumentTitle", "CommonName", "FirstCurrency", "SecondCurrency", "ConcatCcyCode",
    "UnderlyingQuoteRIC", "Leg1UnderlyingRIC", "Leg1UnderlyingName", "ContributorCommonName",
    "ContributorCode", "DTSimpleType", "RCSAssetCategoryLeaf", "RCSTermMaturityName",
    "ScalingFactor", "PriceType", "RCSQuoteTypeLeaf", "IsPrimaryQuote", "HasTimeSeries",
    "ExchangeCode",
]


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


class RequestLog:
    def __init__(self, path: Path):
        self.path = path

    def write(self, record: dict) -> None:
        record = {"library_version": ld.__version__, **record}
        with self.path.open("a") as f:
            f.write(json.dumps(record, default=str) + "\n")


def request_history(ric: str, start: str, end: str, log: RequestLog, block: str, retries: int = 2):
    """One get_history call, retried on exceptions, each attempt logged."""
    for attempt in range(retries + 1):
        record = {"kind": "history", "ric": ric, "block": block, "fields": None, "interval": INTERVAL,
                  "start": start, "end": end, "request_utc": utc_now(), "attempt": attempt + 1}
        try:
            df = ld.get_history(universe=ric, interval=INTERVAL, start=start, end=end)
        except Exception as exc:  # the library raises LDError for unknown or unentitled RICs
            record.update(rows=0, error=f"{type(exc).__name__}: {str(exc)[:500]}")
            log.write(record)
            if "not found" in str(exc).lower() or attempt == retries:
                return None, record["error"]
            time.sleep(2.0 * (attempt + 1))
            continue
        n = 0 if df is None else len(df)
        record.update(
            rows=n,
            columns=[] if df is None else [str(c) for c in df.columns],
            first_date=str(df.index.min().date()) if n else None,
            last_date=str(df.index.max().date()) if n else None,
            error=None if n else "no rows returned",
        )
        log.write(record)
        return (df if n else None), record["error"]
    return None, "unreachable"


def fetch_full_history(ric: str, start: str, end: str, log: RequestLog, block: str, pause: float):
    """Request [start, end] and extend backwards while the response is truncated."""
    chunks, errors = [], []
    lo, hi = pd.Timestamp(start), pd.Timestamp(end)
    while True:
        df, err = request_history(ric, str(lo.date()), str(hi.date()), log, block)
        time.sleep(pause)
        if df is None:
            if err and not chunks:
                errors.append(err)
            break
        chunks.append(df)
        first = pd.Timestamp(df.index.min())
        if first <= lo + pd.Timedelta(days=7):
            break
        hi = first - pd.Timedelta(days=1)
    if not chunks:
        return None, errors
    # Chunks do not overlap by construction; keep the column order of the latest chunk.
    out = pd.concat(chunks[::-1]).sort_index()
    out.index.name = "Date"
    return out, errors


def search_metadata(rics: list[str], log: RequestLog, batch: int = 20) -> pd.DataFrame:
    rows = []
    for i in range(0, len(rics), batch):
        group = rics[i : i + batch]
        filt = " or ".join(f"RIC eq '{r}'" for r in group)
        record = {"kind": "search", "rics": group, "view": "SEARCH_ALL", "select": SEARCH_FIELDS,
                  "request_utc": utc_now()}
        try:
            df = discovery.search(view=discovery.Views.SEARCH_ALL, filter=filt,
                                  select=",".join(SEARCH_FIELDS), top=100)
            record.update(rows=len(df), error=None)
            rows.append(df)
        except Exception as exc:
            record.update(rows=0, error=f"{type(exc).__name__}: {str(exc)[:500]}")
        log.write(record)
        time.sleep(0.2)
    found = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=SEARCH_FIELDS)
    found = found.drop_duplicates(subset="RIC")
    # Some search fields are lists (for example FirstCurrency); store them as "a|b".
    for col in found.columns:
        found[col] = found[col].map(lambda v: "|".join(map(str, v)) if isinstance(v, list) else v)
    out = pd.DataFrame({"requested_ric": rics}).merge(found, how="left", left_on="requested_ric", right_on="RIC")
    out.insert(1, "found", out["RIC"].notna())
    return out


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_manifest(out: Path, args) -> dict:
    files = []
    for path in sorted(p for p in out.rglob("*") if p.is_file() and p.name != "manifest.json"):
        entry = {"file": str(path.relative_to(out)), "sha256": sha256(path)}
        if path.parent.parent.name == "raw":
            df = pd.read_csv(path, index_col=0, parse_dates=True)
            entry.update(rows=len(df), first_date=str(df.index.min().date()) if len(df) else None,
                         last_date=str(df.index.max().date()) if len(df) else None)
        files.append(entry)
    manifest = {
        "retrieval_date": args.retrieval_date, "requested_start": args.start, "requested_end": args.end,
        "interval": INTERVAL, "library": f"lseg-data {ld.__version__}",
        "python": sys.version.split()[0], "written_utc": utc_now(), "files": files,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--start", default="1995-01-01")
    p.add_argument("--end", default="2026-09-22")
    p.add_argument("--retrieval-date", default=dt.date.today().isoformat())
    p.add_argument("--out-root", default=str(ROOT / "data" / "private" / "lseg"))
    p.add_argument("--only", nargs="*", help="restrict to these RICs")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--pause", type=float, default=0.2, help="seconds between requests")
    p.add_argument("--skip-search", action="store_true")
    args = p.parse_args()

    out = Path(args.out_root) / args.retrieval_date
    out.mkdir(parents=True, exist_ok=True)
    log = RequestLog(out / "requests.jsonl")
    catalogue = instrument_catalogue()
    if args.only:
        catalogue = [i for i in catalogue if i.ric in set(args.only)]

    ld.open_session(app_key=KEY_FILE.read_text().strip())
    try:
        failed = {}
        for k, inst in enumerate(catalogue, 1):
            path = raw_path(out, inst)
            if path.exists() and not args.overwrite:
                continue
            df, errors = fetch_full_history(inst.ric, args.start, args.end, log, inst.block, args.pause)
            if df is None:
                failed[inst.ric] = errors[-1] if errors else "no rows returned"
                print(f"[{k}/{len(catalogue)}] {inst.ric}: no data ({failed[inst.ric][:80]})")
                continue
            path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(path)  # pandas writes floats at full (shortest round-trip) precision
            print(f"[{k}/{len(catalogue)}] {inst.ric}: {len(df)} rows, {df.index.min().date()} to {df.index.max().date()}")

        if not args.skip_search:
            meta = search_metadata([i.ric for i in instrument_catalogue()], log)
            meta.to_csv(out / "metadata.csv", index=False)
            print(f"metadata: {int(meta['found'].sum())} of {len(meta)} RICs found by search")
    finally:
        ld.close_session()

    manifest = write_manifest(out, args)
    n_raw = sum(1 for f in manifest["files"] if f["file"].startswith("raw/"))
    print(f"{n_raw} raw files; {len(failed)} RICs without data this run: {sorted(failed)}")


if __name__ == "__main__":
    main()
