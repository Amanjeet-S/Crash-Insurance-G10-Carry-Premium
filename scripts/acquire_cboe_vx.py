"""Download Cboe VIX futures (VX) daily settlement files, one per contract.

Source: Cboe Futures Exchange historical data (cboe.com, US futures market
statistics). Contracts expiring before 2013 are in the archive
(`CFE_<month code><yy>_VX.csv`); later contracts are named by expiry date
(`VX_<yyyy-mm-dd>.csv`). The expiry date follows the contract rule: the
Wednesday 30 days before the third Friday of the following month. If no file
exists for that date, the neighbouring business days are tried.

Redistribution terms are not established, so files are stored in the
git-ignored `data/private/cboe/<retrieval-date>/` with a manifest (URL,
SHA-256, rows, retrieval time).
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = "https://cdn.cboe.com/resources/futures/archive/volume-and-price/CFE_{code}{yy}_VX.csv"
CURRENT = "https://cdn.cboe.com/data/us/futures/market_statistics/historical_data/VX/VX_{date}.csv"
MONTH_CODES = "FGHJKMNQUVXZ"


def rule_expiry(year: int, month: int) -> dt.date:
    """Wednesday 30 days before the third Friday of the following month."""
    y, m = (year + 1, 1) if month == 12 else (year, month + 1)
    first = dt.date(y, m, 1)
    third_friday = first + dt.timedelta(days=(4 - first.weekday()) % 7 + 14)
    return third_friday - dt.timedelta(days=30)


def fetch(url: str) -> bytes | None:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.read() if r.status == 200 else None
    except Exception:
        return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--retrieval-date", default=dt.date.today().isoformat())
    p.add_argument("--first", default="2007-01")
    p.add_argument("--last", default="2026-12")
    args = p.parse_args()
    out = ROOT / "data" / "private" / "cboe" / args.retrieval_date
    (out / "raw").mkdir(parents=True, exist_ok=True)
    y0, m0 = map(int, args.first.split("-"))
    y1, m1 = map(int, args.last.split("-"))
    manifest, missing = [], []
    y, m = y0, m0
    while (y, m) <= (y1, m1):
        candidates = []
        if y <= 2012:
            candidates.append(ARCHIVE.format(code=MONTH_CODES[m - 1], yy=f"{y % 100:02d}"))
        base = rule_expiry(y, m)
        for shift in (0, -1, 1, -2, -7):
            candidates.append(CURRENT.format(date=(base + dt.timedelta(days=shift)).isoformat()))
        data, url = None, None
        for url in candidates:
            data = fetch(url)
            if data and data[:10].lower().startswith(b"trade date"):
                break
            data = None
            time.sleep(0.2)
        if data is None:
            missing.append(f"{y}-{m:02d}")
        else:
            name = url.rsplit("/", 1)[1]
            (out / "raw" / name).write_bytes(data)
            manifest.append({"contract_month": f"{y}-{m:02d}", "file": name, "url": url,
                             "sha256": hashlib.sha256(data).hexdigest(), "rows": data.count(b"\n") - 1,
                             "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")})
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
        time.sleep(0.2)
    (out / "manifest.json").write_text(json.dumps({"source": "Cboe Futures Exchange historical data (VX)",
                                                   "files": manifest, "missing_contract_months": missing}, indent=1))
    print(f"{len(manifest)} contracts saved; missing: {missing}")


if __name__ == "__main__":
    main()
