"""Download the currency portfolios of Lustig, Roussanov and Verdelhan (2011).

    .venv/bin/python scripts/acquire_verdelhan.py [--retrieval-date 2026-09-27]

Source: Adrien Verdelhan's data page (https://web.mit.edu/adrienv/www/Data.html),
file CurrencyPortfolios.xls ("Monthly Currency Excess Returns", portfolios of
Lustig, Roussanov and Verdelhan, 2011). The page states no terms of use, so
redistribution is not established and the file is stored in the git-ignored
data/private/verdelhan/<retrieval-date>/ with a manifest (URL, retrieval time,
HTTP Last-Modified header, size and SHA-256). The research design (section 10)
uses these portfolios for a sign and magnitude check of HML_FX.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
URL = "https://web.mit.edu/adrienv/www/CurrencyPortfolios.xls"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--retrieval-date", default=dt.date.today().isoformat())
    args = p.parse_args()
    out = ROOT / "data" / "private" / "verdelhan" / args.retrieval_date
    out.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        body, last_modified = r.read(), r.headers.get("Last-Modified")
    target = out / "CurrencyPortfolios.xls"
    target.write_bytes(body)
    manifest = {"url": URL, "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                "last_modified": last_modified, "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest(),
                "terms": "none stated on the source page; stored privately"}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"{target.name}: {len(body)} bytes, last modified {last_modified}, sha256 {manifest['sha256'][:12]}...")


if __name__ == "__main__":
    main()
