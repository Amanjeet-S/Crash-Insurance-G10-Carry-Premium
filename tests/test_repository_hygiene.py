"""Repository hygiene: licensed data are never tracked.

LSEG data are licensed for individual study and research and may not be
redistributed (README). These checks fail if anything under data/private/, a
data extract outside the attributed public snapshots in data/public/, or a
notebook with stored outputs is tracked by git.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA_SUFFIXES = {".csv", ".tsv", ".xlsx", ".xls", ".parquet", ".feather", ".h5", ".pkl", ".json", ".jsonl"}


def _tracked() -> list[Path]:
    if shutil.which("git") is None or not (ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
    return [Path(p) for p in out.splitlines()]


def test_no_private_paths_or_data_extracts_tracked():
    bad = [str(p) for p in _tracked()
           if p.parts[:2] == ("data", "private")
           or (p.suffix.lower() in DATA_SUFFIXES and p.parts[:2] != ("data", "public"))]
    assert not bad, f"data files are tracked: {bad}"


def test_tracked_notebooks_have_no_outputs():
    for p in (q for q in _tracked() if q.suffix == ".ipynb"):
        nb = json.loads((ROOT / p).read_text())
        stored = [c for c in nb.get("cells", []) if c.get("cell_type") == "code" and c.get("outputs")]
        assert not stored, f"{p} has stored outputs"
