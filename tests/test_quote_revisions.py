import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

from qef.data.panel import Instrument

SPEC = importlib.util.spec_from_file_location(
    "check_quote_revisions", Path(__file__).resolve().parents[1] / "scripts" / "check_quote_revisions.py")
rev = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rev)

END = pd.Timestamp("2021-03-31")


def frame(seed=0):
    # Synthetic values on New York business days; the month-ends are 29 Jan, 26 Feb and 31 Mar 2021.
    days = pd.bdate_range("2021-01-04", END)
    rng = np.random.default_rng(seed)
    bid = 1.2 + rng.normal(0, 1e-3, len(days)).cumsum()
    return pd.DataFrame({"BID": bid, "ASK": bid + 2e-4, "MID_PRICE": bid + 1e-4, "MATUR_DATE": "2021-05-01"},
                        index=pd.DatetimeIndex(days, name="Date"))


def write(root: Path, inst: Instrument, df: pd.DataFrame):
    path = root / "raw" / inst.block / f"{inst.ric}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path)


def test_identical_retrievals_show_no_revision(tmp_path):
    inst = Instrument("EUR1MO=", "vol_atm", "EUR", "1M", "atm")
    for d in ("old", "new"):
        write(tmp_path / d, inst, frame())
    r = rev.compare(inst, [tmp_path / "old"], tmp_path / "new", END)
    assert r["identical_bytes"] and r["quote_changed"] == r["other_changed"] == 0
    assert r["month_ends"] == 3 and r["month_ends_revised"] == 0


def test_revisions_are_counted_by_kind_and_at_month_ends(tmp_path):
    inst = Instrument("EUR1MRR=", "vol_rr25", "EUR", "1M", "rr25")
    old = frame()
    new = old.copy()
    new.loc["2021-02-26", "BID"] += 1e-4            # revised bid on a month-end
    new.loc["2021-01-12", "MID_PRICE"] = np.nan     # a vendor field that disappeared
    new.loc["2021-03-15", "ASK"] = old.loc["2021-03-15", "ASK"]  # unchanged
    old.loc["2021-03-10", "ASK"] = np.nan           # an ask that appeared in the new retrieval
    new = pd.concat([new, frame().loc[["2021-03-31"]].rename(lambda _: pd.Timestamp("2021-01-02"))]).sort_index()
    write(tmp_path / "old", inst, old)
    write(tmp_path / "new", inst, new)
    r = rev.compare(inst, [tmp_path / "old"], tmp_path / "new", END)
    assert not r["identical_bytes"]
    assert (r["quote_changed"], r["quote_appeared"], r["quote_disappeared"]) == (1, 1, 0)
    assert (r["other_changed"], r["other_appeared"], r["other_disappeared"]) == (0, 0, 1)
    assert (r["dates_only_new"], r["dates_only_old"]) == (1, 0)
    assert r["month_ends_revised"] == 1
    assert r["quote_first_revised"] == pd.Timestamp("2021-02-26") and r["quote_last_revised"] == pd.Timestamp("2021-03-10")


def test_instrument_is_taken_from_the_first_old_retrieval_that_holds_it(tmp_path):
    inst = Instrument("EUR3M=", "forward", "EUR", "3M", "points")
    write(tmp_path / "second", inst, frame(1))
    write(tmp_path / "new", inst, frame(1))
    r = rev.compare(inst, [tmp_path / "first", tmp_path / "second"], tmp_path / "new", END)
    assert r["old_retrieval"] == "second" and r["identical_bytes"]
    missing = rev.compare(Instrument("GBP3M=", "forward", "GBP", "3M", "points"), [tmp_path / "first"], tmp_path / "new", END)
    assert missing["status"] == "missing_both"
