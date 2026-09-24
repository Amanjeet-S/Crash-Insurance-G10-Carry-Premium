import pandas as pd

from qef.data.vix import load_vx, price_on


def _write(path, rows):
    pd.DataFrame(rows, columns=["Trade Date", "Futures", "Settle"]).to_csv(path, index=False)


def test_archive_rescaling_and_expiry_keys(tmp_path):
    # Archive contract spanning the rescaling: old scale (x10) before 26 March 2007.
    _write(tmp_path / "CFE_K07_VX.csv", [["03/22/2007", "K (May 07)", 140.0], ["03/23/2007", "K (May 07)", 141.0],
                                          ["03/26/2007", "K (May 07)", 14.2], ["05/15/2007", "K (May 07)", 15.0]])
    # Archive contract that expired before the rescaling: entirely at the old scale.
    _write(tmp_path / "CFE_H07_VX.csv", [["02/20/2007", "H (Mar 07)", 120.0], ["03/20/2007", "H (Mar 07)", 118.0]])
    # Current-format file: keyed by the expiry in its name, even while still trading.
    _write(tmp_path / "VX_2026-10-21.csv", [["2026-09-22", "V (Oct 2026)", 18.4]])
    contracts, rescaled = load_vx(tmp_path)
    assert ("CFE_K07_VX.csv", "2007-03-26") in rescaled
    k07 = contracts[pd.Timestamp("2007-05-15")]
    assert k07.loc["2007-03-23"] == 14.1 and k07.loc["2007-03-26"] == 14.2
    assert contracts[pd.Timestamp("2007-03-20")].loc["2007-03-20"] == 11.8
    assert pd.Timestamp("2026-10-21") in contracts
    assert price_on(contracts[pd.Timestamp("2026-10-21")], pd.Timestamp("2026-09-24")) == 18.4
