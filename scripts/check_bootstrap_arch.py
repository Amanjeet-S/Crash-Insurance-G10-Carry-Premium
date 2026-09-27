"""The E1 bootstrap intervals recomputed with arch 8.0.0 (Sheppard) as an independent implementation.

    .venv/bin/python scripts/check_bootstrap_arch.py

For each E1 series (primary sample; market and smile readings; 10Δ and 25Δ;
φ, C_skew and FD) and each regime, this compares the expected block length of
qef.stats.bootstrap with that of arch.bootstrap.optimal_block_length, and the
95% percentile interval of the regime difference in means from
qef.stats.bootstrap.bootstrap_distribution (as in estimate_e1.py, 9,999
draws, seed 20260924) with the interval from arch.bootstrap.StationaryBootstrap
applied to each regime separately with arch's own block length (9,999 draws,
seed 20260927). It also records whether any series falls in the degenerate case
in which the kernel estimate of the long-run variance is not positive. The
module docstring of qef.stats.bootstrap states the known differences between the
two implementations. Output: bootstrap_arch_check.md in the private results
folder, with statistics of pooled series only.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from arch.bootstrap import StationaryBootstrap
from arch.bootstrap import optimal_block_length as arch_block_length

from qef.stats.bootstrap import _block_length, bootstrap_distribution, optimal_block_length

ROOT = Path(__file__).resolve().parents[1]
REGIME_BREAK = pd.Timestamp("2022-01-01")
B = 9999


def arch_difference(a, b, seed=20260927):
    """Percentile interval of mean(b*) − mean(a*), each regime resampled by arch with its own block length."""
    draws = []
    for k, x in enumerate((a, b)):
        bl = float(arch_block_length(x)["stationary"].iloc[0])
        bs = StationaryBootstrap(bl, x, seed=seed + k)
        draws.append(bs.apply(np.mean, reps=B).ravel())
    return np.percentile(draws[1] - draws[0], [2.5, 97.5])


def main():
    d = ROOT / "data" / "private" / "results" / "2026-09-23"
    e1 = pd.read_csv(d / "e1_series.csv", parse_dates=["date"])
    e1 = e1[(e1["sample"] == "primary") & (e1.status == "ok")]
    rows = []
    for (reading, delta), g in e1.groupby(["reading", "delta"]):
        for col in ("phi", "C_skew", "FD"):
            s = g.dropna(subset=[col])
            a = s.loc[s.date < REGIME_BREAK, col].to_numpy()
            b = s.loc[s.date >= REGIME_BREAK, col].to_numpy()
            mine = np.percentile(bootstrap_distribution(lambda x, y: y.mean() - x.mean(), [a, b], B=B), [2.5, 97.5])
            theirs = arch_difference(a, b)
            row = {"reading": reading, "delta": delta, "variable": col}
            for name, x in (("zero_rate", a), ("hiking", b)):
                row[f"b_mine_{name}"] = optimal_block_length(x)
                row[f"b_arch_{name}"] = float(arch_block_length(x)["stationary"].iloc[0])
                row[f"b_mine_archconv_{name}"] = _block_length(x, arch_conventions=True)
            scale = 1e4 if col != "phi" else 1.0
            row.update({"lo_mine": scale * mine[0], "hi_mine": scale * mine[1],
                        "lo_arch": scale * theirs[0], "hi_arch": scale * theirs[1]})
            rows.append(row)
    res = pd.DataFrame(rows)
    degenerate = int(((res.filter(like="b_mine_") == 1.0).sum()).sum())
    fmt = lambda v: f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
    lines = ["# E1 bootstrap intervals with arch 8.0.0 (restricted summary)", "",
             "Interval bounds for C_skew and FD in bp per month; φ dimensionless. b: expected block length per regime.", "",
             "```", res.to_string(index=False, float_format=fmt), "```", "",
             f"Largest difference between interval bounds: φ {res.loc[res.variable == 'phi', ['lo_mine', 'hi_mine']].sub(res.loc[res.variable == 'phi', ['lo_arch', 'hi_arch']].values).abs().max().max():.3f}; "
             f"C_skew and FD {res.loc[res.variable != 'phi', ['lo_mine', 'hi_mine']].sub(res.loc[res.variable != 'phi', ['lo_arch', 'hi_arch']].values).abs().max().max():.2f} bp.",
             f"Block lengths equal to 1 (degenerate case): {degenerate}.",
             f"Do the two implementations agree on whether each interval excludes zero? "
             f"{bool(((res.lo_mine > 0) | (res.hi_mine < 0)).eq((res.lo_arch > 0) | (res.hi_arch < 0)).all())}", ""]
    (d / "bootstrap_arch_check.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
