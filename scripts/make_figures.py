"""Figures for the paper (pgfplots code read by the paper with \\input).

Writes four figures to paper/figures/. Three are drawn from the portfolio-level
series published in data/public/fx_carry_portfolio_series/, so they can be
redrawn without an LSEG licence; the third (fig_theta.tex, Figure 3 of the
paper) combines a synthetic model calculation with the distribution across
currency-months of a reported statistic, and is skipped when those private
results are absent. No quote, per-currency series or calibrated parameter is
plotted.

1. fig_e1.tex (E1): (a) the monthly skew price C_skew and forward-discount
   spread FD of the primary sample, with their regime means; (b) C_skew
   against FD month by month, with rays of constant phi = C_skew/FD.
2. fig_oos.tex (E2): the cumulative difference between the squared
   out-of-sample errors of the historical mean and of the phi forecast,
   month by month over the 116 forecasts of Table 3 (expanding window with
   training from the extended sample), and the cumulative Clark-West adjusted
   difference, whose mean is the numerator of the Clark-West statistic.
3. fig_theta.tex (R10): (a) for the Heston model with rho = -0.3 of the
   model validation table (EURUSD-type market), the identified set of the
   third central moment, scaled by the model's variance to the power 3/2, as
   the fraction theta of the largest admissible exponents varies; (b) the
   cumulative distribution across currency-months of the breakdown point
   theta*, for the primary sample and, on the currency-months of the 5-delta
   test where theta* with the 5-delta prices is defined, without and with
   those prices.
4. fig_e3.tex (E3 and the long sample): (a) cumulative unhedged and
   at-the-money-hedged carry returns over 1995-2026; (b) the cumulative
   realised cost of 10-delta protection, HML^H - HML^U, against the cumulative
   ex-ante skew premium, term (iii), over 2013-2026.

Before writing, the script checks each plotted series against the numbers
printed in the paper (means, counts, test statistics, quartiles) and stops
on a mismatch.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from estimate_stage4 import FIRST_OOS, clark_west  # noqa: E402
from validate_moments_models import MARKETS, MODELS  # noqa: E402

from qef.fx.conventions import G10  # noqa: E402
from qef.fx.gk import CALL, PUT, strike_from_delta_smile  # noqa: E402
from qef.fx.identification import SharpMomentSet, boundary, skewness_sign_breakdown  # noqa: E402
from qef.fx.moments import middle_contracts  # noqa: E402
from qef.validation.models import model_smile, reference_moments, tail_exponent  # noqa: E402

PUB = ROOT / "data" / "public" / "fx_carry_portfolio_series"
RES = ROOT / "data" / "private" / "results"
OUT = ROOT / "paper" / "figures"
HIKE = 2022.0
TICKS = r"x tick label style={/pgf/number format/1000 sep=,font=\footnotesize}, y tick label style={font=\footnotesize}"
# Okabe-Ito palette (colour-blind safe); line styles also differ, so the figures remain legible in black and white.
# Vermillion: skew price and premium; blue: forward-discount spread and unhedged carry; green: hedged carry and realised hedge cost.
COLOURS = "\n".join(rf"\definecolor{{{n}}}{{HTML}}{{{h}}}" for n, h in [
    ("oiBlue", "0072B2"), ("oiVerm", "D55E00"), ("oiGreen", "009E73"), ("oiSky", "56B4E9"),
    ("oiOrange", "E69F00"), ("oiPurple", "CC79A7"), ("oiYellow", "F0E442")])


def check(value: float, printed: float, decimals: int, what: str) -> None:
    if round(value, decimals) != printed:
        raise SystemExit(f"{what}: computed {value:.6f}, paper prints {printed}")


def year(dates: pd.Series, months_ahead: int = 0) -> pd.Series:
    """Decimal year at the middle of the month, optionally moved to the month of the return."""
    m = dates.dt.to_period("M") + months_ahead
    return m.dt.year + (m.dt.month - 0.5) / 12


def coords(x, y) -> str:
    return " ".join(f"({a:.3f},{b:.3f})" for a, b in zip(x, y))


def coords_meta(x, y, m) -> str:
    return " ".join(f"({a:.3f},{b:.3f}) [{c:.3f}]" for a, b, c in zip(x, y, m))


def shade(x0: float, x1: float, y0: float, y1: float) -> str:
    return rf"\addplot[fill=gray!13, draw=none, forget plot] coordinates {{({x0},{y0}) ({x1},{y0}) ({x1},{y1}) ({x0},{y1})}} \closedcycle;"


def fig_e1() -> str:
    pr = pd.read_csv(PUB / "primary.csv", parse_dates=["date"])
    cs, fd, x = pr.C_skew_market_10d * 1e4, pr.FD * 1e4, year(pr.date)
    z = pr.regime == "zero_rate"
    for v, p, what in [(cs[z].mean(), 11.24, "C_skew zero-rate"), (cs[~z].mean(), 11.19, "C_skew hiking"),
                       (fd[z].mean(), 18.61, "FD zero-rate"), (fd[~z].mean(), 27.05, "FD hiking")]:
        check(v, p, 2, what)
    check(pr.phi_market_10d[z].mean(), 0.809, 3, "phi zero-rate")
    check(pr.phi_market_10d[~z].mean(), 0.454, 3, "phi hiking")
    if (int(z.sum()), int((~z).sum())) != (104, 56):
        raise SystemExit("regime counts differ from 104 and 56")
    # Statements of Section 5.1 on the two series.
    low = pr.date[fd < 10]
    check(cs.max(), 47.5, 1, "largest skew price")
    if (str(pr.date[cs.idxmax()].date()), len(low), str(low.min().date())[:7], str(low.max().date())[:7]) != ("2020-03-31", 21, "2020-04", "2022-01"):
        raise SystemExit("the statements on March 2020 or on the months with a spread below 10 bp do not hold")
    if set(pr.date[pr.phi_market_10d > 1.5].dt.year) != {2020, 2021}:
        raise SystemExit("a month with phi above 1.5 lies outside 2020-2021")
    x0, x1 = round(float(x.min()) - 0.2, 2), round(float(x.max()) + 0.2, 2)
    seg = lambda a, b, v, style: rf"\addplot[{style}, forget plot] coordinates {{({a},{v:.3f}) ({b},{v:.3f})}};"
    rays = []
    for r, label in [(0.25, "0.25"), (0.5, "0.5"), (1.0, "1"), (2.0, "2")]:
        xe = min(50.0, 50.0 / r)
        rays.append(rf"\addplot[gray!60, thin, forget plot] coordinates {{(0,0) ({xe:.2f},{r * xe:.2f})}};")
        lx, anchor = {0.25: (49.5, "north east"), 0.5: (49.5, "south east"), 1.0: (41.0, "south east"), 2.0: (19.0, "north west")}[r]
        rays.append(rf"\node[font=\scriptsize, text=gray, anchor={anchor}, inner sep=1pt] at (axis cs:{lx},{r * lx:.2f}) {{$\varphi={label}$}};")
    return rf"""% Written by scripts/make_figures.py from data/public/fx_carry_portfolio_series/primary.csv.
{COLOURS}
\begin{{tikzpicture}}
\begin{{axis}}[name=a, width=0.54\linewidth, height=6cm, xmin={x0}, xmax={x1}, ymin=0, ymax=56,
  xtick={{2014,2017,2020,2023,2026}}, {TICKS},
  ylabel={{\footnotesize bp per month}}, title={{\footnotesize (a) Skew price and forward-discount spread}},
  legend style={{font=\footnotesize, draw=none, fill=none, at={{(0.02,0.98)}}, anchor=north west}}, legend cell align=left]
{shade(HIKE, x1, 0, 56)}
\addplot[oiVerm, thick] coordinates {{{coords(x, cs)}}};
\addlegendentry{{Skew price $\Cskew_t$}}
\addplot[oiBlue, thick, densely dashed] coordinates {{{coords(x, fd)}}};
\addlegendentry{{Forward-discount spread $\FD_t$}}
{seg(x0, HIKE, cs[z].mean(), 'oiVerm, thin')}
{seg(HIKE, x1, cs[~z].mean(), 'oiVerm, thin')}
{seg(x0, HIKE, fd[z].mean(), 'oiBlue, thin')}
{seg(HIKE, x1, fd[~z].mean(), 'oiBlue, thin')}
\end{{axis}}
\begin{{axis}}[at={{(a.east)}}, anchor=west, xshift=1.2cm, width=0.36\linewidth, height=6cm, xmin=0, xmax=50, ymin=0, ymax=50,
  {TICKS}, xlabel={{\footnotesize $\FD_t$ (bp)}}, ylabel={{\footnotesize $\Cskew_t$ (bp)}},
  title={{\footnotesize (b) $\Cskew_t$ against $\FD_t$}},
  colormap name=viridis, point meta min=2013, point meta max=2027, colorbar,
  colorbar style={{width=0.18cm, ytick={{2014,2018,2022,2026}}, yticklabel style={{font=\scriptsize, /pgf/number format/1000 sep=}}}},
  legend style={{font=\footnotesize, draw=none, at={{(0.5,-0.3)}}, anchor=north, legend columns=2, /tikz/every even column/.append style={{column sep=0.3cm}}}}, legend cell align=left]
{chr(10).join(rays)}
\addplot[scatter, only marks, mark=*, mark size=1.3pt, point meta=explicit, scatter/use mapped color={{draw=black!55, fill=mapped color}}, forget plot] coordinates {{{coords_meta(fd[z], cs[z], x[z])}}};
\addplot[scatter, only marks, mark=square*, mark size=1.4pt, point meta=explicit, scatter/use mapped color={{draw=black!55, fill=mapped color}}, forget plot] coordinates {{{coords_meta(fd[~z], cs[~z], x[~z])}}};
\addlegendimage{{only marks, mark=*, mark size=1.3pt, draw=black!55, fill=gray!40}}
\addlegendentry{{Zero-rate regime}}
\addlegendimage{{only marks, mark=square*, mark size=1.4pt, draw=black!55, fill=gray!40}}
\addlegendentry{{Hiking regime}}
\end{{axis}}
\end{{tikzpicture}}
"""


def fig_oos() -> str:
    pr = pd.read_csv(PUB / "primary.csv", parse_dates=["date"])
    ex = pd.read_csv(PUB / "extended.csv", parse_dates=["date"])
    allm = pd.concat([ex[ex.status == "ok"][["date", "U", "phi"]],
                      pr[pr.status == "ok"][["date", "U", "phi_market_10d"]].rename(columns={"phi_market_10d": "phi"})]).sort_values("date")
    d, y, x = allm.date.to_numpy(), allm.U.to_numpy(), allm.phi.to_numpy()
    cw = clark_west(d, y, x)
    check(cw["cw_t"], 1.639, 3, "Clark-West statistic")
    if cw["n_forecasts"] != 116:
        raise SystemExit(f"{cw['n_forecasts']} forecasts, the paper 116")
    rows = []
    for i in np.where(d >= FIRST_OOS.to_datetime64())[0]:
        a, b = np.linalg.lstsq(np.column_stack([np.ones(i), x[:i]]), y[:i], rcond=None)[0]
        yhat, ybar = a + b * x[i], y[:i].mean()
        rows.append((pd.Timestamp(d[i]), (y[i] - ybar) ** 2 - (y[i] - yhat) ** 2,
                     (y[i] - ybar) ** 2 - ((y[i] - yhat) ** 2 - (ybar - yhat) ** 2)))
    o = pd.DataFrame(rows, columns=["date", "dsq", "f"])
    if len(o) != 116 or not np.isclose(o.f.mean(), cw["cw_mean"], rtol=1e-12, atol=0):
        raise SystemExit("the Clark-West differences do not reproduce the statistic's numerator")
    o["x"] = year(o.date, 1)
    # Statements of Section 5.2 on where the out-of-sample gain comes from (post hoc).
    ret_year = (o.date.dt.to_period("M") + 1).dt.year
    in_2020 = ret_year == 2020
    check(o.dsq[in_2020].sum() / o.dsq.sum(), 0.89, 2, "share of the squared-error difference from 2020")
    check(o.f[in_2020].sum() / o.f.sum(), 0.79, 2, "share of the Clark-West numerator from 2020")
    rest = o.f[~in_2020]
    check(rest.mean() / (rest.std(ddof=1) / np.sqrt(len(rest))), 1.24, 2, "Clark-West statistic without 2020")
    if len(rest) != 104 or o.dsq[ret_year >= 2022].sum() >= 0:
        raise SystemExit("the statements on the forecasts outside 2020 or from 2022 do not hold")
    o["cum"], o["cumf"] = o.dsq.cumsum() * 1e4, o.f.cumsum() * 1e4  # per cent squared
    lo = math.floor(min(o.cum.min(), o.cumf.min(), 0.0) - 1)
    hi = math.ceil(max(o.cum.max(), o.cumf.max()) * 1.5 + 1)
    x0, x1 = 2017.0, round(float(o.x.max()) + 0.1, 2)
    return rf"""% Written by scripts/make_figures.py from data/public/fx_carry_portfolio_series/ (primary.csv, extended.csv).
{COLOURS}
\begin{{tikzpicture}}
\begin{{axis}}[width=0.8\linewidth, height=5.2cm, xmin={x0}, xmax={x1}, ymin={lo}, ymax={hi}, {TICKS},
  ylabel={{\footnotesize Cumulative sum (per cent squared)}},
  legend style={{font=\footnotesize, draw=none, fill=none, at={{(0.02,0.98)}}, anchor=north west}}, legend cell align=left]
{shade(HIKE, x1, lo, hi)}
\addplot[gray, densely dotted, forget plot] coordinates {{({x0},0) ({x1},0)}};
\addplot[oiBlue, thick, const plot] coordinates {{{coords(o.x, o.cum)}}};
\addlegendentry{{Squared error of the historical mean less that of the $\varphi$ forecast}}
\addplot[oiPurple, thick, densely dashed, const plot] coordinates {{{coords(o.x, o.cumf)}}};
\addlegendentry{{Clark--West adjusted difference}}
\end{{axis}}
\end{{tikzpicture}}
"""


def heston_band():
    label = next(k for k in MODELS if k.startswith("Heston, ρ −0.3"))
    model, (ccy, market) = MODELS[label], MARKETS["EURUSD"]
    F, tau, usd = market.forward, market.tau, G10[ccy].usd_base
    smile, conv = model_smile(model, market), G10[ccy].delta
    K_min = strike_from_delta_smile(-0.10, F, tau, PUT, conv, smile, market.df_base)
    K_max = strike_from_delta_smile(0.10, F, tau, CALL, conv, smile, market.df_base)
    bd = boundary(F, tau, smile, K_min, K_max)
    m = reference_moments(model, tau, usd)
    var, k3 = m[1] - m[0] ** 2, m[2] - 3 * m[0] * m[1] + 2 * m[0] ** 3
    check(k3 / var ** 1.5, -0.638, 3, "Heston skewness")
    theta_star, _ = skewness_sign_breakdown(F, tau, smile, usd, bd)
    check(theta_star, 0.49, 2, "Heston theta*")
    lo = tail_exponent(model, tau, math.log(K_min / F), "lower")["inf"] / bd.gamma_bar
    hi = (tail_exponent(model, tau, math.log(K_max / F), "upper")["inf"] - 1) / (bd.eta_bar - 1)
    check(lo, 0.97, 2, "Heston lower-tail fraction")
    check(hi, 0.96, 2, "Heston upper-tail fraction")
    mid = middle_contracts(F, K_min, K_max, tau, smile, usd)["values"]
    band = []
    for th in np.round(np.arange(0.05, 1.0001, 0.025), 3):
        g, e = bd.exponents(min(float(th), 0.9999))
        a, b = SharpMomentSet(F, tau, smile, usd, bd, g, e, mid=mid).third_central_range()
        band.append((float(th), a / var ** 1.5, b / var ** 1.5))
    return band, k3 / var ** 1.5, theta_star, min(lo, hi)


def ecdf(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    return 100 * np.searchsorted(np.sort(values), grid, side="right") / values.size


def fig_theta() -> str:
    band, skew, theta_star, theta_model = heston_band()
    ident = pd.read_csv(RES / "2026-09-23" / "identification_1m.csv")
    ts = ident["theta_star"].dropna().to_numpy()
    if ts.size != 1440:
        raise SystemExit(f"theta* defined in {ts.size} currency-months, the paper uses 1,440")
    for q, p in zip([0.25, 0.5, 0.75], [0.49, 0.63, 0.80]):
        check(float(np.quantile(ts, q)), p, 2, f"theta* quartile {q}")
    w = pd.read_csv(RES / "2026-09-27" / "wing_test_1m.csv")
    same = w[w["theta_star_5"].notna()]
    if len(same) != 394:
        raise SystemExit(f"theta* with the 5-delta prices defined in {len(same)} currency-months, the paper 394")
    for col, printed in [("theta_star_5", (0.36, 0.45, 0.57)), ("theta_star_no5d", (0.53, 0.62, 0.73))]:
        for q, p in zip([0.25, 0.5, 0.75], printed):
            check(float(same[col].quantile(q)), p, 2, f"{col} quartile {q}")
    grid = np.round(np.linspace(0, 1, 201), 3)
    line = lambda v: coords(grid, ecdf(v, grid))
    clip = lambda v: min(max(v, -3.2), 1.6)
    th = [t for t, _, _ in band]
    lower, upper = coords(th, [clip(a) for _, a, _ in band]), coords(th, [clip(b) for _, _, b in band])
    return rf"""% Written by scripts/make_figures.py: (a) synthetic Heston model; (b) distributions across currency-months.
{COLOURS}
\begin{{tikzpicture}}
\begin{{axis}}[name=a, width=0.5\linewidth, height=5.4cm, xmin=0, xmax=1, ymin=-3.2, ymax=1.6, {TICKS},
  xlabel={{\footnotesize Fraction $\theta$ of the largest admissible exponents}}, ylabel={{\footnotesize Third central moment/variance$^{{3/2}}$}},
  title={{\footnotesize (a) Heston model: identified set}}]
\addplot[name path=lo, oiBlue, thick] coordinates {{{lower}}};
\addplot[name path=hi, oiBlue, thick] coordinates {{{upper}}};
\addplot[oiSky, fill opacity=0.35] fill between[of=lo and hi];
\addplot[gray, densely dotted] coordinates {{(0,0) (1,0)}};
\addplot[black, densely dashed] coordinates {{(0,{skew:.4f}) (1,{skew:.4f})}} node[pos=0.03, anchor=north west, font=\scriptsize] {{true value}};
\addplot[oiVerm, thick, densely dotted] coordinates {{({theta_star:.4f},-3.2) ({theta_star:.4f},1.6)}} node[pos=0.97, anchor=north west, font=\scriptsize, text=oiVerm] {{$\theta^*$}};
\addplot[oiGreen, thick, densely dotted] coordinates {{({theta_model:.4f},-3.2) ({theta_model:.4f},1.6)}} node[pos=0.97, anchor=north east, font=\scriptsize, text=oiGreen] {{model}};
\end{{axis}}
\begin{{axis}}[at={{(a.east)}}, anchor=west, xshift=1.6cm, width=0.5\linewidth, height=5.4cm, xmin=0, xmax=1, ymin=0, ymax=100, {TICKS},
  xlabel={{\footnotesize Fraction $\theta$ of the largest admissible exponents}}, ylabel={{\footnotesize Currency-months with $\theta^*\le\theta$ (\%)}},
  title={{\footnotesize (b) G10 currency-months: breakdown points}},
  legend style={{font=\footnotesize, draw=none, at={{(-0.3,-0.3)}}, anchor=north, legend columns=1}}, legend cell align=left]
\addplot[fill=oiYellow, fill opacity=0.3, draw=none, forget plot] coordinates {{(0.58,0) (0.97,0) (0.97,100) (0.58,100)}} \closedcycle;
\addplot[oiBlue, thick] coordinates {{{line(ts)}}};
\addlegendentry{{Primary sample, $10\D$ boundaries (1,440 currency-months)}}
\addplot[oiGreen, thick, densely dashed] coordinates {{{line(same['theta_star_no5d'].to_numpy())}}};
\addlegendentry{{$5\D$ test sample, without the $5\D$ prices (394)}}
\addplot[oiVerm, thick, densely dotted] coordinates {{{line(same['theta_star_5'].to_numpy())}}};
\addlegendentry{{$5\D$ test sample, with the $5\D$ prices as constraints (394)}}
\end{{axis}}
\end{{tikzpicture}}
"""


def fig_e3() -> str:
    la = pd.read_csv(PUB / "long_atm.csv", parse_dates=["date"])
    if len(la) != 379:
        raise SystemExit(f"long sample has {len(la)} months, the paper 379")
    check(la.U.mean() * 1e4, 37.2, 1, "long-sample unhedged mean")
    check(la.Hatm.mean() * 1e4, 16.3, 1, "long-sample hedged mean")
    pr = pd.read_csv(PUB / "primary.csv", parse_dates=["date"])
    pp = pr[pr.status == "ok"].copy()
    if len(pp) != 159:
        raise SystemExit(f"{len(pp)} primary return months, the paper 159")
    cost, prem = (pp.H10 - pp.U) * 1e4, pp.c3 * 1e4
    check(cost.mean(), -9.3, 1, "mean hedge cost")
    check(prem.mean(), -11.2, 1, "mean skew term")
    xl, xp = year(la.date, 1), year(pp.date, 1)
    cu, ch = la.U.cumsum() * 100, la.Hatm.cumsum() * 100
    ymin, ymax = math.floor(min(cu.min(), ch.min()) / 20 - 1) * 20, math.ceil(cu.max() / 20 + 1.5) * 20
    cmin = math.floor(min(cost.cumsum().min(), prem.cumsum().min()) / 250 - 1) * 250
    return rf"""% Written by scripts/make_figures.py from data/public/fx_carry_portfolio_series/ (long_atm.csv, primary.csv).
{COLOURS}
\begin{{tikzpicture}}
\begin{{axis}}[name=a, width=0.48\linewidth, height=5.6cm, xmin=1995, xmax=2027, ymin={ymin}, ymax={ymax}, {TICKS},
  ylabel={{\footnotesize Cumulative return (per cent)}}, title={{\footnotesize (a) Carry, 1995--2026}},
  legend style={{font=\footnotesize, draw=none, fill=none, at={{(0.02,0.98)}}, anchor=north west}}, legend cell align=left]
{shade(HIKE, 2027, ymin, ymax)}
\addplot[oiBlue, thick] coordinates {{{coords(xl, cu)}}};
\addlegendentry{{Unhedged}}
\addplot[oiGreen, thick, densely dashed] coordinates {{{coords(xl, ch)}}};
\addlegendentry{{At-the-money hedged}}
\end{{axis}}
\begin{{axis}}[at={{(a.east)}}, anchor=west, xshift=2.1cm, width=0.48\linewidth, height=5.6cm, xmin=2013.4, xmax=2026.9, ymin={cmin}, ymax=250, {TICKS},
  xtick={{2014,2016,2018,2020,2022,2024,2026}}, ylabel={{\footnotesize Cumulative (bp)}}, title={{\footnotesize (b) Cost of $10\D$ protection, 2013--2026}},
  legend style={{font=\footnotesize, draw=none, fill=none, at={{(0.98,0.98)}}, anchor=north east}}, legend cell align=left]
{shade(HIKE, 2026.9, cmin, 250)}
\addplot[oiVerm, thick, densely dashed] coordinates {{{coords(xp, prem.cumsum())}}};
\addlegendentry{{Skew premium, term (iii)}}
\addplot[oiGreen, thick] coordinates {{{coords(xp, cost.cumsum())}}};
\addlegendentry{{Realised cost}}
\end{{axis}}
\end{{tikzpicture}}
"""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    private = all((RES / d / f).exists() for d, f in [("2026-09-23", "identification_1m.csv"), ("2026-09-27", "wing_test_1m.csv")])
    for name, make in [("fig_e1.tex", fig_e1), ("fig_oos.tex", fig_oos), ("fig_theta.tex", fig_theta), ("fig_e3.tex", fig_e3)]:
        if make is fig_theta and not private:
            print(f"skipped paper/figures/{name}: it needs the per-currency results in data/private/results/")
            continue
        (OUT / name).write_text(make())
        print(f"wrote paper/figures/{name}")


if __name__ == "__main__":
    main()
