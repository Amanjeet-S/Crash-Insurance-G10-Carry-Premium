"""Timings of the C++ kernel against the Python reference, written to cpp/README.md.

Builds the kernel with cpp/build.sh, times each kernel and its Python
reference on a realistic workload, measures parity on the same inputs, and
replaces the block between the benchmark markers of cpp/README.md with the
build command, the coverage, the parity tolerances, the timings and the
parity statistics. All inputs are synthetic and seeded (seed 20260928); no
licensed data are read, so every output may be published.

Workloads, modelled on the 1,440 currency-months of the one-month moment
estimation (nine currencies, 160 months):

- premium and delta: 1,000,000 options on the grid of tests/test_cpp_kernel.py
  but with τ exactly one or three months (the tests also draw 28 to 31 or 89
  to 92 days over 365), σ on [0.03, 0.4], strikes from 5-delta to
  95-delta, both option types and D_b on [0.97, 1.01]; the reference is
  vectorised NumPy, so this compares array code with array code;
- SABR: one smile at 1,000,000 strikes (both vectorised), and 1,440 smiles at
  513 strikes each (the reference takes scalar parameters, so it loops over
  smiles, while the kernel takes all smiles in one call);
- strikes: the 10-delta put and call strikes of the 1,440 smiles at σ = α in
  each currency's convention, 1,280 pips and 1,600 premium-adjusted; the
  reference is scalar and loops;
- middle part: 1,440 SABR smiles, τ = 1/12, α on [0.05, 0.20], ρ on
  [−0.5, 0.5], ν on [0.2, 2.0], F a rounded public level per currency times
  exp(U(−0.2, 0.2)), D_b = 1, K_min and K_max the smile's 10-delta put and
  call strikes (strike_from_delta_smile, as in implied_moment_intervals).
  Adaptive: the reference's doubling from n0 = 32 to a relative change below
  1e-10; fixed: every smile at the median converged n, against the
  reference's inner integrate(n), which this script reproduces line for
  line and checks bitwise against middle_contracts.

Rules fixed before computing the results: each timing is the median of five
repetitions of one call (three for the reference's Python loops over smiles
or strikes), after one untimed warm-up call, single-threaded, on an otherwise
idle machine; the speed-up is the reference's median over the kernel's.
Parity statistics are pooled over each workload: the largest relative error,
the largest normwise error where the tests use one (premium, cubic contract),
and the share of results that are bitwise equal.
"""

from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import datetime  # noqa: E402
import math  # noqa: E402
import platform  # noqa: E402
import statistics  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import scipy  # noqa: E402
from scipy import special  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qef import kernel as kk  # noqa: E402
from qef.fx import gk, moments, sabr  # noqa: E402
from qef.fx.conventions import G10  # noqa: E402
from qef.fx.smile import SabrSmile  # noqa: E402

README = ROOT / "cpp" / "README.md"
BEGIN, END = "<!-- benchmark:begin -->", "<!-- benchmark:end -->"
SEED = 20260928
LEVELS = {"EUR": 1.10, "GBP": 1.30, "AUD": 0.70, "NZD": 0.65, "JPY": 130.0, "CHF": 0.95, "CAD": 1.30,
          "NOK": 10.0, "SEK": 10.0}
MONTHS = 160
REPS, REPS_LOOP = 5, 3
CONVENTIONS = [gk.DeltaConvention(s, p) for s in (True, False) for p in (False, True)]

COVERAGE = [
    ("`forward_premium`", "`qef.fx.gk.forward_premium`", "calls and puts; any broadcast shape",
     "`test_forward_premium_*`"),
    ("`delta`", "`qef.fx.gk.delta`", "spot and forward, pips and premium-adjusted; the nine G10 conventions",
     "`test_delta_*`"),
    ("`sabr_vol`", "`qef.fx.sabr.sabr_vol`, β = 1", "both z/x(z) branches and the branch point, ν = 0, "
     "ρ near ±1, deep out-of-the-money strikes, NaN and infinite inputs", "`test_sabr_vol_*`"),
    ("`pa_call_delta_maximiser`", "`qef.fx.gk.pa_call_delta_maximiser`", "bracket expansion and root",
     "`test_strike_near_maximal_*`"),
    ("`strike_from_delta`", "`qef.fx.gk.strike_from_delta`", "pips closed form; premium-adjusted put root; "
     "premium-adjusted call root on [K*, ∞); every check and error of the reference",
     "`test_strike_from_delta_*`, `test_strike_near_maximal_*`"),
    ("`middle_contracts`, `middle_contracts_fixed`, `middle_contracts_sabr`",
     "`qef.fx.moments.middle_contracts` for a SABR smile", "k = 1, 2, 3; EURUSD-type and USD-base pairs; "
     "fixed n and the adaptive doubling; many smiles in one call", "`test_middle_contracts_*`"),
]


def timed(fn, reps):
    fn()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return statistics.median(ts)


def rel(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.abs(a - b) / np.abs(b)
    return np.where(a == b, 0.0, r)


def build():
    out = subprocess.run(["bash", str(ROOT / "cpp" / "build.sh")], check=True, capture_output=True, text=True).stdout
    lines = dict(line.split(": ", 1) for line in out.splitlines() if ": " in line)
    return lines["compiler"], lines["command"]


def option_grid(rng, n):
    tau = rng.choice([1 / 12, 0.25], n)
    sigma = rng.uniform(0.03, 0.4, n)
    F = np.array(list(LEVELS.values()))[rng.integers(0, len(LEVELS), n)] * np.exp(rng.uniform(-0.3, 0.3, n))
    s = sigma * np.sqrt(tau)
    K = F * np.exp(0.5 * s**2 - s * special.ndtri(rng.uniform(0.05, 0.95, n)))
    return F, K, sigma, tau, rng.choice([-1.0, 1.0], n), rng.uniform(0.97, 1.01, n)


def smile_panel(rng):
    rows = []
    for ccy in LEVELS:
        for _ in range(MONTHS):
            F = LEVELS[ccy] * math.exp(rng.uniform(-0.2, 0.2))
            smile = SabrSmile(F, 1 / 12, rng.uniform(0.05, 0.20), rng.uniform(-0.5, 0.5), rng.uniform(0.2, 2.0))
            conv = G10[ccy]
            K_min = gk.strike_from_delta_smile(-0.10, F, smile.tau, gk.PUT, conv.delta, smile.vol)
            K_max = gk.strike_from_delta_smile(0.10, F, smile.tau, gk.CALL, conv.delta, smile.vol)
            rows.append((ccy, smile, K_min, K_max, conv.usd_base))
    return rows


def reference_integrate(F, K_min, K_max, tau, vol_of_strike, usd_base, n):
    """The inner integrate(n) of qef.fx.moments.middle_contracts, line for line."""
    total = np.zeros(len(moments.CONTRACTS))
    for a, b in ((K_min, F), (F, K_max)):
        K = np.linspace(a, b, n + 1)
        p = moments.otm_price(K, F, tau, vol_of_strike)
        total += [moments._simpson(moments.contract_weight(k, K, F, usd_base) * p, (b - a) / n)
                  for k in moments.CONTRACTS]
    return total


def abs_integral(F, K_min, K_max, tau, vol_of_strike, usd_base, n):
    """Simpson's rule for |h_k'' P| on the reference's nodes (the normwise scale of the tests)."""
    total = np.zeros(len(moments.CONTRACTS))
    for a, b in ((K_min, F), (F, K_max)):
        K = np.linspace(a, b, n + 1)
        p = moments.otm_price(K, F, tau, vol_of_strike)
        total += [moments._simpson(np.abs(moments.contract_weight(k, K, F, usd_base) * p), (b - a) / n)
                  for k in moments.CONTRACTS]
    return total


def exact_check(F, K, phi, d1, d2, ref, out, r, top=20):
    """Errors of both premia against 50 digits at the ``top`` points of largest relative disagreement.

    Both implementations form d± identically in double precision, so the
    50-digit value φ[F Φ(φd+) − K Φ(φd−)] at those d± isolates the rounding
    of each implementation's Φ, products and difference.
    """
    try:
        import mpmath
    except ImportError:
        return "mpmath not installed, 50-digit check skipped"
    mpmath.mp.dps = 50
    ek, er = [], []
    for i in np.argsort(-r)[:top]:
        ph = mpmath.mpf(float(phi[i]))
        exact = ph * (mpmath.mpf(float(F[i])) * mpmath.ncdf(ph * mpmath.mpf(float(d1[i])))
                      - mpmath.mpf(float(K[i])) * mpmath.ncdf(ph * mpmath.mpf(float(d2[i]))))
        ek.append(abs(float((mpmath.mpf(float(out[i])) - exact) / exact)))
        er.append(abs(float((mpmath.mpf(float(ref[i])) - exact) / exact)))
    ek, er = np.array(ek), np.array(er)
    return (f"at the {top} largest disagreements the kernel is closer to a 50-digit value at {np.sum(ek < er)}, "
            f"with largest relative errors {fmt_e(ek.max())} (kernel) and {fmt_e(er.max())} (reference)")


def fmt_t(t):
    return f"{t * 1e3:,.1f} ms" if t >= 1e-3 else f"{t * 1e6:,.0f} µs"


def fmt_e(x):
    return "0" if x == 0 else f"{x:.1e}"


def main():
    compiler, command = build()
    if not kk.available():
        raise SystemExit("the kernel did not build")
    rng = np.random.default_rng(SEED)
    timings, parity = [], []

    # Premium and delta, 1,000,000 options.
    F, K, sigma, tau, phi, df = option_grid(rng, 1_000_000)
    ref = gk.forward_premium(F, K, sigma, tau, phi)
    out = kk.forward_premium(F, K, sigma, tau, phi)
    d1, d2 = gk.d_plus_minus(F, K, sigma, tau)
    scale = F * special.ndtr(phi * d1) + K * special.ndtr(phi * d2)
    r = rel(out, ref)
    timings.append(("`forward_premium`", "1,000,000 options (both vectorised)",
                    timed(lambda: gk.forward_premium(F, K, sigma, tau, phi), REPS),
                    timed(lambda: kk.forward_premium(F, K, sigma, tau, phi), REPS)))
    parity.append(("`forward_premium`", "1,000,000", fmt_e(r.max()),
                   f"normwise {fmt_e(np.max(np.abs(out - ref) / scale))}; relative ≤ 1e-13 for "
                   f"{np.mean(r <= 1e-13):.2%}, median {fmt_e(np.median(r))}; {exact_check(F, K, phi, d1, d2, ref, out, r)}",
                   f"{np.mean(out == ref):.1%}"))
    for conv in CONVENTIONS:
        name = f"`delta` ({'spot' if conv.spot else 'forward'}, {'premium-adjusted' if conv.premium_adjusted else 'pips'})"
        ref, out = gk.delta(F, K, sigma, tau, phi, conv, df), kk.delta(F, K, sigma, tau, phi, conv, df)
        timings.append((name, "1,000,000 options (both vectorised)",
                        timed(lambda: gk.delta(F, K, sigma, tau, phi, conv, df), REPS),
                        timed(lambda: kk.delta(F, K, sigma, tau, phi, conv, df), REPS)))
        parity.append((name, "1,000,000", fmt_e(rel(out, ref).max()), "", f"{np.mean(out == ref):.1%}"))

    # SABR.
    Ks = 1.1 * np.exp(rng.uniform(-0.25, 0.25, 1_000_000))
    one = (Ks, 1.1, 1 / 12, 0.1, -0.3, 1.2)
    ref, out = sabr.sabr_vol(*one), kk.sabr_vol(*one)
    timings.append(("`sabr_vol`", "one smile, 1,000,000 strikes (both vectorised)",
                    timed(lambda: sabr.sabr_vol(*one), REPS), timed(lambda: kk.sabr_vol(*one), REPS)))
    sabr_rel, sabr_eq = [rel(out, ref).max()], [np.mean(out == ref)]
    panel = smile_panel(rng)
    Kg = np.array([np.linspace(K_min, K_max, 513) for _, _, K_min, K_max, _ in panel])
    pars = [np.array([getattr(s, a) for _, s, *_ in panel])[:, None] for a in ("F", "tau", "alpha", "rho", "nu")]
    loop = lambda: [sabr.sabr_vol(Kg[i], s.F, s.tau, s.alpha, s.rho, s.nu) for i, (_, s, *_) in enumerate(panel)]
    ref, out = np.array(loop()), kk.sabr_vol(Kg, *pars)
    timings.append(("`sabr_vol`", "1,440 smiles × 513 strikes (reference loops over smiles)",
                    timed(loop, REPS), timed(lambda: kk.sabr_vol(Kg, *pars), REPS)))
    sabr_rel.append(rel(out, ref).max())
    sabr_eq.append(np.mean(out == ref))
    parity.append(("`sabr_vol`", f"{Ks.size + Kg.size:,}", fmt_e(max(sabr_rel)), "",
                   f"{min(sabr_eq):.1%}"))

    # Strikes from delta: the 10-delta put and call of every smile, per convention type.
    for pa in (False, True):
        rows = [(s, G10[c].delta) for c, s, *_ in panel if G10[c].delta.premium_adjusted == pa]
        conv = rows[0][1]
        t = np.tile([-0.10, 0.10], len(rows))
        ph = np.sign(t)
        Fs = np.repeat([s.F for s, _ in rows], 2)
        sig = np.repeat([s.alpha for s, _ in rows], 2)
        loop = lambda: [gk.strike_from_delta(t[i], Fs[i], sig[i], 1 / 12, ph[i], conv) for i in range(t.size)]
        ref, out = np.array(loop()), kk.strike_from_delta(t, Fs, sig, 1 / 12, ph, conv)
        name = f"`strike_from_delta` ({'premium-adjusted' if pa else 'pips'})"
        timings.append((name, f"{t.size:,} strikes (reference loops)", timed(loop, REPS_LOOP),
                        timed(lambda: kk.strike_from_delta(t, Fs, sig, 1 / 12, ph, conv), REPS)))
        extra = ""
        if pa:
            bound = kk.PARITY_TOLERANCES["root_xtol_per_unit_forward"] * Fs + (
                kk.PARITY_TOLERANCES["root_rtol"] + kk.PARITY_TOLERANCES["root_slack"]) * np.abs(ref)
            extra = f"largest ΔK as a share of the root bound {np.max(np.abs(out - ref) / bound):.2f}"
        parity.append((name, f"{t.size:,}", fmt_e(rel(out, ref).max()), extra, f"{np.mean(out == ref):.1%}"))

    # Middle part, 1,440 smiles.
    args = [(s.F, K_min, K_max, s.tau, s, usd) for _, s, K_min, K_max, usd in panel]
    batch = [np.array(col) for col in zip(*[(s.F, K_min, K_max, s.tau, s.alpha, s.rho, s.nu, usd)
                                            for _, s, K_min, K_max, usd in panel])]
    ref_loop = lambda: [moments.middle_contracts(F_, lo, hi, t_, s.vol, u) for F_, lo, hi, t_, s, u in args]
    refs = ref_loop()
    out = kk.middle_contracts_sabr(*batch)
    n_ref = np.array([r_["n"] for r_ in refs])
    n_fixed = int(np.median(n_ref))
    t_ref = timed(ref_loop, REPS_LOOP)
    timings.append(("middle part, adaptive", "1,440 SABR smiles, one kernel call", t_ref,
                    timed(lambda: kk.middle_contracts_sabr(*batch), REPS)))
    timings.append(("middle part, adaptive", "1,440 SABR smiles, one kernel call per smile", t_ref,
                    timed(lambda: [kk.middle_contracts(F_, lo, hi, t_, s.vol, u) for F_, lo, hi, t_, s, u in args],
                          REPS)))
    vals_ref = np.array([r_["values"] for r_ in refs])
    A = np.array([abs_integral(F_, lo, hi, t_, s.vol, u, n) for (F_, lo, hi, t_, s, u), n in zip(args, n_ref)])
    same_n = np.mean(out["n"] == n_ref)
    same_conv = np.mean(out["converged"] == np.array([r_["converged"] for r_ in refs]))
    parity.append(("middle part, adaptive", "1,440 smiles × 3 contracts",
                   f"k = 1, 2: {fmt_e(rel(out['values'][:, :2], vals_ref[:, :2]).max())}; "
                   f"k = 3: {fmt_e(rel(out['values'][:, 2], vals_ref[:, 2]).max())}",
                   f"normwise {fmt_e(np.max(np.abs(out['values'] - vals_ref) / A))}; same n in {same_n:.0%} and "
                   f"same convergence flag in {same_conv:.0%} of smiles", f"{np.mean(out['values'] == vals_ref):.1%}"))
    for a in args[:20]:  # the reproduced integrate(n) is the reference's
        assert np.array_equal(reference_integrate(*a[:4], a[4].vol, a[5], n_fixed),
                              moments.middle_contracts(*a[:4], a[4].vol, a[5], n0=n_fixed // 2, n_max=n_fixed)["values"])
    fixed_loop = lambda: [reference_integrate(F_, lo, hi, t_, s.vol, u, n_fixed) for F_, lo, hi, t_, s, u in args]
    vals_ref = np.array(fixed_loop())
    out = kk.middle_contracts_sabr(*batch, n=n_fixed)
    A = np.array([abs_integral(F_, lo, hi, t_, s.vol, u, n_fixed) for F_, lo, hi, t_, s, u in args])
    timings.append(("middle part, fixed n", f"1,440 SABR smiles at n = {n_fixed} (reference integrate(n) loop)",
                    timed(fixed_loop, REPS_LOOP), timed(lambda: kk.middle_contracts_sabr(*batch, n=n_fixed), REPS)))
    parity.append((f"middle part, n = {n_fixed}", "1,440 smiles × 3 contracts",
                   f"k = 1, 2: {fmt_e(rel(out['values'][:, :2], vals_ref[:, :2]).max())}; "
                   f"k = 3: {fmt_e(rel(out['values'][:, 2], vals_ref[:, 2]).max())}",
                   f"normwise {fmt_e(np.max(np.abs(out['values'] - vals_ref) / A))}",
                   f"{np.mean(out['values'] == vals_ref):.1%}"))

    n_q = np.quantile(n_ref, [0.0, 0.5, 1.0]).astype(int)
    write_readme(compiler, command, timings, parity, n_q)
    for name, work, t_ref, t_k in timings:
        print(f"{name:55s} {work:62s} {fmt_t(t_ref):>12s} {fmt_t(t_k):>12s} {t_ref / t_k:9.1f}x")


def write_readme(compiler, command, timings, parity, n_q):
    cpu = platform.processor() or platform.machine()
    try:
        cpu = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True,
                             check=True).stdout.strip() or cpu
    except (OSError, subprocess.CalledProcessError):
        pass
    tol = kk.PARITY_TOLERANCES
    lines = [
        BEGIN,
        "",
        f"Generated by `scripts/benchmark_kernel.py` on {datetime.date.today():%d %B %Y}; do not edit by hand.",
        "",
        "### Build",
        "",
        f"Compiler: {compiler}. Command, as run by `cpp/build.sh`:",
        "",
        "```",
        command,
        "```",
        "",
        "### Coverage",
        "",
        "| Kernel wrapper (`qef.kernel`) | Python reference | Branches and conventions covered | Parity tests |",
        "|---|---|---|---|",
        *[f"| {a} | {b} | {c} | {d} |" for a, b, c, d in COVERAGE],
        "",
        "### Parity tolerances",
        "",
        "From `qef.kernel.PARITY_TOLERANCES`; `tests/test_cpp_kernel.py` states why each looser bound holds.",
        "",
        "| Quantity | Criterion |",
        "|---|---|",
        f"| normal distribution function, deltas, SABR volatilities, pips strikes | relative error ≤ {tol['relative']:g} |",
        f"| inverse normal | absolute error ≤ {tol['relative']:g} max(1, \\|x\\|) |",
        f"| premium | \\|ΔC\\| ≤ {tol['premium_normwise']:g} [F Φ(φd+) + K Φ(φd−)]; relative ≤ {tol['relative']:g} "
        f"where that sum over \\|C\\| is at most {tol['premium_condition_for_relative']:g} |",
        f"| premium-adjusted strikes | \\|ΔK\\| ≤ {tol['root_xtol_per_unit_forward']:g} F + {tol['root_rtol']:g} "
        f"\\|K\\| + {tol['root_slack']:g} \\|K\\| (the reference's brentq tolerance plus 4.5 ulp); near the maximal "
        "call delta, the reference's delta at the kernel's strike within 1e-14 of the target |",
        f"| middle part | \\|ΔI_k\\| ≤ {tol['middle_normwise']:g} A_k, A_k Simpson's rule for \\|h_k'' P\\|; relative "
        f"≤ {tol['relative']:g} where A_k/\\|I_k\\| ≤ {tol['premium_condition_for_relative']:g} (always for k = 1, 2) |",
        "",
        "### Timings",
        "",
        f"Machine: {cpu}, {platform.system()} {platform.release()}; Python {platform.python_version()}, "
        f"NumPy {np.__version__}, SciPy {scipy.__version__}; one thread. Median of {REPS} repetitions "
        f"({REPS_LOOP} for the reference's Python loops) after a warm-up call. In the adaptive middle part the "
        f"reference converged at n from {n_q[0]} to {n_q[2]} subintervals per side (median {n_q[1]}).",
        "",
        "| Kernel | Workload | Python reference | C++ kernel | Speed-up |",
        "|---|---|---:|---:|---:|",
        *[f"| {a} | {w} | {fmt_t(r)} | {fmt_t(k)} | {r / k:,.1f}× |" for a, w, r, k in timings],
        "",
        "### Parity on the benchmark workloads",
        "",
        "| Kernel | Results | Largest relative error | Other criteria | Bitwise equal |",
        "|---|---:|---|---|---:|",
        *[f"| {a} | {b} | {c} | {d} | {e} |" for a, b, c, d, e in parity],
        "",
        END,
    ]
    text = README.read_text() if README.exists() else ""
    block = "\n".join(lines)
    if BEGIN in text and END in text:
        text = text[: text.index(BEGIN)] + block + text[text.index(END) + len(END):]
    else:
        text = text.rstrip() + "\n\n## Benchmark\n\n" + block + "\n"
    README.write_text(text)


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        main()
