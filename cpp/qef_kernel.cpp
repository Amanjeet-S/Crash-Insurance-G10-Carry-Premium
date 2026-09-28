/*
 * C++ kernel for the pricing, delta, SABR and middle-part moment code of
 * project 01 (research design, Stage 5: "the C++ kernel with parity against
 * the Python reference"). The C interface is declared in qef_kernel.h and
 * wrapped for Python by src/qef/kernel.py through ctypes.
 *
 * Each routine mirrors one Python function operation by operation, in the
 * same order and with the same floating-point expressions, so that the two
 * differ only where a library function differs:
 *
 *   qef_forward_premium  qef.fx.gk.forward_premium
 *       phi [F Phi(phi d+) - K Phi(phi d-)], d+- = [ln(F/K) +- s^2/2] / s,
 *       s = sigma sqrt(tau);
 *   qef_delta            qef.fx.gk.delta
 *       pips: phi D Phi(phi d+); premium-adjusted: phi D (K/F) Phi(phi d-),
 *       with D = D_b for spot deltas and D = 1 for forward deltas;
 *   qef_sabr_vol         qef.fx.sabr.sabr_vol with beta = 1
 *       alpha (z/x(z)) [1 + (rho nu alpha/4 + (2 - 3 rho^2) nu^2/24) tau],
 *       z = (nu/alpha) ln(F/K), x(z) = ln[(sqrt(1 - 2 rho z + z^2) + z - rho)
 *       / (1 - rho)], z/x(z) = 1 - rho z/2 + (2 - 3 rho^2) z^2/12 for
 *       |z| < 1e-7, and z/x(z) = 1 when nu is not positive;
 *   qef_pa_call_delta_maximiser, qef_strike_from_delta
 *                        qef.fx.gk.pa_call_delta_maximiser, strike_from_delta
 *       pips: K = F exp(s^2/2 - phi s Phi^{-1}(|delta|/D)) (result R5(a));
 *       premium-adjusted puts: the root on the reference's bracket; calls:
 *       the root on [K*, inf), K* the maximiser of the call delta (R5(b));
 *   qef_middle_contracts_fixed, qef_middle_contracts_adaptive
 *                        qef.fx.moments.middle_contracts for a SABR smile
 *       composite Simpson's rule for the integrals of h_k'' P over [K_min, F]
 *       and of h_k'' C over [F, K_max], k = 1, 2, 3, with n subintervals per
 *       side (fixed), or the reference's doubling of n from n0 (adaptive).
 *
 * Sources. The formulas are those of the Python reference modules
 * src/qef/fx/gk.py, sabr.py and moments.py as of 27 September 2026, which I
 * read line by line for this port. Their citations are inherited and were
 * not consulted for the kernel: Garman and Kohlhagen (1983) and Reiswich and
 * Wystup (2012) for pricing and delta conventions, Hagan, Kumar, Lesniewski
 * and Woodward (2002) for the SABR expansion; full references are in
 * projects/01_fx_crash_insurance_carry/references.md. Library behaviour that
 * the port reproduces was read in the installed sources: numpy 2.3.5
 * (numpy.linspace; the numpy.sum docstring on partial pairwise summation) and
 * scipy 1.16.3 (scipy.stats.norm evaluates Phi by scipy.special.ndtr, its
 * inverse by scipy.special.ndtri and phi as np.exp(-x**2/2.0) divided by
 * np.sqrt(2*np.pi); the scipy.optimize.brentq docstring states
 * that its root x0 satisfies |x - x0| <= xtol + rtol |x0| for the exact root
 * x). The blocking of numpy's pairwise sum (below) is not documented; I
 * verified it bitwise against numpy.sum on strided arrays of 2 to 40,000
 * elements and did not consult the numpy C source.
 *
 * Normal distribution. Phi(x) = erfc(-x/sqrt(2))/2 with the C library erfc,
 * accurate in both tails. scipy's ndtr forms the same argument but uses its
 * own erfc (for x < -1 it equals erfc(-x/sqrt(2))/2 with scipy.special.erfc
 * bitwise on the grid below). Measured on this machine, the two functions
 * differ by about 20 units in the last place for |x| < 6: at most 19 on
 * 2,000,001 evenly spaced points of [-6, 6] and 21 on 5 million uniform
 * draws, largest for x < -4 and at most 2 for x >= 0. These are sample
 * maxima, not a bound, and the difference is not the error of either
 * function alone. In the lower tail the condition number of Phi,
 * |x| phi(x)/Phi(x), grows like x^2, so rounding the argument alone costs up
 * to tens of units in the last place, and against 40-digit values both
 * functions are accurate there only to tens of units in the last place (at
 * most 54 each on 400 random points of [-6, -4]). That difference, amplified
 * by the cancellation in an out-of-the-money premium, is the one source of
 * disagreement larger than a few ulps (see tests/test_cpp_kernel.py).
 * The density is formed as scipy.stats.norm.pdf forms it,
 * phi(x) = exp(-x^2/2) / np.sqrt(2 pi), whose divisor is one unit in the
 * last place below the double nearest sqrt(2 pi) that the tail expansion
 * below uses; I checked it bitwise against norm.pdf on 5 million points.
 * (Corrected on 27 September 2026 after an independent review: the kernel
 * multiplied by 1/sqrt(2 pi), which differs from norm.pdf by one unit in the
 * last place at about a quarter of arguments, and it quoted the difference
 * of Phi and ndtr as at most 13 units in the last place.)
 *
 * Phi^{-1} is computed here rather than taken from a rational approximation:
 * for 0 < p <= 1/2, Newton's method on h(x) = ln Phi(x) - ln p from
 * x0 = -sqrt(-2 ln p), then two Newton steps on Phi(x) - p; for p > 1/2,
 * Phi^{-1}(p) = -Phi^{-1}(1 - p), where 1 - p is exact. Convergence without
 * overshoot: h is increasing and concave, since
 * (ln Phi)'' = -phi(x) [x Phi(x) + phi(x)] / Phi(x)^2 and
 * x Phi(x) + phi(x) = int_{-inf}^x Phi(t) dt > 0, so Newton's iterates rise
 * monotonically to the root from any point on its left; and x0 is on its left
 * because Phi(-u) = P(Z >= u) <= E exp(u(Z - u)) = exp(-u^2/2) = p for
 * u = sqrt(-2 ln p) (Markov's inequality applied to exp(uZ)).
 * For subnormal p (below DBL_MIN, so x < -37.5) Phi(x) itself is subnormal
 * and ln Phi(x) loses digits, so Newton's method runs on ln Phi from the
 * asymptotic tail expansion Phi(-u) = phi(u)/u [1 - u^-2 + 3u^-4 - 15u^-6 +
 * 105u^-8 - 945u^-10 + R], 0 <= R <= 10395 u^-12, which at u > 37.5 is
 * below 2e-15 relative; phi/Phi = u / [1 - u^-2 + ...] is the derivative. The
 * bound on R is mine: Phi(-u)/phi(u) = int_0^inf exp(-ut - t^2/2) dt, and
 * exp(-y) for y >= 0 lies between consecutive partial sums of its Taylor
 * series, so integrating term by term brackets the ratio between
 * consecutive partial sums of the expansion.
 *
 * Root finding. The premium-adjusted strike is a root of delta(K) - target on
 * the reference's own bracket, with the reference's checks in the order of
 * scipy.optimize.brentq: the argument check, NaN function values (raised by
 * its Python wrapper, which I read), an end value of zero returned as the
 * root, then no sign change on the bracket, tested by comparing sign bits.
 * The order of the last two I established by calling brentq: an end value
 * of -0.0 is returned as the root, and two tiny end values of one sign
 * raise. (Corrected on 27 September 2026 after an independent review: the
 * kernel tested the product of the end values, which underflows to zero when
 * both are tiny, and then returned a point that was not a root.) The kernel
 * then runs Newton's method safeguarded by bisection on that bracket, with
 * the derivative d delta/dK = phi D [Phi(phi d-) - phi phi(d-)/s] / F (my
 * derivation: d/dK Phi(phi d-) = -phi phi(d-) / (K s)), until the step or the
 * bracket is within a few units in the last place. The kernel's root is thus
 * accurate to a few units in the last place, up to the conditioning of the
 * delta equation, whereas the reference's is guaranteed only to
 * xtol = 1e-14 F plus rtol = 1e-14 relative. On the first 150 retained
 * premium-adjusted targets of each premium-adjusted convention of the parity
 * grid of test_strike_from_delta_parity (300 in all), against 40-digit roots
 * in units of np.spacing of the root, the kernel's errors are at most 4.6
 * (median 0.6) and the reference's at most 36 (median 0.5), although the
 * kernel's root is the further from the exact one at 108 of the 300 points;
 * these are sample maxima, not bounds. The maximiser K*
 * solves g(d) = s Phi(d) - phi(d) = 0 with g'(d) = phi(d)(s + d), on the
 * reference's bracket, and K* = F exp(-d* s - s^2/2).
 * Both this root and the reference's lose accuracy as s grows: near the
 * root s + d* is about 1/s, so g'(d*) is about phi(d*)/s, while rounding the
 * arguments of erfc and exp perturbs g by about s^2 units in the last place
 * of phi(d*); the error in ln K* is s times that in d*, of order s^4 units in
 * the last place. Against 50-digit values at 300 geometrically spaced values
 * of s from 0.01 to 20, both relative errors are below 1e-14 for s < 2 and
 * reach 1.3e-11 (kernel) and 1.25e-11 (reference) near s = 20; in this
 * project s is below 0.3. (Remeasured on 27 September 2026 after the change
 * to phi above; the kernel's figure near s = 20 was 1.7e-11 before it.)
 *
 * Summation. Simpson's rule sums the odd and even interior nodes with numpy's
 * partial pairwise summation (fewer than 8 terms: a running sum; up to 128:
 * eight interleaved accumulators combined as ((r0+r1)+(r2+r3))+((r4+r5)+
 * (r6+r7)), then the remainder; above 128: split at the largest multiple of 8
 * not above half and recurse), and nodes are placed as numpy.linspace places
 * them (i * step + a, the last node exactly b). Between successive levels of
 * the adaptive rule the nodes of level n are the even nodes of level 2n
 * bitwise. Level 2n forms its step as delta/(2n), not as the step of level n
 * halved, but when the result is a normal number the two are the same double,
 * because scaling by 2 commutes with rounding; 2i times that step is then the
 * same real number as i times delta/n, so the two products round alike and
 * node 2i of level 2n equals node i of level n. The kernel therefore reuses
 * their integrand values instead of recomputing them (fill_side checks that
 * the step is normal).
 *
 * Floating point. The file must be compiled without -ffast-math and without
 * contraction of a * b + c into fused multiply-adds (build.sh passes
 * -ffp-contract=off; the pragma below does the same for clang), because the
 * reference rounds every product and sum separately.
 */

#include "qef_kernel.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
#include <new>
#include <vector>

#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#endif

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kSqrtHalf = 0.70710678118654752440;    // 1/sqrt(2), as M_SQRT1_2
constexpr double kSqrt2Pi = 2.50662827463100050242;     // sqrt(2 pi), the nearest double (tail expansion)
constexpr double kNumpySqrt2Pi = 0x1.40d931ff62705p+1;  // np.sqrt(2 * np.pi), scipy's divisor in norm.pdf
constexpr int kMaxExpansions = 200;                     // qef.fx.gk.MAX_EXPANSIONS
constexpr int kMaxIterations = 500;                     // bisection alone needs about 260 on the widest bracket
constexpr int64_t kMaxSubintervals = int64_t(1) << 22;  // bounds n, n0 and n_max, so a level has at most 2^23 subintervals

// Strided input: increment 0 broadcasts a scalar, increment 1 walks an array.
struct In {
    const double *p;
    int64_t s;
    double operator[](int64_t i) const { return p[i * s]; }
};

inline double norm_cdf(double x) { return 0.5 * std::erfc(-(x * kSqrtHalf)); }

// scipy.stats.norm.pdf: np.exp(-x**2/2.0) / np.sqrt(2*np.pi).
inline double norm_pdf(double x) { return std::exp(-(x * x) / 2.0) / kNumpySqrt2Pi; }

// ln Phi(x) and phi(x)/Phi(x) where Phi(x) is below the smallest normal
// number (x < -37.5), from the tail expansion (header comment).
inline void log_norm_cdf_tail(double x, double &log_cdf, double &pdf_over_cdf) {
    const double u = -x, v = 1.0 / (u * u);
    const double series = 1.0 + v * (-1.0 + v * (3.0 + v * (-15.0 + v * (105.0 - v * 945.0))));
    log_cdf = -0.5 * (x * x) - std::log(u * kSqrt2Pi) + std::log(series);
    pdf_over_cdf = u / series;
}

// Phi^{-1}(p) for 0 < p <= 1/2 (header comment).
double norm_ppf_lower(double p) {
    if (p == 0.5) return 0.0;
    const double lp = std::log(p);
    double x = -std::sqrt(-2.0 * lp);
    if (p < DBL_MIN) {  // subnormal p: Phi(x) is subnormal too, so ln Phi comes from the tail expansion
        for (int it = 0; it < 100; ++it) {
            double lc, r;
            log_norm_cdf_tail(x, lc, r);
            const double step = (lc - lp) / r;
            x -= step;
            if (std::fabs(step) <= 1e-15 * std::fabs(x)) break;
        }
        return x;
    }
    for (int it = 0; it < 100; ++it) {
        const double c = norm_cdf(x), d = norm_pdf(x);
        if (!(c > 0.0 && d > 0.0)) break;  // only below the smallest normal p
        const double step = (std::log(c) - lp) * c / d;
        x -= step;
        if (std::fabs(step) <= 1e-12 * (1.0 + std::fabs(x))) break;
    }
    for (int k = 0; k < 2; ++k) {
        const double d = norm_pdf(x);
        if (!(d > 0.0)) break;
        x -= (norm_cdf(x) - p) / d;
    }
    return x;
}

double norm_ppf(double p) {
    if (std::isnan(p) || p < 0.0 || p > 1.0) return kNaN;
    if (p == 0.0) return -kInf;
    if (p == 1.0) return kInf;
    if (p > 0.5) return -norm_ppf_lower(1.0 - p);
    return norm_ppf_lower(p);
}

// gk.d_plus_minus with log_fk = ln(F/K) and s = sigma sqrt(tau) already formed.
inline void d_plus_minus_log(double log_fk, double s, double &d1, double &d2) {
    d1 = (log_fk + 0.5 * (s * s)) / s;
    d2 = d1 - s;
}

inline void d_plus_minus(double F, double K, double s, double &d1, double &d2) {
    d_plus_minus_log(std::log(F / K), s, d1, d2);
}

inline double forward_premium_log(double F, double K, double log_fk, double sigma, double sqrt_tau, double phi) {
    double d1, d2;
    d_plus_minus_log(log_fk, sigma * sqrt_tau, d1, d2);
    return phi * (F * norm_cdf(phi * d1) - K * norm_cdf(phi * d2));
}

inline double forward_premium_1(double F, double K, double sigma, double sqrt_tau, double phi) {
    return forward_premium_log(F, K, std::log(F / K), sigma, sqrt_tau, phi);
}

inline double delta_1(double F, double K, double sigma, double sqrt_tau, double phi, double scale, bool pa) {
    double d1, d2;
    d_plus_minus(F, K, sigma * sqrt_tau, d1, d2);
    if (pa) return phi * scale * (K / F) * norm_cdf(phi * d2);
    return phi * scale * norm_cdf(phi * d1);
}

// sabr._z_over_x.
inline double z_over_x(double z, double rho) {
    if (std::fabs(z) < 1e-7) return 1.0 - 0.5 * rho * z + (2.0 - 3.0 * (rho * rho)) * (z * z) / 12.0;
    const double x = std::log((std::sqrt(1.0 - 2.0 * rho * z + z * z) + z - rho) / (1.0 - rho));
    return z / x;
}

// sabr.sabr_vol at beta = 1, where (F K)^{(1 - beta)/2} = 1. The terms with
// factor (1 - beta)^2 = 0 are kept because the reference forms them: they are
// zero for finite arguments and NaN for infinite ones, as in the reference.
inline double sabr_vol_log(double log_fk, double tau, double alpha, double rho, double nu) {
    const double l2 = log_fk * log_fk;
    const double denom = 1.0 + 0.0 * l2 + 0.0 * (l2 * l2);
    const double ratio = nu > 0 ? z_over_x((nu / alpha) * log_fk, rho) : 1.0;
    const double correction =
        1.0 + (0.0 * (alpha * alpha) + 0.25 * rho * nu * alpha + (2.0 - 3.0 * (rho * rho)) / 24.0 * (nu * nu)) * tau;
    return alpha / denom * ratio * correction;
}

inline double sabr_vol_1(double K, double F, double tau, double alpha, double rho, double nu) {
    return sabr_vol_log(std::log(F / K), tau, alpha, rho, nu);
}

inline bool vol_ok(double sigma, double tau) {  // gk._check_vol
    return std::isfinite(sigma) && std::isfinite(tau) && sigma > 0 && tau > 0;
}

// Root of f on [a, b] (see the header comment). fd(x, f, df) evaluates f and
// its derivative; x0 is a starting point, replaced by the midpoint if it is
// not inside the bracket.
template <class FD>
int32_t solve_bracketed(const FD &fd, double a, double b, double xtol, double x0, double &root) {
    if (!(xtol > 0.0)) return QEF_ERR_XTOL;
    double fa, fb, unused;
    fd(a, fa, unused);
    fd(b, fb, unused);
    if (std::isnan(fa) || std::isnan(fb)) return QEF_ERR_NAN;
    // brentq's order: an end value of zero is the root, then the signs are
    // compared. signbit cannot underflow as the product fa * fb can; the zero
    // tests come first because signbit(-0.0) is true.
    if (fa == 0.0) {
        root = a;
        return QEF_OK;
    }
    if (fb == 0.0) {
        root = b;
        return QEF_OK;
    }
    if (std::signbit(fa) == std::signbit(fb)) return QEF_ERR_SAME_SIGN;
    const bool fa_negative = fa < 0.0;
    double x = x0;
    if (!(x > std::min(a, b) && x < std::max(a, b))) x = 0.5 * (a + b);
    for (int it = 0; it < kMaxIterations; ++it) {
        double fx, dfx;
        fd(x, fx, dfx);
        if (std::isnan(fx)) return QEF_ERR_NAN;
        if (fx == 0.0) {
            root = x;
            return QEF_OK;
        }
        if ((fx < 0.0) == fa_negative) a = x; else b = x;
        const double lo = std::min(a, b), hi = std::max(a, b);
        double xn = x - fx / dfx;
        if (!(xn > lo && xn < hi)) xn = 0.5 * (lo + hi);
        const double step = std::fabs(xn - x);
        x = xn;
        if (step <= 2.0 * DBL_EPSILON * std::fabs(x) || hi - lo <= 4.0 * DBL_EPSILON * std::fabs(hi)) {
            root = x;
            return QEF_OK;
        }
    }
    root = x;
    return QEF_ERR_NO_CONVERGENCE;
}

int32_t pa_call_delta_maximiser_1(double F, double sigma, double tau, double &k_star) {
    if (!vol_ok(sigma, tau)) return QEF_ERR_VOL;
    const double s = sigma * std::sqrt(tau);
    auto g = [s](double d) { return s * norm_cdf(d) - norm_pdf(d); };
    double lo = -1.0, hi = 1.0;
    for (int i = 0; i < kMaxExpansions; ++i) {
        if (g(lo) <= 0) break;
        lo *= 2.0;
    }
    for (int i = 0; i < kMaxExpansions; ++i) {
        if (g(hi) >= 0) break;
        hi *= 2.0;
    }
    auto fd = [s](double d, double &f, double &df) {
        const double pdf = norm_pdf(d);
        f = s * norm_cdf(d) - pdf;
        df = pdf * (s + d);
    };
    double d_star = 0.0;
    const int32_t st = solve_bracketed(fd, lo, hi, 1e-14, 0.5 * (lo + hi), d_star);
    if (st != QEF_OK) return st;
    k_star = F * std::exp(-d_star * s - 0.5 * (s * s));
    return QEF_OK;
}

// gk.strike_from_delta for one element; the checks come in the reference's order.
int32_t strike_from_delta_1(double target, double F, double sigma, double tau, double phi, double df_base, bool spot,
                            bool pa, double &K) {
    K = kNaN;
    const double sign = target > 0 ? 1.0 : (target < 0 ? -1.0 : (target == 0 ? 0.0 : kNaN));  // np.sign
    if (!(sign == phi) || target == 0) return QEF_ERR_SIGN;
    if (!vol_ok(sigma, tau)) return QEF_ERR_VOL;
    const double scale = spot ? df_base : 1.0;
    const double sqrt_tau = std::sqrt(tau);
    const double s = sigma * sqrt_tau;
    if (!pa) {
        const double a = std::fabs(target) / scale;
        if (!(0 < a && a < 1)) return QEF_ERR_PIPS;
        K = F * std::exp(0.5 * (s * s) - phi * s * norm_ppf(a));
        return QEF_OK;
    }
    auto f = [&](double k) { return delta_1(F, k, sigma, sqrt_tau, phi, scale, true) - target; };
    auto fd = [&](double k, double &fk, double &dfk) {
        double d1, d2;
        d_plus_minus(F, k, s, d1, d2);
        const double cdf = norm_cdf(phi * d2);
        fk = phi * scale * (k / F) * cdf - target;
        dfk = phi * scale * (cdf - phi * norm_pdf(d2) / s) / F;
    };
    // Starting point: the strike at which phi D Phi(phi d-) = target, that is
    // the premium-adjusted delta without its factor K/F (exact as s -> 0).
    const double a = std::fabs(target) / scale;
    const double x0 = (0 < a && a < 1) ? F * std::exp(-0.5 * (s * s) - phi * s * norm_ppf(a)) : kNaN;
    const double xtol = 1e-14 * F;
    if (phi < 0) {
        const double lo = F * std::exp(-10 * s - 1e-8);
        double hi = F;
        for (int i = 0; i < kMaxExpansions; ++i) {
            if (f(hi) <= 0) break;
            hi *= 2.0;
        }
        return solve_bracketed(fd, lo, hi, xtol, x0, K);
    }
    double k_star = 0.0;
    const int32_t st = pa_call_delta_maximiser_1(F, sigma, tau, k_star);
    if (st != QEF_OK) return st;
    if (f(k_star) < 0) return QEF_ERR_PA_MAX;
    double hi = k_star * 2.0;
    for (int i = 0; i < kMaxExpansions; ++i) {
        if (f(hi) <= 0) break;
        hi *= 2.0;
    }
    return solve_bracketed(fd, k_star, hi, xtol, x0, K);
}

// numpy's partial pairwise summation of n doubles at a stride (header comment).
double pairwise_sum(const double *a, int64_t n, int64_t stride) {
    if (n < 8) {
        double res = 0.0;
        for (int64_t i = 0; i < n; ++i) res += a[i * stride];
        return res;
    }
    if (n <= 128) {
        double r[8];
        for (int j = 0; j < 8; ++j) r[j] = a[j * stride];
        int64_t i = 8;
        for (; i < n - (n % 8); i += 8)
            for (int j = 0; j < 8; ++j) r[j] += a[(i + j) * stride];
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; ++i) res += a[i * stride];
        return res;
    }
    int64_t n2 = n / 2;
    n2 -= n2 % 8;
    return pairwise_sum(a, n2, stride) + pairwise_sum(a + n2 * stride, n - n2, stride);
}

// moments._simpson on v[0..n]: h/3 [v0 + vn + 4 sum(v[1:-1:2]) + 2 sum(v[2:-1:2])].
double simpson(const double *v, int64_t n, double h) {
    const double s_odd = pairwise_sum(v + 1, n / 2, 2);
    const double s_even = n >= 3 ? pairwise_sum(v + 2, (n - 1) / 2, 2) : 0.0;
    return h / 3.0 * (v[0] + v[n] + 4.0 * s_odd + 2.0 * s_even);
}

struct Smile {
    double F, tau, sqrt_tau;  // integration forward and expiry
    bool usd_base;
    double sF, stau, alpha, rho, nu;  // SABR smile
    bool same_forward;                // sF == F: the smile's ln(F/K) is the premium's, computed once
};

// h_k''(K) P(K) for k = 1, 2, 3 (moments.contract_weight times moments.otm_price).
inline void integrand(const Smile &c, double K, double &w1, double &w2, double &w3) {
    const double log_fk = std::log(c.F / K);
    const double vol = c.same_forward ? sabr_vol_log(log_fk, c.stau, c.alpha, c.rho, c.nu)
                                      : sabr_vol_1(K, c.sF, c.stau, c.alpha, c.rho, c.nu);
    const double phi = K < c.F ? -1.0 : 1.0;
    const double p = forward_premium_log(c.F, K, log_fk, vol, c.sqrt_tau, phi);
    const double u = std::log(K / c.F);
    const double y = c.usd_base ? -u : u;
    const double den = c.usd_base ? c.F * K : K * K;
    w1 = (-1.0 / den) * p;
    w2 = ((2.0 - 2.0 * y) / den) * p;
    w3 = ((6.0 * y - 3.0 * (y * y)) / den) * p;
}

// Integrand values on one side at the current level; n == 0 means none yet.
struct Side {
    std::vector<double> v[3];
    std::vector<double> scratch[3];
    int64_t n = 0;
};

// Fills side.v with the integrand at numpy.linspace(a, b, n + 1), reusing
// the previous level's values when n doubles it.
void fill_side(const Smile &c, double a, double b, int64_t n, Side &side) {
    const double delta = b - a;
    const double step = delta / static_cast<double>(n);
    const bool reuse = side.n > 0 && 2 * side.n == n && std::fabs(step) >= DBL_MIN;
    for (auto &w : side.scratch) w.resize(n + 1);
    for (int64_t i = 0; i <= n; ++i) {
        if (reuse && i % 2 == 0) {
            for (int k = 0; k < 3; ++k) side.scratch[k][i] = side.v[k][i / 2];
            continue;
        }
        double K;
        if (i == n) K = b;
        else if (step == 0.0) K = static_cast<double>(i) / static_cast<double>(n) * delta + a;  // numpy's denormal branch
        else K = static_cast<double>(i) * step + a;
        integrand(c, K, side.scratch[0][i], side.scratch[1][i], side.scratch[2][i]);
    }
    for (int k = 0; k < 3; ++k) side.v[k].swap(side.scratch[k]);
    side.n = n;
}

// The reference's integrate(n): left side, then right side, added to zeros.
void integrate(const Smile &c, double K_min, double K_max, int64_t n, Side sides[2], double out[3]) {
    double total[3] = {0.0, 0.0, 0.0};
    const double ends[2][2] = {{K_min, c.F}, {c.F, K_max}};
    for (int s = 0; s < 2; ++s) {
        const double a = ends[s][0], b = ends[s][1];
        fill_side(c, a, b, n, sides[s]);
        const double h = (b - a) / static_cast<double>(n);
        for (int k = 0; k < 3; ++k) total[k] += simpson(sides[s].v[k].data(), n, h);
    }
    for (int k = 0; k < 3; ++k) out[k] = total[k];
}

// float(np.max(np.abs(cur - prev) / np.maximum(np.abs(cur), tiny))), NaN-propagating.
double max_rel_change(const double cur[3], const double prev[3]) {
    double change = -kInf;
    for (int k = 0; k < 3; ++k) {
        const double a = std::fabs(cur[k]);
        const double den = std::isnan(a) ? a : std::max(a, DBL_MIN);
        const double r = std::fabs(cur[k] - prev[k]) / den;
        if (std::isnan(r)) return kNaN;
        change = std::max(change, r);
    }
    return change;
}

Smile make_smile(int64_t i, const In &F, const In &tau, const In &usd, const In &sF, const In &stau, const In &alpha,
                 const In &rho, const In &nu) {
    Smile c;
    c.F = F[i];
    c.tau = tau[i];
    c.sqrt_tau = std::sqrt(c.tau);
    c.usd_base = usd[i] != 0.0;
    c.sF = sF[i];
    c.stau = stau[i];
    c.alpha = alpha[i];
    c.rho = rho[i];
    c.nu = nu[i];
    c.same_forward = c.sF == c.F;
    return c;
}

}  // namespace

extern "C" {

int32_t qef_abi_version(void) { return QEF_KERNEL_ABI_VERSION; }

void qef_norm_cdf(int64_t n, const double *x, const int64_t *inc, double *out) {
    const In X{x, inc[0]};
    for (int64_t i = 0; i < n; ++i) out[i] = norm_cdf(X[i]);
}

void qef_norm_ppf(int64_t n, const double *p, const int64_t *inc, double *out) {
    const In P{p, inc[0]};
    for (int64_t i = 0; i < n; ++i) out[i] = norm_ppf(P[i]);
}

void qef_forward_premium(int64_t n, const double *F, const double *K, const double *sigma, const double *tau,
                         const double *phi, const int64_t *inc, double *out) {
    const In f{F, inc[0]}, k{K, inc[1]}, s{sigma, inc[2]}, t{tau, inc[3]}, p{phi, inc[4]};
    for (int64_t i = 0; i < n; ++i) out[i] = forward_premium_1(f[i], k[i], s[i], std::sqrt(t[i]), p[i]);
}

void qef_delta(int64_t n, const double *F, const double *K, const double *sigma, const double *tau, const double *phi,
               const double *df_base, const int64_t *inc, int32_t spot, int32_t premium_adjusted, double *out) {
    const In f{F, inc[0]}, k{K, inc[1]}, s{sigma, inc[2]}, t{tau, inc[3]}, p{phi, inc[4]}, d{df_base, inc[5]};
    for (int64_t i = 0; i < n; ++i)
        out[i] = delta_1(f[i], k[i], s[i], std::sqrt(t[i]), p[i], spot ? d[i] : 1.0, premium_adjusted != 0);
}

void qef_sabr_vol(int64_t n, const double *K, const double *F, const double *tau, const double *alpha,
                  const double *rho, const double *nu, const int64_t *inc, double *out) {
    const In k{K, inc[0]}, f{F, inc[1]}, t{tau, inc[2]}, a{alpha, inc[3]}, r{rho, inc[4]}, v{nu, inc[5]};
    for (int64_t i = 0; i < n; ++i) out[i] = sabr_vol_1(k[i], f[i], t[i], a[i], r[i], v[i]);
}

void qef_pa_call_delta_maximiser(int64_t n, const double *F, const double *sigma, const double *tau,
                                 const int64_t *inc, double *out, int32_t *status) {
    const In f{F, inc[0]}, s{sigma, inc[1]}, t{tau, inc[2]};
    for (int64_t i = 0; i < n; ++i) {
        double k_star = kNaN;
        status[i] = pa_call_delta_maximiser_1(f[i], s[i], t[i], k_star);
        out[i] = status[i] == QEF_OK ? k_star : kNaN;
    }
}

void qef_strike_from_delta(int64_t n, const double *target, const double *F, const double *sigma, const double *tau,
                           const double *phi, const double *df_base, const int64_t *inc, int32_t spot,
                           int32_t premium_adjusted, double *out, int32_t *status) {
    const In g{target, inc[0]}, f{F, inc[1]}, s{sigma, inc[2]}, t{tau, inc[3]}, p{phi, inc[4]}, d{df_base, inc[5]};
    for (int64_t i = 0; i < n; ++i) {
        double K = kNaN;
        status[i] = strike_from_delta_1(g[i], f[i], s[i], t[i], p[i], d[i], spot != 0, premium_adjusted != 0, K);
        out[i] = status[i] == QEF_OK ? K : kNaN;
    }
}

void qef_middle_contracts_fixed(int64_t m, const double *F, const double *K_min, const double *K_max,
                                const double *tau, const double *usd_base, const double *smile_F,
                                const double *smile_tau, const double *alpha, const double *rho, const double *nu,
                                const int64_t *inc, int64_t n, double *values, int32_t *status) {
    const In f{F, inc[0]}, lo{K_min, inc[1]}, hi{K_max, inc[2]}, t{tau, inc[3]}, u{usd_base, inc[4]},
        sf{smile_F, inc[5]}, st{smile_tau, inc[6]}, a{alpha, inc[7]}, r{rho, inc[8]}, v{nu, inc[9]};
    Side sides[2];
    for (int64_t i = 0; i < m; ++i) {
        double *out = values + 3 * i;
        out[0] = out[1] = out[2] = kNaN;
        const Smile c = make_smile(i, f, t, u, sf, st, a, r, v);
        if (!(lo[i] < c.F && c.F < hi[i])) {
            status[i] = QEF_ERR_MIDDLE;
            continue;
        }
        if (n < 1 || n > kMaxSubintervals) {
            status[i] = QEF_ERR_SUBINTERVALS;
            continue;
        }
        sides[0].n = sides[1].n = 0;
        try {  // no C++ exception may cross the C interface
            integrate(c, lo[i], hi[i], n, sides, out);
        } catch (const std::bad_alloc &) {
            status[i] = QEF_ERR_MEMORY;
            continue;
        }
        status[i] = QEF_OK;
    }
}

void qef_middle_contracts_adaptive(int64_t m, const double *F, const double *K_min, const double *K_max,
                                   const double *tau, const double *usd_base, const double *smile_F,
                                   const double *smile_tau, const double *alpha, const double *rho, const double *nu,
                                   const int64_t *inc, int64_t n0, double rtol, int64_t n_max, double *values,
                                   int64_t *n_out, double *rel_change, int32_t *converged, int32_t *status) {
    const In f{F, inc[0]}, lo{K_min, inc[1]}, hi{K_max, inc[2]}, t{tau, inc[3]}, u{usd_base, inc[4]},
        sf{smile_F, inc[5]}, st{smile_tau, inc[6]}, a{alpha, inc[7]}, r{rho, inc[8]}, v{nu, inc[9]};
    Side sides[2];
    for (int64_t i = 0; i < m; ++i) {
        double *out = values + 3 * i;
        out[0] = out[1] = out[2] = kNaN;
        n_out[i] = 0;
        rel_change[i] = kNaN;
        converged[i] = 0;
        const Smile c = make_smile(i, f, t, u, sf, st, a, r, v);
        if (!(lo[i] < c.F && c.F < hi[i])) {
            status[i] = QEF_ERR_MIDDLE;
            continue;
        }
        // n doubles from n0 until n >= n_max, so every level is at most
        // 2 max(n0, n_max) <= 2^23 (2 n0 when n0 >= n_max).
        if (n0 < 1 || n0 > kMaxSubintervals || n_max > kMaxSubintervals) {
            status[i] = QEF_ERR_SUBINTERVALS;
            continue;
        }
        sides[0].n = sides[1].n = 0;
        double prev[3], cur[3];
        int64_t n = n0;
        try {  // no C++ exception may cross the C interface
            integrate(c, lo[i], hi[i], n, sides, prev);
            for (;;) {
                n *= 2;
                integrate(c, lo[i], hi[i], n, sides, cur);
                const double change = max_rel_change(cur, prev);
                if (change < rtol || n >= n_max) {
                    for (int k = 0; k < 3; ++k) out[k] = cur[k];
                    n_out[i] = n;
                    rel_change[i] = change;
                    converged[i] = change < rtol;
                    break;
                }
                for (int k = 0; k < 3; ++k) prev[k] = cur[k];
            }
        } catch (const std::bad_alloc &) {
            status[i] = QEF_ERR_MEMORY;
            continue;
        }
        status[i] = QEF_OK;
    }
}

}  // extern "C"
