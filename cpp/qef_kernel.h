/*
 * C interface of the C++ kernel (research design, Stage 5), loaded from Python
 * by qef.kernel through ctypes. qef_kernel.cpp documents the formulas, the
 * Python functions each routine mirrors and the sources.
 *
 * Array conventions. Every routine works element by element on arrays of
 * length n (m for the middle part, one element per smile). Each double input
 * x has an increment inc[j] (j its position among the double inputs), 0 or 1:
 * element i reads x[i * inc[j]], so a scalar is passed as one element with
 * increment 0. Outputs are contiguous arrays of length n (the middle part
 * writes m rows of three contracts). Flags are int32 (0 false, otherwise true).
 *
 * Status codes (int32, one per element or smile) name the exception that the
 * Python reference raises for the same input; qef.kernel raises it.
 */
#ifndef QEF_KERNEL_H
#define QEF_KERNEL_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define QEF_KERNEL_ABI_VERSION 1

enum qef_status {
    QEF_OK = 0,
    QEF_ERR_SIGN = 1,           /* ValueError: target sign does not match the option type, or target is zero */
    QEF_ERR_VOL = 2,            /* ValueError: volatility or time to expiry not positive and finite */
    QEF_ERR_PIPS = 3,           /* ValueError: pips delta outside (0, D_b) */
    QEF_ERR_PA_MAX = 4,         /* ValueError: target above the maximal premium-adjusted call delta */
    QEF_ERR_SAME_SIGN = 5,      /* ValueError: no sign change on the bracket (scipy brentq) */
    QEF_ERR_NAN = 6,            /* ValueError: NaN function value in the root search (scipy brentq) */
    QEF_ERR_XTOL = 7,           /* ValueError: absolute root tolerance not positive (scipy brentq) */
    QEF_ERR_NO_CONVERGENCE = 8, /* RuntimeError: root search did not converge */
    QEF_ERR_MIDDLE = 9,         /* ValueError: not K_min < F < K_max */
    QEF_ERR_SUBINTERVALS = 10,  /* ValueError: n or n0 not in [1, 2^22], or n_max above 2^22 */
    QEF_ERR_MEMORY = 11         /* MemoryError: the node buffers could not be allocated */
};

int32_t qef_abi_version(void);

/* Standard normal distribution function and its inverse (inputs: x or p). */
void qef_norm_cdf(int64_t n, const double *x, const int64_t *inc, double *out);
void qef_norm_ppf(int64_t n, const double *p, const int64_t *inc, double *out);

/* qef.fx.gk.forward_premium. Inputs: F, K, sigma, tau, phi. */
void qef_forward_premium(int64_t n, const double *F, const double *K, const double *sigma, const double *tau,
                         const double *phi, const int64_t *inc, double *out);

/* qef.fx.gk.delta. Inputs: F, K, sigma, tau, phi, df_base (ignored unless spot). */
void qef_delta(int64_t n, const double *F, const double *K, const double *sigma, const double *tau, const double *phi,
               const double *df_base, const int64_t *inc, int32_t spot, int32_t premium_adjusted, double *out);

/* qef.fx.sabr.sabr_vol with beta = 1. Inputs: K, F, tau, alpha, rho, nu. */
void qef_sabr_vol(int64_t n, const double *K, const double *F, const double *tau, const double *alpha,
                  const double *rho, const double *nu, const int64_t *inc, double *out);

/* qef.fx.gk.pa_call_delta_maximiser. Inputs: F, sigma, tau. */
void qef_pa_call_delta_maximiser(int64_t n, const double *F, const double *sigma, const double *tau,
                                 const int64_t *inc, double *out, int32_t *status);

/* qef.fx.gk.strike_from_delta. Inputs: target, F, sigma, tau, phi, df_base. */
void qef_strike_from_delta(int64_t n, const double *target, const double *F, const double *sigma, const double *tau,
                           const double *phi, const double *df_base, const int64_t *inc, int32_t spot,
                           int32_t premium_adjusted, double *out, int32_t *status);

/*
 * Middle part of qef.fx.moments.middle_contracts for the SABR (beta = 1) smile
 * sabr_vol(K, smile_F, smile_tau, alpha, rho, nu). Inputs, in order: F, K_min,
 * K_max, tau, usd_base (non-zero for a USD-base pair), smile_F, smile_tau,
 * alpha, rho, nu. values has m rows of the contracts k = 1, 2, 3.
 *
 * _fixed: composite Simpson's rule at n subintervals per side (the reference's
 * inner integrate(n)). _adaptive: the reference's doubling loop from n0 until
 * the largest relative change is below rtol or n >= n_max.
 */
void qef_middle_contracts_fixed(int64_t m, const double *F, const double *K_min, const double *K_max,
                                const double *tau, const double *usd_base, const double *smile_F,
                                const double *smile_tau, const double *alpha, const double *rho, const double *nu,
                                const int64_t *inc, int64_t n, double *values, int32_t *status);
void qef_middle_contracts_adaptive(int64_t m, const double *F, const double *K_min, const double *K_max,
                                   const double *tau, const double *usd_base, const double *smile_F,
                                   const double *smile_tau, const double *alpha, const double *rho, const double *nu,
                                   const int64_t *inc, int64_t n0, double rtol, int64_t n_max, double *values,
                                   int64_t *n_out, double *rel_change, int32_t *converged, int32_t *status);

#ifdef __cplusplus
}
#endif

#endif /* QEF_KERNEL_H */
