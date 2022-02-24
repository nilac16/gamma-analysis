/** Gamma 2D analysis
 * 
 *  ANSI C implementation of the gamma analysis algorithm commonly used in 
 *  medical physics to quantify the similarity of two dose distributions.
 * 
 *  Copyright (C) 2022 Calin M Reamy
 *  
 *  THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS 
 *  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
 *  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 *  DEALINGS IN THE SOFTWARE.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include <stb/stb_image_write.h>
#include "gamma.h"

#define LINE_SEARCH_THRESHOLD 0.000001

/* Explicit typecasting, conspicuously stolen from C++ */
#define STATIC_CAST(type, expr) (type)(expr)

#if __STDC_VERSION__ < 199901L
#   if _MSC_VER
#       define inline __inline
#   else
#       define inline __inline__
#   endif
/* I think pretty much every compiler supports vendor-extended restricted
pointer semantics */
#   define restrict __restrict
#define fma(a, b, c) (((a) * (b)) + (c))
#endif

/* Array declarator qualifiers starting C11? */
#if __STDC_VERSION__ < 201112L
#   define _q(qualifiers)
#else
#   define _q(qualifiers) qualifiers
#endif


/** An inner product using (diagonal) metric tensor S */
static double gamma2d_metric(const double X[_q(static 2)],
                             const double S[_q(static 2)])
{
    return S[0] * X[0] * X[0] + S[1] * X[1] * X[1];
}

/** Floors vec into zs */
static void gamma2d_floorz_vec(const double vec[_q(static 2)],
                               long zs[_q(static 2)])
{
    zs[0] = STATIC_CAST(long, floor(vec[0]));
    zs[1] = STATIC_CAST(long, floor(vec[1]));
}

static inline void gamma2d_add_vec(double augend[_q(static 2)],
                                   const double addend[_q(static 2)])
{
    augend[0] += addend[0];
    augend[1] += addend[1];
}

static inline void gamma2d_sub_vec(double minuend[_q(static 2)],
                                   const double subtrahend[_q(static 2)])
{
    minuend[0] -= subtrahend[0];
    minuend[1] -= subtrahend[1];
}

static inline void gamma2d_mul_vec(double dst[_q(static 2)],
                                   const double src[_q(static 2)])
{
    dst[0] *= src[0];
    dst[1] *= src[1];
}

/** Fused-multiply add: a * b + c */
static inline void gamma2d_fma_vec(double a[_q(static 2)],
                                   const double b[_q(static 2)],
                                   const double c[_q(static 2)])
{
    a[0] = fma(a[0], b[0], c[0]);
    a[1] = fma(a[1], b[1], c[1]);
}

/** Casts z1 and z2 to doubles and stores them in order in dst */
static inline void gamma2d_loadz_vec(double dst[_q(static 2)], long z1, long z2)
{
    dst[0] = STATIC_CAST(double, z1);
    dst[1] = STATIC_CAST(double, z2);
}

/** Loads dst with doubles x1 and x2, in order */
static inline void gamma2d_loadpd_vec(double dst[_q(static 2)], double x1, double x2)
{
    dst[0] = x1;
    dst[1] = x2;
}

static inline void gamma2d_copy_vec(double dst[_q(static 2)], const double src[_q(static 2)])
{
    memcpy(dst, src, sizeof *dst * 2);
}

/** Initializes the scaling vector for converting the dose coordinates to lattice coordinates */
static void gamma2d_init_transform(double dst[_q(static 2)],
                                   const struct dose_distribution *restrict meas,
                                   const struct dose_distribution *restrict calc)
{
    dst[0] = meas->px_spacing[0] / calc->px_spacing[0];
    dst[1] = meas->px_spacing[1] / calc->px_spacing[1];
}

/** Initializes the offset vector for dose coordinate conversion */
static void gamma2d_init_shift(double dst[_q(static 2)],
                               const struct dose_distribution *restrict meas,
                               const struct dose_distribution *restrict calc)
{
    dst[0] = (meas->top_left[0] - calc->top_left[0]) / calc->px_spacing[0];
    dst[1] = (meas->top_left[1] - calc->top_left[1]) / calc->px_spacing[1];
}

static int gamma2d_square_in_radius(const double Xm[_q(static 2)],
                                    const double S[_q(static 2)],
                                    double rmax)
{
    double r00, r10, r01, r11;
    r00 = gamma2d_metric(Xm, S);
    if (r00 < rmax) {
        return 1;
    }
    r10 = fma(S[0], 1 + 2 * Xm[0], r00);
    if (r10 < rmax) {
        return 1;
    }
    r01 = fma(S[1], 1 + 2 * Xm[1], r00);
    if (r01 < rmax) {
        return 1;
    }
    r11 = r01 + r10 - r00;
    return r11 < rmax;
}

static double gamma2d_dist_bounded_access(const struct dose_distribution *dist,
                                          long i, long j)
{
    /* Index in mathematically unsupported region */
    const int unsup_x = (i < 0) || (i >= dist->px_dim[0]);
    const int unsup_y = (j < 0) || (j >= dist->px_dim[1]);
    if (unsup_x || unsup_y) {
        return 0.0;
    } else {
        return dist->data[j * dist->px_dim[0] + i];
    }
}

static void gamma2d_load_interpolant(const struct gamma_analysis *gamma,
                                     const struct dose_distribution *calc,
                                     const long cpoint[_q(static 2)],
                                     double interp[_q(static 4)],
                                     double mdose)
{
    double mult = gamma->dta / gamma->diff;
    if (gamma->normalization_mode == GAMMA_NORM_LOCAL) {
        mult /= mdose;
    }
    interp[0] = gamma2d_dist_bounded_access(calc, cpoint[0], cpoint[1]);
    interp[1] = gamma2d_dist_bounded_access(calc, cpoint[0] + 1, cpoint[1]);
    interp[2] = gamma2d_dist_bounded_access(calc, cpoint[0], cpoint[1] + 1);
    interp[3] = gamma2d_dist_bounded_access(calc, cpoint[0] + 1, cpoint[1] + 1);
    interp[1] -= interp[0];
    interp[3] -= interp[1] + interp[2];
    interp[2] -= interp[0];
    interp[0] -= mdose;
    interp[0] *= mult;
    interp[1] *= mult;
    interp[2] *= mult;
    interp[3] *= mult;
}

static double gamma2d_objective_function(const double X[_q(static 2)],
                                         const double S[_q(static 2)],
                                         const double Xm[_q(static 2)],
                                         const double interp[_q(static 4)])
{
    double r[2], rsqr;
    const double f = fma(
        X[1],
        fma(X[0], interp[3], interp[2]),
        fma(X[0], interp[1], interp[0])
    );
    gamma2d_copy_vec(r, X);
    gamma2d_add_vec(r, Xm);
    rsqr = gamma2d_metric(r, S);
    return rsqr + f * f;
}                                         

static double gamma2d_objective_xpartial_zero(const double X[_q(static 2)],
                                              const double S[_q(static 2)],
                                              const double Xm[_q(static 2)],
                                              const double interp[_q(static 4)])
{
    const double inter1 = fma(X[1], interp[3], interp[1]);
    const double inter2 = fma(X[1], interp[2], interp[0]);
    const double numer = S[0] * Xm[0] + inter1 * inter2;
    const double denom = S[0] + inter1 * inter1;
    return -numer / denom;
}

static double gamma2d_objective_ypartial_zero(const double X[_q(static 2)],
                                              const double S[_q(static 2)],
                                              const double Xm[_q(static 2)],
                                              const double interp[_q(static 4)])
{
    const double inter1 = fma(X[0], interp[3], interp[2]);
    const double inter2 = fma(X[0], interp[1], interp[0]);
    const double numer = S[1] * Xm[1] + inter1 * inter2;
    const double denom = S[1] + inter1 * inter1;
    return -numer / denom;
}

static void gamma2d_clamp_unity(double *x)
{
    if (*x < 0.0) {
        *x = 0.0;
    } else if (*x > 1.0) {
        *x = 1.0;
    }
}

static void gamma2d_line_search_move(double X[_q(static 2)],
                                     const double S[_q(static 2)],
                                     const double Xm[_q(static 2)],
                                     const double interp[_q(static 2)])
{
    X[0] = gamma2d_objective_xpartial_zero(X, S, Xm, interp);
    gamma2d_clamp_unity(X);
    X[1] = gamma2d_objective_ypartial_zero(X, S, Xm, interp);
    gamma2d_clamp_unity(X + 1);
}

static double gamma2d_line_search(const double S[_q(static 2)],
                                  const double Xm[_q(static 2)],
                                  const double interp[_q(static 4)])
{
    double X[2];
    double obj, last;
    gamma2d_loadpd_vec(X, 0.0, 0.0);
    obj = gamma2d_objective_function(X, S, Xm, interp);
    do {
        last = obj;
        gamma2d_line_search_move(X, S, Xm, interp);
        obj = gamma2d_objective_function(X, S, Xm, interp);
    } while ((last - obj) > LINE_SEARCH_THRESHOLD);
    return obj;
}

static inline void gamma2d_valcmp(double val, double *tval)
{
    if (val < *tval) {
        *tval = val;
    }
}

static int gamma2d_minimize_square(const struct gamma_analysis *gamma,
                                   const struct dose_distribution *calc,
                                   double *rmax,
                                   const long cpoint[_q(static 2)],
                                   const double r_m[_q(static 2)],
                                   double mdose)
{
    double Xm[2];
    double interp[4];
    double S[2];
    gamma2d_loadz_vec(Xm, cpoint[0], cpoint[1]);
    gamma2d_sub_vec(Xm, r_m);
    gamma2d_copy_vec(S, calc->px_spacing);
    gamma2d_mul_vec(S, S);
    if (!gamma2d_square_in_radius(Xm, S, *rmax)) {
        return 0;
    }
    gamma2d_load_interpolant(gamma, calc, cpoint, interp, mdose);
    gamma2d_valcmp(gamma2d_line_search(S, Xm, interp), rmax);
    return 1;
}

static int gamma2d_fan_horizontal(const struct gamma_analysis *gamma,
                                  const struct dose_distribution *calc,
                                  double *rmax,
                                  long cpoint[_q(static 2)],
                                  const double r_m[_q(static 2)],
                                  double mdose)
{
    int n = 0, aug;
    const long x_c = cpoint[0];
    long offset = 0;
    aug = gamma2d_minimize_square(gamma, calc, rmax, cpoint, r_m, mdose);
    if (!aug) {
        return 0;
    }
    do {
        n += aug;
        offset++;
        cpoint[0] = x_c + offset;
        aug = gamma2d_minimize_square(gamma, calc, rmax, cpoint, r_m, mdose);
        cpoint[0] = x_c - offset;
        aug += gamma2d_minimize_square(gamma, calc, rmax, cpoint, r_m, mdose);
    } while (aug);
    cpoint[0] = x_c;
    return n;
}

static void gamma2d_fan_vertical(const struct gamma_analysis *gamma,
                                 const struct dose_distribution *calc,
                                 double *rmax,
                                 long cpoint[_q(static 2)],
                                 const double r_m[_q(static 2)],
                                 double mdose)
{
    int n;
    const long y_c = cpoint[1];
    long offset = 0;
    gamma2d_fan_horizontal(gamma, calc, rmax, cpoint, r_m, mdose);
    do {
        offset++;
        cpoint[1] = y_c + offset;
        n = gamma2d_fan_horizontal(gamma, calc, rmax, cpoint, r_m, mdose);
        cpoint[1] = y_c - offset;
        n += gamma2d_fan_horizontal(gamma, calc, rmax, cpoint, r_m, mdose);
    } while (n);
}

static double gamma2d_compute_point(const struct gamma_analysis *gamma,
                                    const struct dose_distribution *calc,
                                    const double r_m[_q(static 2)],
                                    double mdose)
{
    long origin[2];
    double rmax = HUGE_VAL;
    gamma2d_floorz_vec(r_m, origin);
    gamma2d_fan_vertical(gamma, calc, &rmax, origin, r_m, mdose);
    return sqrt(rmax / (gamma->dta * gamma->dta));
}

static int gamma2d_fan_horizontal_pass_only(const struct gamma_analysis *gamma,
                                            const struct dose_distribution *calc,
                                            double *rmax,
                                            long cpoint[_q(static 2)],
                                            const double r_m[_q(static 2)],
                                            double mdose,
                                            jmp_buf env)
{
    const double rnorm = gamma->dta * gamma->dta;
    int n = 0, aug;
    const long x_c = cpoint[0];
    long offset = 0;
    aug = gamma2d_minimize_square(gamma, calc, rmax, cpoint, r_m, mdose);
    if (!aug) {
        return 0;
    }
    if (*rmax < rnorm) {
        longjmp(env, 1);
    }
    do {
        n += aug;
        offset++;
        cpoint[0] = x_c + offset;
        aug = gamma2d_minimize_square(gamma, calc, rmax, cpoint, r_m, mdose);
        if (*rmax < rnorm) {
            longjmp(env, 1);
        }
        cpoint[0] = x_c - offset;
        aug += gamma2d_minimize_square(gamma, calc, rmax, cpoint, r_m, mdose);
        if (*rmax < rnorm) {
            longjmp(env, 1);
        }
    } while (aug);
    cpoint[0] = x_c;
    return n;
}

static void gamma2d_fan_vertical_pass_only(const struct gamma_analysis *gamma,
                                           const struct dose_distribution *calc,
                                           double *rmax,
                                           long cpoint[_q(static 2)],
                                           const double r_m[_q(static 2)],
                                           double mdose,
                                           jmp_buf env)
{
    int n;
    const long y_c = cpoint[1];
    long offset = 0;
    gamma2d_fan_horizontal_pass_only(gamma, calc, rmax, cpoint, r_m, mdose, env);
    do {
        offset++;
        cpoint[1] = y_c + offset;
        n = gamma2d_fan_horizontal_pass_only(gamma, calc, rmax, cpoint, r_m, mdose, env);
        cpoint[1] = y_c - offset;
        n += gamma2d_fan_horizontal_pass_only(gamma, calc, rmax, cpoint, r_m, mdose, env);
    } while (n);
}

/** Uses a setjmp call to jump back here as soon as it finds a passing value
 */
static int gamma2d_point_passes(const struct gamma_analysis *gamma,
                                const struct dose_distribution *calc,
                                const double r_m[_q(static 2)],
                                double mdose)
/** There is currently no need to access locals after a setjmp invocation, 
 *  so they are NOT declared volatile.
 *  volatile is commented below to remind me of this if I ever decide to 
 *  change this behavior
 */
{
    /* volatile */ long origin[2];
    /* volatile */ double rmax = HUGE_VAL;
    /* volatile */ jmp_buf env;
    gamma2d_floorz_vec(r_m, origin);
    if (setjmp(env)) {
        return 1;
    }
    gamma2d_fan_vertical_pass_only(gamma, calc, &rmax, origin, r_m, mdose, env);
    return 0;
}

static void gamma2d_compute_dist(struct gamma_analysis *gamma,
                                 const struct dose_distribution *restrict meas,
                                 const struct dose_distribution *restrict calc,
                                 double *gptr)
{
    long i, j;
    const double *mptr = meas->data;
    double trns[2], shft[2], a[2];
    gamma2d_init_transform(trns, meas, calc);
    gamma2d_init_shift(shft, meas, calc);
    for (j = 0; j < meas->px_dim[1]; j++) {
        for (i = 0; i < meas->px_dim[0]; i++, mptr++, gptr++) {
            if (*mptr < gamma->threshold) {
                *gptr = GAMMA_SIG_VAL;
                continue;
            }
            gamma2d_loadz_vec(a, i, j);
            gamma2d_fma_vec(a, trns, shft);
            *gptr = gamma2d_compute_point(gamma, calc, a, *mptr);
        }
    }
}

static double gamma2d_compute_pass_rate(struct gamma_analysis *gamma,
                                        const struct dose_distribution *restrict meas,
                                        const struct dose_distribution *restrict calc)
{
    long i, j;
    const double *mptr = meas->data;
    double trns[2], shft[2], a[2];
    long passed = 0, npts = 0;
    gamma2d_init_transform(trns, meas, calc);
    gamma2d_init_shift(shft, meas, calc);
    for (j = 0; j < meas->px_dim[1]; j++) {
        for (i = 0; i < meas->px_dim[0]; i++, mptr++) {
            if (*mptr < gamma->threshold) {
                continue;
            }
            gamma2d_loadz_vec(a, i, j);
            gamma2d_fma_vec(a, trns, shft);
            if (gamma2d_point_passes(gamma, calc, a, *mptr)) {
                passed++;
            }
            npts++;
        }
    }
    if (!npts) {
        return GAMMA_SIG_VAL;
    } else {
        return STATIC_CAST(double, passed) / STATIC_CAST(double, npts);
    }
}

static double gamma2d_find_max(long N, const double dptr[_q(static N)])
{
    double x = -HUGE_VAL;
    do {
        N--;
        if (*dptr > x) {
            x = *dptr;
        }
        dptr++;
    } while (N);
    return x;
}

static double gamma2d_pass_rate_from_dist(long N, const double gptr[_q(static N)])
{
    long passed = 0, npts = 0;
    do {
        double g = *gptr;
        if (g != GAMMA_SIG_VAL) {
            if (g < 1.0) {
                passed++;
            }
            npts++;
        }
        gptr++;
        N--;
    } while (N);
    if (!npts) {
        return GAMMA_SIG_VAL;
    } else {
        return STATIC_CAST(double, passed) / STATIC_CAST(double, npts);
    }
}

double gamma2d_compute(struct gamma_analysis *gamma,
                       const struct dose_distribution *restrict meas,
                       const struct dose_distribution *restrict calc,
                       double *gdist)
{
    const double threshold = gamma->threshold;
    const double diff = gamma->diff;
    const double calcmax = gamma2d_find_max(calc->px_dim[0] * calc->px_dim[1], calc->data);
    double pass_rate;

    if (!calcmax) {
        return GAMMA_SIG_VAL;
    }

    if (gamma->normalization_mode == GAMMA_NORM_GLOBAL) {
        gamma->diff *= calcmax;
    }
    gamma->threshold *= calcmax;

    if (gdist) {
        gamma2d_compute_dist(gamma, meas, calc, gdist);
        pass_rate = gamma2d_pass_rate_from_dist(meas->px_dim[0] * meas->px_dim[1], gdist);
    } else {
        pass_rate = gamma2d_compute_pass_rate(gamma, meas, calc);
    }

    gamma->diff = diff;
    gamma->threshold = threshold;
    return pass_rate;
}
