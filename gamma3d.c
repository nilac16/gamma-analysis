/** Gamma 3D analysis
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
#include "gamma.h"

#define LINE_SEARCH_THRESHOLD 0.000001

/* Explicit typecasting, conspicuously stolen from C++ */
#define STATIC_CAST(type, expr) (type)(expr)

#if __STDC_VERSION__ < 199901L
#   if MSVC
#       define inline __inline
#   else
#       define inline __inline__
#   endif
#   define restrict __restrict
#define fma(a, b, c) (((a) * (b)) + (c))
#define round(x) floor((x) + 0.5)
#endif

/* Array declarator qualifiers starting C11? */
#if __STDC_VERSION__ < 201112L
#   define _q(qualifiers)
#else
#   define _q(qualifiers) qualifiers
#endif


/** An inner product using (diagonal) metric tensor S */
static double gamma3d_metric(const double X[_q(static 3)],
                             const double S[_q(static 3)])
{
    return S[0] * X[0] * X[0] + S[1] * X[1] * X[1] + S[2] * X[2] * X[2];
}

/** Floors vec into zs */
static void gamma3d_vec_floorz(const double vec[_q(static 3)],
                               long zs[_q(static 3)])
{
    zs[0] = STATIC_CAST(long, floor(vec[0]));
    zs[1] = STATIC_CAST(long, floor(vec[1]));
    zs[2] = STATIC_CAST(long, floor(vec[2]));
}

static inline void gamma3d_vec_add(double augend[_q(static 3)],
                                   const double addend[_q(static 3)])
{
    augend[0] += addend[0];
    augend[1] += addend[1];
    augend[2] += addend[2];
}

static inline void gamma3d_vec_sub(double minuend[_q(static 3)],
                                   const double subtrahend[_q(static 3)])
{
    minuend[0] -= subtrahend[0];
    minuend[1] -= subtrahend[1];
    minuend[2] -= subtrahend[2];
}

static inline void gamma3d_vec_mul(double dst[_q(static 3)],
                                   const double src[_q(static 3)])
{
    dst[0] *= src[0];
    dst[1] *= src[1];
    dst[2] *= src[2];
}

/** Fused-multiply add: a * b + c */
static inline void gamma3d_vec_fma(double a[_q(static 3)],
                                   const double b[_q(static 3)],
                                   const double c[_q(static 3)])
{
    a[0] = fma(a[0], b[0], c[0]);
    a[1] = fma(a[1], b[1], c[1]);
    a[2] = fma(a[2], b[2], c[2]);
}

/** Casts z1 and z2 to doubles and stores them in order in dst */
static inline void gamma3d_vec_loadz(double dst[_q(static 3)], long z1, long z2, long z3)
{
    dst[0] = STATIC_CAST(double, z1);
    dst[1] = STATIC_CAST(double, z2);
    dst[2] = STATIC_CAST(double, z3);
}

static inline void gamma3d_vec_copy(double dst[_q(static 3)], const double src[_q(static 3)])
{
    memcpy(dst, src, sizeof *dst * 3);
}

/** Initializes the scaling vector for converting the dose coordinates to lattice coordinates */
static void gamma3d_init_transform(double dst[_q(static 3)],
                                   const struct dose_distribution *restrict meas,
                                   const struct dose_distribution *restrict calc)
{
    dst[0] = meas->px_spacing[0] / calc->px_spacing[0];
    dst[1] = meas->px_spacing[1] / calc->px_spacing[1];
    dst[2] = meas->px_spacing[2] / calc->px_spacing[2];
}

/** Initializes the offset vector for dose coordinate conversion */
static void gamma3d_init_shift(double dst[_q(static 3)],
                               const struct dose_distribution *restrict meas,
                               const struct dose_distribution *restrict calc)
{
    dst[0] = (meas->top_left[0] - calc->top_left[0]) / calc->px_spacing[0];
    dst[1] = (meas->top_left[1] - calc->top_left[1]) / calc->px_spacing[1];
    dst[2] = (meas->top_left[2] - calc->top_left[2]) / calc->px_spacing[2];
}

/** There has to be a better way to do this */
static int gamma3d_cube_in_radius(const double Xm[_q(static 3)],
                                  const double S[_q(static 3)],
                                  double rmax)
{
    double r000, rtest1, rtest2, rtest3, clbr;
    r000 = gamma3d_metric(Xm, S);
    if (r000 < rmax) {
        return 1;
    }
    rtest1 = fma(S[0], 2 * Xm[0] + 1, r000);
    if (rtest1 < rmax) {
        return 1;
    }
    rtest2 = fma(S[1], 2 * Xm[1] + 1, r000);
    if (rtest2 < rmax) {
        return 1;
    }
    rtest3 = fma(S[2], 2 * Xm[2] + 1, r000);
    if (rtest3 < rmax) {
        return 1;
    }
    clbr = rtest1;
    rtest1 += rtest2 - r000;
    if (rtest1 < rmax) {
        return 1;
    }
    rtest2 += rtest3 - r000;
    if (rtest2 < rmax) {
        return 1;
    }
    rtest3 += clbr - r000;
    if (rtest3 < rmax) {
        return 1;
    }
    rtest2 += clbr - r000;
    return rtest2 < rmax;
}

static double gamma3d_dist_bounded_access(const struct dose_distribution *dist,
                                          long i, long j, long k)
{
    /* Index in mathematically unsupported region */
    const int unsup_x = (i < 0) || (i >= dist->px_dim[0]);
    const int unsup_y = (j < 0) || (j >= dist->px_dim[1]);
    const int unsup_z = (k < 0) || (k >= dist->px_dim[2]);
    if (unsup_x || unsup_y || unsup_z) {
        return 0.0;
    } else {
        return dist->data[(k * dist->px_dim[1] + j) * dist->px_dim[0] + i];
    }
}

static void gamma3d_load_interpolant(const struct gamma_analysis *gamma,
                                     const struct dose_distribution *calc,
                                     const long cpoint[_q(static 3)],
                                     double interp[_q(static 8)],
                                     double mdose)
{
    double mult = gamma->dta / gamma->diff;
    if (gamma->normalization_mode == GAMMA_NORM_LOCAL) {
        mult /= mdose;
    }
    interp[0] = gamma3d_dist_bounded_access(calc, cpoint[0],     cpoint[1],     cpoint[2]);
    interp[1] = gamma3d_dist_bounded_access(calc, cpoint[0] + 1, cpoint[1],     cpoint[2]);
    interp[2] = gamma3d_dist_bounded_access(calc, cpoint[0],     cpoint[1] + 1, cpoint[2]);
    interp[3] = gamma3d_dist_bounded_access(calc, cpoint[0] + 1, cpoint[1] + 1, cpoint[2]);
    interp[4] = gamma3d_dist_bounded_access(calc, cpoint[0],     cpoint[1],     cpoint[2] + 1);
    interp[5] = gamma3d_dist_bounded_access(calc, cpoint[0] + 1, cpoint[1],     cpoint[2] + 1);
    interp[6] = gamma3d_dist_bounded_access(calc, cpoint[0],     cpoint[1] + 1, cpoint[2] + 1);
    interp[7] = gamma3d_dist_bounded_access(calc, cpoint[0] + 1, cpoint[1] + 1, cpoint[2] + 1);

    /* Just trust me dude */
    interp[1] -= interp[0];
    interp[3] -= interp[2] + interp[1];
    interp[7] -= interp[3] - interp[4] + interp[5] + interp[6];
    interp[5] -= interp[4] + interp[1];
    interp[6] -= interp[2] + interp[1];
    interp[2] -= interp[0];
    interp[4] -= interp[0];
    interp[0] -= mdose;

    interp[0] *= mult;
    interp[1] *= mult;
    interp[2] *= mult;
    interp[3] *= mult;
    interp[4] *= mult;
    interp[5] *= mult;
    interp[6] *= mult;
    interp[7] *= mult;
}

static double gamma3d_objective_function(const double X[_q(static 3)],
                                         const double S[_q(static 3)],
                                         const double Xm[_q(static 3)],
                                         const double interp[_q(static 8)])
{
    double r[3], rsqr;
    const double f = fma(
        X[2],
        fma(
            X[1],
            fma(X[0], interp[7], interp[6]),
            fma(X[0], interp[5], interp[4])
        ),
        fma(
            X[1],
            fma(X[0], interp[3], interp[2]),
            fma(X[0], interp[1], interp[0])
        )
    );
    gamma3d_vec_copy(r, X);
    gamma3d_vec_add(r, Xm);
    rsqr = gamma3d_metric(r, S);
    return rsqr + f * f;
}                                         

static double gamma3d_objective_xpartial_zero(const double X[_q(static 3)],
                                              const double S[_q(static 3)],
                                              const double Xm[_q(static 3)],
                                              const double interp[_q(static 8)])
{
    const double deriv = fma(
        X[2],
        fma(X[1], interp[7], interp[5]),
        fma(X[1], interp[3], interp[1])
    );
    const double comp = fma(
        X[2],
        fma(X[1], interp[6], interp[4]),
        fma(X[1], interp[2], interp[0])
    );
    const double numer = S[0] * Xm[0] + comp * deriv;
    const double denom = S[0] + deriv * deriv;
    return -numer / denom;
}

static double gamma3d_objective_ypartial_zero(const double X[_q(static 3)],
                                              const double S[_q(static 3)],
                                              const double Xm[_q(static 3)],
                                              const double interp[_q(static 8)])
{
    const double deriv = fma(
        X[2],
        fma(X[0], interp[7], interp[6]),
        fma(X[0], interp[3], interp[2])
    );
    const double comp = fma(
        X[2],
        fma(X[1], interp[5], interp[4]),
        fma(X[1], interp[1], interp[0])
    );
    const double numer = S[1] * Xm[1] + comp * deriv;
    const double denom = S[1] + deriv * deriv;
    return -numer / denom;
}

static double gamma3d_objective_zpartial_zero(const double X[_q(static 3)],
                                              const double S[_q(static 3)],
                                              const double Xm[_q(static 3)],
                                              const double interp[_q(static 8)])
{
    const double deriv = fma(
        X[1],
        fma(X[0], interp[7], interp[6]),
        fma(X[0], interp[5], interp[4])
    );
    const double comp = fma(
        X[2],
        fma(X[1], interp[3], interp[2]),
        fma(X[1], interp[1], interp[0])
    );
    const double numer = S[2] * Xm[2] + comp * deriv;
    const double denom = S[2] + deriv * deriv;
    return -numer / denom;
}

static void gamma3d_clamp_unity(double *x)
{
    if (*x < 0.0) {
        *x = 0.0;
    } else if (*x > 1.0) {
        *x = 1.0;
    }
}

static void gamma3d_line_search_move(double X[_q(static 3)],
                                     const double S[_q(static 3)],
                                     const double Xm[_q(static 3)],
                                     const double interp[_q(static 3)])
{
    X[0] = gamma3d_objective_xpartial_zero(X, S, Xm, interp);
    gamma3d_clamp_unity(X);
    X[1] = gamma3d_objective_ypartial_zero(X, S, Xm, interp);
    gamma3d_clamp_unity(X + 1);
    X[2] = gamma3d_objective_zpartial_zero(X, S, Xm, interp);
    gamma3d_clamp_unity(X + 2);
}

static double gamma3d_line_search(const double S[_q(static 3)],
                                  const double Xm[_q(static 3)],
                                  const double interp[_q(static 4)])
{
    double X[3];
    double obj, last;
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = 0.0;
    obj = gamma3d_objective_function(X, S, Xm, interp);
    do {
        last = obj;
        gamma3d_line_search_move(X, S, Xm, interp);
        obj = gamma3d_objective_function(X, S, Xm, interp);
    } while ((last - obj) > LINE_SEARCH_THRESHOLD);
    return obj;
}

static int gamma3d_minimize_square(const struct gamma_analysis *gamma,
                                   const struct dose_distribution *calc,
                                   double *rmax,
                                   const long cpoint[_q(static 3)],
                                   const double r_m[_q(static 3)],
                                   double mdose)
{
    double rmin;
    double Xm[3];
    double interp[8];
    double S[3];
    gamma3d_vec_loadz(Xm, cpoint[0], cpoint[1], cpoint[2]);
    gamma3d_vec_sub(Xm, r_m);
    gamma3d_vec_copy(S, calc->px_spacing);
    gamma3d_vec_mul(S, S);
    if (!gamma3d_cube_in_radius(Xm, S, *rmax)) {
        return 0;
    }
    gamma3d_load_interpolant(gamma, calc, cpoint, interp, mdose);
    rmin = gamma3d_line_search(S, Xm, interp);
    if (rmin < *rmax) {
        *rmax = rmin;
    }
    return 1;
}

static int gamma3d_fan_horizontal(const struct gamma_analysis *gamma,
                                  const struct dose_distribution *calc,
                                  double *rmax,
                                  long cpoint[_q(static 3)],
                                  const double r_m[_q(static 3)],
                                  double mdose)
{
    int n = 0, aug;
    const long x_c = cpoint[0];
    long offset = 0;
    aug = gamma3d_minimize_square(gamma, calc, rmax, cpoint, r_m, mdose);
    if (!aug) {
        return 0;
    }
    do {
        n += aug;
        offset++;
        cpoint[0] = x_c + offset;
        aug = gamma3d_minimize_square(gamma, calc, rmax, cpoint, r_m, mdose);
        cpoint[0] = x_c - offset;
        aug += gamma3d_minimize_square(gamma, calc, rmax, cpoint, r_m, mdose);
    } while (aug);
    cpoint[0] = x_c;
    return n;
}

static int gamma3d_fan_vertical(const struct gamma_analysis *gamma,
                                 const struct dose_distribution *calc,
                                 double *rmax,
                                 long cpoint[_q(static 3)],
                                 const double r_m[_q(static 3)],
                                 double mdose)
{
    int n;
    const long y_c = cpoint[1];
    long offset = 0;
    if (!gamma3d_fan_horizontal(gamma, calc, rmax, cpoint, r_m, mdose)) {
        return 0;
    }
    do {
        offset++;
        cpoint[1] = y_c + offset;
        n = gamma3d_fan_horizontal(gamma, calc, rmax, cpoint, r_m, mdose);
        cpoint[1] = y_c - offset;
        n += gamma3d_fan_horizontal(gamma, calc, rmax, cpoint, r_m, mdose);
    } while (n);
    cpoint[1] = y_c;
    return 1;
}

static void gamma3d_fan_internal(const struct gamma_analysis *gamma,
                                 const struct dose_distribution *calc,
                                 double *rmax,
                                 long cpoint[_q(static 3)],
                                 const double r_m[_q(static 3)],
                                 double mdose)
{
    int n;
    const long z_c = cpoint[2];
    long offset = 0;
    gamma3d_fan_vertical(gamma, calc, rmax, cpoint, r_m, mdose);
    do {
        offset++;
        cpoint[2] = z_c + offset;
        n = gamma3d_fan_vertical(gamma, calc, rmax, cpoint, r_m, mdose);
        cpoint[2] = z_c - offset;
        n += gamma3d_fan_vertical(gamma, calc, rmax, cpoint, r_m, mdose);
    } while (n);
}

static double gamma3d_compute_point(const struct gamma_analysis *gamma,
                                    const struct dose_distribution *calc,
                                    const double r_m[_q(static 3)],
                                    double mdose)
{
    long origin[3];
    double rmax = HUGE_VAL;
    gamma3d_vec_floorz(r_m, origin);
    gamma3d_fan_internal(gamma, calc, &rmax, origin, r_m, mdose);
    return sqrt(rmax / (gamma->dta * gamma->dta));
}

static int gamma3d_fan_horizontal_pass_only(const struct gamma_analysis *gamma,
                                            const struct dose_distribution *calc,
                                            double *rmax,
                                            long cpoint[_q(static 3)],
                                            const double r_m[_q(static 3)],
                                            double mdose,
                                            jmp_buf env)
{
    const double rnorm = gamma->dta * gamma->dta;
    int n = 0, aug;
    const long x_c = cpoint[0];
    long offset = 0;
    aug = gamma3d_minimize_square(gamma, calc, rmax, cpoint, r_m, mdose);
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
        aug = gamma3d_minimize_square(gamma, calc, rmax, cpoint, r_m, mdose);
        if (*rmax < rnorm) {
            longjmp(env, 1);
        }
        cpoint[0] = x_c - offset;
        aug += gamma3d_minimize_square(gamma, calc, rmax, cpoint, r_m, mdose);
        if (*rmax < rnorm) {
            longjmp(env, 1);
        }
    } while (aug);
    cpoint[0] = x_c;
    return n;
}

static int gamma3d_fan_vertical_pass_only(const struct gamma_analysis *gamma,
                                          const struct dose_distribution *calc,
                                          double *rmax,
                                          long cpoint[_q(static 3)],
                                          const double r_m[_q(static 3)],
                                          double mdose,
                                          jmp_buf env)
{
    int n;
    const long y_c = cpoint[1];
    long offset = 0;
    if (!gamma3d_fan_horizontal_pass_only(gamma, calc, rmax, cpoint, r_m, mdose, env)) {
        return 0;
    }
    do {
        offset++;
        cpoint[1] = y_c + offset;
        n = gamma3d_fan_horizontal_pass_only(gamma, calc, rmax, cpoint, r_m, mdose, env);
        cpoint[1] = y_c - offset;
        n += gamma3d_fan_horizontal_pass_only(gamma, calc, rmax, cpoint, r_m, mdose, env);
    } while (n);
    cpoint[1] = y_c;
    return 1;
}

static void gamma3d_fan_internal_pass_only(const struct gamma_analysis *gamma,
                                           const struct dose_distribution *calc,
                                           double *rmax,
                                           long cpoint[_q(static 3)],
                                           const double r_m[_q(static 3)],
                                           double mdose,
                                           jmp_buf env)
{
    int n;
    const long z_c = cpoint[2];
    long offset = 0;
    gamma3d_fan_vertical_pass_only(gamma, calc, rmax, cpoint, r_m, mdose, env);
    do {
        offset++;
        cpoint[2] = z_c + offset;
        n = gamma3d_fan_vertical_pass_only(gamma, calc, rmax, cpoint, r_m, mdose, env);
        cpoint[2] = z_c - offset;
        n += gamma3d_fan_vertical_pass_only(gamma, calc, rmax, cpoint, r_m, mdose, env);
    } while (n);
}

static int gamma3d_point_passes(const struct gamma_analysis *gamma,
                                const struct dose_distribution *calc,
                                const double r_m[_q(static 3)],
                                double mdose)
{
    /* volatile */ long origin[3];
    /* volatile */ double rmax = HUGE_VAL;
    /* volatile */ jmp_buf env;
    gamma3d_vec_floorz(r_m, origin);
    if (setjmp(env)) {
        return 1;
    }
    gamma3d_fan_internal_pass_only(gamma, calc, &rmax, origin, r_m, mdose, env);
    return 0;
}

static void gamma3d_compute_dist(struct gamma_analysis *gamma,
                                 const struct dose_distribution *restrict meas,
                                 const struct dose_distribution *restrict calc,
                                 double *gptr)
{
    long i, j, k;
    const double *mptr = meas->data;
    double trns[3], shft[3], a[3];
    gamma3d_init_transform(trns, meas, calc);
    gamma3d_init_shift(shft, meas, calc);
    for (k = 0; k < meas->px_dim[2]; k++) {
        for (j = 0; j < meas->px_dim[1]; j++) {
            for (i = 0; i < meas->px_dim[0]; i++, mptr++, gptr++) {
                if (*mptr < gamma->threshold) {
                    *gptr = GAMMA_SIG_VAL;
                    continue;
                }
                gamma3d_vec_loadz(a, i, j, k);
                gamma3d_vec_fma(a, trns, shft);
                *gptr = gamma3d_compute_point(gamma, calc, a, *mptr);
            }
        }
    }
}

static double gamma3d_compute_pass_rate(struct gamma_analysis *gamma,
                                        const struct dose_distribution *restrict meas,
                                        const struct dose_distribution *restrict calc)
{
    long i, j, k;
    const double *mptr = meas->data;
    double trns[3], shft[3], a[3];
    long passed = 0, npts = 0;
    gamma3d_init_transform(trns, meas, calc);
    gamma3d_init_shift(shft, meas, calc);
    for (k = 0; k < meas->px_dim[2]; k++) {
        for (j = 0; j < meas->px_dim[1]; j++) {
            for (i = 0; i < meas->px_dim[0]; i++, mptr++) {
                if (*mptr < gamma->threshold) {
                    continue;
                }
                gamma3d_vec_loadz(a, i, j, k);
                gamma3d_vec_fma(a, trns, shft);
                if (gamma3d_point_passes(gamma, calc, a, *mptr)) {
                    passed++;
                }
                npts++;
            }
        }
    }
    if (!npts) {
        return GAMMA_SIG_VAL;
    } else {
        return STATIC_CAST(double, passed) / STATIC_CAST(double, npts);
    }
}

static double gamma3d_find_max(long N, const double dptr[_q(static N)])
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

static double gamma3d_pass_rate_from_dist(long N, const double gptr[_q(static N)])
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

double gamma3d_compute(struct gamma_analysis *gamma,
                       const struct dose_distribution *restrict meas,
                       const struct dose_distribution *restrict calc,
                       double *gdist)
{
    const double threshold = gamma->threshold;
    const double diff = gamma->diff;
    const double calcmax = gamma3d_find_max(calc->px_dim[0] * calc->px_dim[1] * calc->px_dim[2], calc->data);
    double pass_rate;

    if (!calcmax) {
        return GAMMA_SIG_VAL;
    }

    if (gamma->normalization_mode == GAMMA_NORM_GLOBAL) {
        gamma->diff *= calcmax;
    }
    gamma->threshold *= calcmax;

    if (gdist) {
        gamma3d_compute_dist(gamma, meas, calc, gdist);
        pass_rate = gamma3d_pass_rate_from_dist(meas->px_dim[0] * meas->px_dim[1] * meas->px_dim[2], gdist);
    } else {
        pass_rate = gamma3d_compute_pass_rate(gamma, meas, calc);
    }

    gamma->diff = diff;
    gamma->threshold = threshold;
    return pass_rate;
}
