/** See source files for copyright information */
#pragma once

#ifndef GAMMA_ND_H
#define GAMMA_ND_H

#if __cplusplus
extern "C" {
#endif

/** This value is assigned to points that fall below the low-dose
 *  threshold. It is also returned as the pass rate if there was an error
 *  during calculation
 */
#define GAMMA_SIG_VAL -1.0


/** Only the first two elements of each member array are accessed by
 *  gamma2d
 */
struct dose_distribution {
    double px_spacing[3];
    double top_left[3];
    long px_dim[3];
    double *data;
};


struct gamma_analysis {
    double diff;        /* Percentage (<= 1.0, so 3% = 0.03), unless absolute 
                           normalization is selected, in which case it is an
                           absolute dose value */
    double dta;         /* Same units as the input pixel spacings */
    double threshold;   /* Percentage (<= 1.0) */
    enum {
        GAMMA_NORM_LOCAL,   /* Measured dose is used to norm dose diff */
        GAMMA_NORM_GLOBAL,  /* Global maximum dose is used to norm dose diff */ 
        GAMMA_NORM_ABSOLUTE /* The value of diff is used to norm dose diff */
    } normalization_mode;
};


/** Computes the 2D gamma index pass rate for dose distributions @p meas and 
 *  @p calc, optionally computing the gamma distribution.
 * 
 *  \param gamma
 *      Gamma analysis parameter struct
 *  \param meas
 *      Measured dose distribution
 *  \param calc
 *      Calculated/TPS dose distribution
 *  \param gdist
 *      Optional pointer to preallocated buffer intended to store the gamma 
 *      distribution. This buffer's dimensions and physical properties must 
 *      be identical to those of @p meas. If this argument is NULL, only the 
 *      pass rate will be computed. 
 *      
 *  \returns Gamma pass rate, or GAMMA_SIG_VAL if an error occurred. Signal 
 *      value is returned if either the calculated distribution is zero 
 *      everywhere, or if no points in the measured distribution were above 
 *      the low-dose threshold.
 */
double gamma2d_compute(struct gamma_analysis *gamma,
                       const struct dose_distribution *meas,
                       const struct dose_distribution *calc,
                       double *gdist);


/** See gamma2d, call profile is identical */
double gamma3d_compute(struct gamma_analysis *gamma,
                       const struct dose_distribution *meas,
                       const struct dose_distribution *calc,
                       double *gdist);


#ifdef __cplusplus
}
#endif

#endif /* GAMMA_ND_H */
