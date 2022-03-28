# Gamma Analysis

<p align='center'>
    <image src='example.png'>
</p>

ANSI C implementation of the gamma analysis algorithm commonly used in medical physics. Comes in 2- and 3-dimensional versions, and optionally computes the gamma distribution.

Running time is _Ω_(_M_<sup>2</sup>_γ<sup>k</sup>_) to calculate the distribution, or _Θ_(_M_<sup>2</sup>_γ<sup>k</sup>_) to calculate the pass rate, where
- _M_ is the number of measured dose points that are above the low-dose threshold.
- _γ_ is the average pointwise gamma value.
- _k_ is the dimensionality of the problem.

The algorithm splits the calculated dose into lattice cells, each bounded by 2<sup>_k_</sup> dose points, and linearly interpolates the TPS dose within. The _Γ_ function for each lattice cell is then minimized by simple [coordinate descent](https://en.wikipedia.org/wiki/Coordinate_descent), and the minimum of all lattice cells is taken as the result. 

Note that despite the poor asymptotic bounds on the time to calculate the distribution (_Ω_), this runs quite fast in practice if the doses mostly match and the average pointwise gamma is low.

This implementation makes no allocations, and uses a fixed amount of stack space. The source code is written entirely in ANSI C, with feature test macros selectively enabling features from C99 and C11.


## Usage

The header defines two structs, first of which is the dose distribution:

### Dose distribution

```C
struct dose_distribution {
    double px_spacing[3];
    double top_left[3];
    long px_dim[3];
    double *data;
};
```
where
- ``px_spacing`` are the spacings between adjacent pixels in the image, in physical coordinates.
- ``top_left`` are the coordinates to the 0-th pixel in the image, in the same units as the pixel spacing.
- ``px_dim`` are the dimensions of the image, in pixels.
- ``data`` is a pointer to the raw dose buffer. It is not modified by the algorithm.

Each member array contains space for three elements, but only the first two are accessed by the 2D algorithm. The dose data is completely owned by the caller, and must be stored linearly in row-major order. Iliffe vectors are unsupported.

The other struct contains the gamma parameters:

### Gamma parameters

```C
struct gamma_analysis {
    double diff;
    double dta;
    double threshold;
    enum {
        GAMMA_NORM_LOCAL,
        GAMMA_NORM_GLOBAL,
        GAMMA_NORM_ABSOLUTE
    } normalization_mode;
};
```
where
- ``diff`` is the dose difference criterion. This will most likely be a percentage, but it may also be an absolute dose value (see ``normalization_mode`` below).
- ``dta`` is the distance-to-agreement criterion. This must be in the same physical units as the input doses' physical dimensions.
- ``threshold`` is the low-dose threshold, strictly a percentage.
- ``normalization_mode`` can be one of the following:
    * ``GAMMA_NORM_LOCAL`` normalizes the dose difference to the pointwise measured dose.
    * ``GAMMA_NORM_GLOBAL`` normalizes to the maximum TPS dose.
    * ``GAMMA_NORM_ABSOLUTE`` normalizes to an absolute dose value supplied by the caller (as ``diff``).

Please note that <b>percentages are unitized (e.g. 75% = 0.75)</b>.

After setting up the relevant structs, a call to one of two functions can be made:

### Functions

```C
double gamma2d_compute(struct gamma_analysis *gamma,
                       const struct dose_distribution *meas,
                       const struct dose_distribution *calc,
                       double *gdist);

double gamma3d_compute(struct gamma_analysis *gamma,
                       const struct dose_distribution *meas,
                       const struct dose_distribution *calc,
                       double *gdist);
```
The return value is the gamma pass rate (%GP). The ``gamma_analysis`` struct is non-``const``, but is unmodified upon return.

If ``gdist`` is not ``NULL``, then the program will write the gamma distribution to the location it points to. This region _must_ be at least as large as the buffer containing the measured dose.

The header #defines the quantity ``GAMMA_SIG_VAL``. This is a signal value that is either
- returned by the compute function in case of error.
- assigned to points in the gamma distribution whose measured dose was below the low-dose threshold.


## Limitations
- Coordinate descent is <b>not guaranteed to converge to the absolute minimum</b>. [It is very unlikely to converge to a saddle point](https://doi.org/10.1007/s10107-019-01374-3), though.

- There is currently no provision for <b>rotating</b> a dose. Shifting is possible by changing the origin coordinates.

- If a point fails, the function will not report whether it failed high or low. This is reasonably easy for the caller to implement, but I may add this functionality in the future.
