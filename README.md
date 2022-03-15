# Gamma Analysis

ANSI C implementation of the gamma analysis algorithm commonly used in medical physics. Comes in 2- and 3-dimensional versions, and optionally computes the gamma distribution.

Running time is _Ω_(_M_<sup>2</sup>_γ<sup>k</sup>_) to calculate the distribution, or _Θ_(_M_<sup>2</sup>_γ<sup>k</sup>_) to calculate the pass rate, where:
- _M_ is the number of measured dose points that are above the low-dose threshold.
- _γ_ is the average pointwise gamma value.
- _k_ is the dimensionality of the problem.

The algorithm splits the calculated dose into lattice cells, each bounded by 2<sup>_k_</sup> dose points, and interpolates the TPS dose within. The _Γ_ function for each lattice cell is then minimized by simple [coordinate descent](https://en.wikipedia.org/wiki/Coordinate_descent), and the minimum of all lattice cells is taken as the result.

Note that despite the poor asymptotic bounds to calculate the distribution (_Ω_), this runs quite fast in practice if the doses mostly match and the average pointwise gamma is low.
