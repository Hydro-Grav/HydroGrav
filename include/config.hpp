/**
 * @file config.hpp
 * @brief Global numerical configuration parameters used in the sound shell model implementation.
 *
 * All tunable constants that control quadrature accuracy, spline grid
 * bounds and integration limits are collected here.  Modify these values
 * to trade accuracy for speed (or vice-versa); the defaults marked below
 * are the values used in the accompanying publication.
 *
 * @note Do **not** change these values unless you understand the effect on
 *       numerical accuracy.  The "safe" alternatives noted in the comments
 *       reduce runtime at the cost of some precision.
 */
#ifndef INCLUDE_CONFIG_HPP_H
#define INCLUDE_CONFIG_HPP_H

/**
 * @namespace config
 * @brief Compile-time numerical configuration constants.
 */
namespace config {
    /**
     * @brief Tolerance for residual functions in FluidProfile class.
     *
     * For EoS with p(T), e(T) ~ O(1e+8), need residuals to converge
     * to at least O(1e+10). Not used for bag/mu-nu model.
     */
    const static double minimiser_tol = 1e-9;
    
    /**
     * @brief Momentum threshold χ above which Filon quadrature replaces
     *        Gauss–Legendre for the fluid-profile integrals.
     *
     * For χ > chi_threshold the integrand oscillates rapidly and
     * Filon quadrature is more efficient; below this value the
     * standard Boost Gauss–Legendre rule is used.
     */
    constexpr double chi_threshold = 1e1; // default = 1e1

    /// Polynomial order used in the Filon quadrature for oscillatory integrals.
    constexpr int filon_polynomial_order = 16; // default = 16

    /// Number of sample points per subinterval for Gauss–Legendre quadrature
    /// used in the fluid-profile integrals \f$f'(\chi)\f$ and \f$l(\chi)\f$.
    constexpr int fd_l_gauss_legendre_samples = 15; // default = 15

    /**
     * @name Gauss–Kronrod parameters for kinetic energy (Ekin) integration
     * @{
     */
    /// Number of Gauss–Kronrod points for Ekin integration.
    constexpr int Ekin_samples = 31; // default = 61, safe = 31
    /// Maximum number of adaptive refinements for Ekin integration.
    constexpr int Ekin_max_refinements = 5; // default = 10, safe = 5
    /// Absolute tolerance for Ekin integration convergence.
    constexpr double Ekin_tolerance = 1e-6; // default = 1e-12, safe = 1e-6
    /** @} */

    /**
     * @name Spline grid for |A+|^2(chi)
     * @{
     */
    /// Lower bound of the χ grid used to fit the \f$|A_+|^2(\chi)\f$ spline.
    constexpr double chi_min = 1e-3; // default = 1e-3
    /// Upper bound of the χ grid used to fit the \f$|A_+|^2(\chi)\f$ spline.
    constexpr double chi_max = 1e4; // default = 1e4
    /// Number of sample points on the χ grid.
    constexpr int chi_points = 5000; // default = 5000
    /** @} */

    /**
     * @name Spline grid for the momentum variable pRs
     * @{
     */
    /// Minimum value of p·R_* used in the kinetic-spectrum spline and integration.
    constexpr double pRs_minimum = 1e-3; // default = 1e-3
    /// Maximum value of p·R_* used in the kinetic-spectrum spline and integration.
    constexpr double pRs_maximum = 1e+3; // default = 1e+3
    /// Number of points in the p·R_* grid.
    constexpr int n_pRs = 500; // default = 500
    /** @} */

    /**
     * @name Momentum grid for the non-linear timescale (tau_nl) calculation
     * @{
     */
    /// Minimum k·R_* value for the Ekin integration used in tau_nl.
    constexpr double kRs_minimum = 1e-4; // default = 1e-4
    /// Maximum k·R_* value for the Ekin integration used in tau_nl.
    constexpr double kRs_maximum = 1e+5; // default = 1e+5
    /// Number of k·R_* grid points for tau_nl.
    constexpr int n_kRs = 5000; // default = 5000
    /** @} */

    /**
     * @brief Fractional extension factor for the zetaKin spline range.
     *
     * The spline bounds are extended by this fraction beyond the min/max
     * p·R_* values encountered during GWSpec to avoid edge artefacts.
     */
    constexpr double kinetic_spectrum_spline_factor = 0.01; // default = 1%

    /// Number of evaluation points for the zetaKin spline constructed inside GWSpec.
    constexpr int kinetic_spectrum_spline_points = 2 * n_pRs; // default = 2 * n_pRs

    /**
     * @name Gauss–Kronrod parameters for the z integral in GWSpec
     * @{
     */
    /// Number of Gauss–Kronrod points for the z (angle) integral in GWSpec.
    constexpr int z_samples = 31; // default = 31, safe = 31
    /// Maximum adaptive refinements for the z integral.
    constexpr int z_max_refinements = 5; // default = 6, safe = 5
    /// Convergence tolerance for the z integral.
    constexpr double z_tolerance = 1e-6; // default = 1e-8, safe = 1e-6
    /** @} */

    /**
     * @name Gauss–Kronrod parameters for the pRs integral in GWSpec
     * @{
     */
    /// Number of Gauss–Kronrod points for the p·R_* integral in GWSpec.
    constexpr int pRs_samples = 15; // default = 31, safe = 15
    /// Maximum adaptive refinements for the p·R_* integral.
    constexpr int pRs_max_refinements = 5; // default = 6, safe = 5
    /// Convergence tolerance for the p·R_* integral.
    constexpr double pRs_tolerance = 1e-6; // default = 1e-8, safe = 1e-6
    /** @} */

}

#endif // INCLUDE_CONFIG_HPP_H
