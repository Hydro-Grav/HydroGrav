#ifndef INCLUDE_CONFIG_HPP_H
#define INCLUDE_CONFIG_HPP_H

namespace config {
    
    /* 
        Chi cut-off where we switch from Boost's Gauss-Legendre quadrature to 
        an in-house Filon quadrature to handle large chi. 
    */
    constexpr double chi_threshold = 1e1; // default = 1e1

    /* Polynomial order used in Filon quadrature */
    constexpr int filon_polynomial_order = 16; // default = 16

    /* Number of sample points for Gauss-Legendre quadrature */
    constexpr int fd_l_gauss_legendre_samples = 15; // default = 15

    /* 
        Parameters for boost::math::quadrature::gauss_kronrod 
        used in Ekin integrations
    */
    constexpr int Ekin_samples = 31; // default = 61, safe = 31
    constexpr int Ekin_max_refinements = 5; // default = 10, safe = 5
    constexpr double Ekin_tolerance = 1e-6; // default = 1e-12, safe = 1e-6

    /* Bounds for ApSQ(chi) spline fitting and integration */
    constexpr double chi_min = 1e-3; // default = 1e-3
    constexpr double chi_max = 1e4; // default = 1e4
    constexpr int chi_points = 5000; // default = 5000

    /* Bounds for pRs spline fitting and integration */
    constexpr double pRs_minimum = 1e-3; // default = 1e-3
    constexpr double pRs_maximum = 1e+3; // default = 1e+3
    constexpr int n_pRs = 500; // default = 500

    /* Bounds for Ekin integration (tau_nl) */
    constexpr double kRs_minimum = 1e-4; // default = 1e-4
    constexpr double kRs_maximum = 1e+5; // default = 1e+5
    constexpr int n_kRs = 5000; // default = 5000

    /*  
        Factor for extending kinetic spectrum spline bounds
        beyond min/max pRs values found in GWSpec
    */
    constexpr double kinetic_spectrum_spline_factor = 0.01; // default = 1%

    /* Number of points for zetaKin spline in GWSpec */
    constexpr int kinetic_spectrum_spline_points = 2 * n_pRs; // default = 2 * n_pRs

    /*
        Parameters for boost::math::quadrature::gauss_kronrod 
        used in z integral for GWSpec
    */
    constexpr int z_samples = 31; // default = 31, safe = 31
    constexpr int z_max_refinements = 5; // default = 6, safe = 5
    constexpr double z_tolerance = 1e-6; // default = 1e-8, safe = 1e-6

    /*
        Parameters for boost::math::quadrature::gauss_kronrod 
        used in pRs integral for GWSpec
    */
    constexpr int pRs_samples = 15; // default = 31, safe = 15
    constexpr int pRs_max_refinements = 5; // default = 6, safe = 5
    constexpr double pRs_tolerance = 1e-6; // default = 1e-8, safe = 1e-6

}

#endif // INCLUDE_CONFIG_HPP_H
