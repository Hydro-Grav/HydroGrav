// hydrodynamics.hpp
#ifndef INCLUDE_HYDRODYNAMICS_HPP_H
#define INCLUDE_HYDRODYNAMICS_HPP_H


#include <array>
#include <complex>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include <gsl/gsl_integration.h>

#include "interpolation.h"

#include "profile.hpp"

namespace Hydrodynamics {

/**
 * @brief Returns a function that computes the lifetime distribution for the specified nucleation type.
 *
 * @param nuc_type Nucleation type ("exp" for exponential, "sim" for simultaneous)
 * 
 * @return Function that takes a double Ttilde and returns the lifetime distribution value.
 */
std::function<double(double)> lifetime_distribution_function(const std::string& nuc_type);

/**
 * @brief Computes the integrated profile functions f'(χ) and l(χ) defined in Eq. (30), (31) of Pol, Procacci, Caprini (2024).
 *
 * @param chi χ=k*T_n (k = momentum, T_n = lifetime of n'th bubble)
 * @param prof FluidProfile object
 * 
 * @return Pair of vector of the integrated profile functions f' and l.
 */
std::pair<std::vector<double>, std::vector<double>> fluid_profile_integrals(const std::vector<double>& chi_vals, const FluidProfile& prof);

void create_fluid_integrand_splines(const FluidProfile& prof, 
    alglib::spline1dinterpolant& f_sin_spline, 
    alglib::spline1dinterpolant& f_cos_spline, 
    alglib::spline1dinterpolant& l_sin_spline);

// double integrand_qawo(double xi, void *params);

// const double compute_gsl_QAWO(const double& chi, const alglib::spline1dinterpolant& v_spline, const gsl_integration_qawo_enum type);

/**
 * @brief Computes |A₊(χ)|², defined in Eq. (29) of Pol, Procacci, Caprini (2024).
 *
 * @param chi χ=k*T_n (k = momentum, T_n = lifetime of n'th bubble)
 * @param prof Fluid profile object.
 * 
 * @return Vector of squared modulus |A₊(χ)|² values
 */
std::vector<double> Ap_sq(const std::vector<double>& chi_vals, const FluidProfile& prof);

} // namespace Hydrodynamics

#endif // INCLUDE_HYDRODYNAMICS_HPP_H
