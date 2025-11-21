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

std::function<double(double)> lifetime_distribution_function(const std::string& nuc_type);

static void build_nodes_and_samples(const FluidProfile &prof, std::vector<double>& x, std::vector<double>& f_sin, std::vector<double>& f_cos, std::vector<double>& l_sin);

void create_fluid_integrand_splines(const FluidProfile& prof, 
    alglib::spline1dinterpolant& f_sin_spline, 
    alglib::spline1dinterpolant& f_cos_spline, 
    alglib::spline1dinterpolant& l_sin_spline);

static std::vector<double> trapezoid_weights(const std::vector<double>& x);

std::pair<std::vector<double>, std::vector<double>> fluid_profile_integrals(const std::vector<double>& chi_vals, const FluidProfile& prof);

std::vector<double> Ap_sq(const std::vector<double>& chi_vals, const FluidProfile& prof);

} // namespace Hydrodynamics

#endif // INCLUDE_HYDRODYNAMICS_HPP_H
