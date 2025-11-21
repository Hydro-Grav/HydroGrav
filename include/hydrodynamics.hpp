// hydrodynamics.hpp
#ifndef INCLUDE_HYDRODYNAMICS_HPP_H
#define INCLUDE_HYDRODYNAMICS_HPP_H


#include <array>
#include <complex>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include "interpolation.h"

#include "profile.hpp"

namespace Hydrodynamics {

std::function<double(double)> lifetime_distribution_function(const std::string& nuc_type);

std::pair<std::vector<double>, std::vector<double>> fluid_profile_integrals(const std::vector<double>& chi_vals, const FluidProfile& prof);

std::vector<double> Ap_sq(const std::vector<double>& chi_vals, const FluidProfile& prof);

} // namespace Hydrodynamics

#endif // INCLUDE_HYDRODYNAMICS_HPP_H
