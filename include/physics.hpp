// physics.hpp
#ifndef INCLUDE_PHYSICS_HPP_H
#define INCLUDE_PHYSICS_HPP_H

#include <vector>
#include <string>

/**
 * @file physics.hpp
 * @brief Contains basic relativistic physics utility functions.
 */

static constexpr double kB = 8.61733326e-14; // Boltzmann constant (GeV/K)
static constexpr double mP = 1.2209e19; // Planck mass (GeV)

/**
 * @brief Computes the square of the Lorentz factor γ² for a given velocity.
 * 
 * @param v The velocity (as a fraction of the speed of light, 0 ≤ v < 1).
 * @return The square of the Lorentz factor γ².
 */
double gammaSq(double v);

#endif // INCLUDE_PHYSICS_HPP_H
