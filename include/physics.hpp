// physics.hpp
#ifndef INCLUDE_PHYSICS_HPP_H
#define INCLUDE_PHYSICS_HPP_H

/**
 * @file physics.hpp
 * @brief Contains basic relativistic physics utility functions.
 */

static constexpr double kB = 8.61733326e-14; // Boltzmann constant (GeV/K)

/**
 * @brief Computes the Lorentz factor γ for a given velocity.
 * 
 * @param v The velocity (as a fraction of the speed of light, 0 ≤ v < 1).
 * @return The Lorentz factor γ.
 */
// double gamma(double v);

/**
 * @brief Computes the square of the Lorentz factor γ² for a given velocity.
 * 
 * @param v The velocity (as a fraction of the speed of light, 0 ≤ v < 1).
 * @return The square of the Lorentz factor γ².
 */
double gammaSq(double v);

#endif // INCLUDE_PHYSICS_HPP_H
