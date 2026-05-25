/**
 * @file constants.hpp
 * @brief Physical and mathematical constants used throughout HydroGrav.
 */
#ifndef INCLUDE_CONSTANTS_HPP_H
#define INCLUDE_CONSTANTS_HPP_H

/// Speed of light in vacuum (m/s)
constexpr double c = 299792458.0;
/// Newton's gravitational constant (N·m²/kg²)
constexpr double G = 6.67430e-11;
/// Reduced Planck constant (J·s)
constexpr double hbar = 1.05e-34;

/// Bag model speed of sound in fluid
constexpr double csq_bag = 1/3;

/// Euler–Mascheroni constant
constexpr double gamma_euler = 0.5772156649015328606;

#endif // INCLUDE_CONSTANTS_HPP_H
