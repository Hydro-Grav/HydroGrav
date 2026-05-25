/**
 * @file physics.hpp
 * @brief Contains basic relativistic physics utility functions.
 */
#ifndef INCLUDE_PHYSICS_HPP_H
#define INCLUDE_PHYSICS_HPP_H

#include <vector>
#include <string>

static constexpr double kB = 8.61733326e-14; // Boltzmann constant (GeV/K)
static constexpr double mP = 1.2209e19; // Planck mass (GeV)
static constexpr double H_0 = 67.8 / 3.086e19;  // Hubble constant in SI units
static constexpr double h = 0.678;  // Reduced Hubble constant

/**
 * @brief Computes the square of the Lorentz factor γ² for a given velocity.
 * 
 * @param v The velocity (as a fraction of the speed of light, 0 ≤ v < 1).
 * @return The square of the Lorentz factor γ².
 */
double gammaSq(double v);

/**
 * @brief Optical Metrology System (OMS) noise power spectral density for LISA.
 *
 * @param f Frequency in Hz.
 * @return OMS noise PSD in m²/Hz.
 */
double P_oms(double f);

/**
 * @brief Acceleration noise power spectral density for LISA.
 *
 * @param f Frequency in Hz.
 * @return Acceleration noise PSD in (m/s²)²/Hz.
 */
double P_acc(double f);

/**
 * @brief LISA sensitivity curve expressed as an energy density fraction.
 *
 * Returns the effective noise level \f$\Omega_{\rm GW} h^2\f$ of the LISA
 * detector at frequency @p f, combining OMS and acceleration noise contributions.
 *
 * @param f Frequency in Hz.
 * @return \f$\Omega_{\rm GW} h^2\f$ sensitivity at @p f.
 */
double get_LISA_omegahsq(double f);

/**
 * @brief Galactic binary confusion noise in gravitational-wave strain.
 *
 * @param f Frequency in Hz.
 * @return Strain noise spectral density (Hz⁻¹).
 */
double gb_S(double f);

/**
 * @brief Galactic binary confusion noise as an energy density fraction.
 *
 * @param f Frequency in Hz.
 * @return \f$\Omega_{\rm GW} h^2\f$ contribution from galactic binaries.
 */
double gb_omegahsq(double f);

/**
 * @brief Extragalactic binary confusion noise as an energy density fraction.
 *
 * @param f    Frequency in Hz.
 * @param flag Amplitude variant: @c "central" (default), @c "upper", or @c "lower".
 * @return \f$\Omega_{\rm GW} h^2\f$ contribution from extragalactic binaries.
 */
double eb_omegahsq(double f, const std::string& flag = "central");

/**
 * Calculate LISA Signal-to-Noise Ratio
 * 
 * @param freqVals Vector of frequency values in Hz
 * @param ampVals Vector of GW amplitude Omega_GW h^2 values
 * @param Tyear Observation time in years (default: 4)
 * @return SNR value
 */
double LISA_snr(const std::vector<double>& freqVals, 
                const std::vector<double>& ampVals,
                double Tyear = 4.0);

#endif // INCLUDE_PHYSICS_HPP_H
