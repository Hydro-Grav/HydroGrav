#ifndef INCLUDE_SNR_HPP_H
#define INCLUDE_SNR_HPP_H

#include <vector>
#include <string>

static constexpr double H_0 = 67.8 / 3.086e19;  // Hubble constant in SI units
static constexpr double h = 0.678;  // Reduced Hubble constant

// LISA SNR calculation

double P_oms(double f);

// Acceleration noise power spectral density
double P_acc(double f);

// LISA sensitivity curve in Omega_GW h^2
double get_LISA_omegahsq(double f);

// Galactic binary confusion noise in strain
double gb_S(double f);

// Galactic binary confusion noise in Omega_GW h^2
double gb_omegahsq(double f);

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

#endif // INCLUDE_SNR_HPP_H