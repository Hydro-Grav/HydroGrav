#ifndef INCLUDE_SNR_HPP_H
#define INCLUDE_SNR_HPP_H

#include <cmath>
#include <vector>
#include <string>

#include "detectors.hpp"

struct SNRResult {
    std::string detector_name;
    double snr;
};

/**
 * Calculate SNR for a single detector.
 *
 * @param freqVals  Frequency values in Hz
 * @param ampVals   Signal amplitude in Omega_GW h^2
 * @param detector  Detector instance providing the noise curve
 * @param Tyear     Observation time in years
 * @return SNR value
 */
double calculate_snr(
    const std::vector<double>& freqVals,
    const std::vector<double>& ampVals,
    const Detector& detector,
    double Tyear = 4.0);

/**
 * Calculate SNR for multiple detectors over the same signal.
 *
 * @param freqVals  Frequency values in Hz
 * @param ampVals   Signal amplitude in Omega_GW h^2
 * @param detectors List of detector instances
 * @param Tyear     Observation time in years
 * @return Vector of SNRResult, one per detector, in input order
 */
std::vector<SNRResult> calculate_all_snrs(
    const std::vector<double>& freqVals,
    const std::vector<double>& ampVals,
    const std::vector<const Detector*>& detectors,
    double Tyear = 4.0);

/**
 * Convenience function to calculate SNR for all standard detectors (LISA, BBO, DECIGO, ET).
 * 
 * @param freqVals  Frequency values in Hz
 * @param ampVals   Signal amplitude in Omega_GW h^2
 * @param Tyear     Observation time in years
 * @return Vector of SNRResult for LISA, BBO, DECIGO, and ET, in that order
 */
std::vector<SNRResult> get_SNR(
    const std::vector<double>& freqVals,
    const std::vector<double>& ampVals,
    double Tyear = 4.0);

/**
 * Convenience function to calculate SNR for LISA only.
 *
 * @param freqVals  Frequency values in Hz
 * @param ampVals   Signal amplitude in Omega_GW h^2
 * @param Tyear     Observation time in years
 * @return SNR value for LISA
 */
double LISA_snr(
    const std::vector<double>& freqVals,
    const std::vector<double>& ampVals,
    double Tyear = 4.0);

/**
 * Convenience function to calculate SNR for DECIGO only.
 *
 * @param freqVals  Frequency values in Hz
 * @param ampVals   Signal amplitude in Omega_GW h^2
 * @param Tyear     Observation time in years
 * @return SNR value for DECIGO
 */
double DECIGO_snr(
    const std::vector<double>& freqVals,
    const std::vector<double>& ampVals,
    double Tyear = 4.0);

/**
 * Convenience function to calculate SNR for BBO only.
 *
 * @param freqVals  Frequency values in Hz
 * @param ampVals   Signal amplitude in Omega_GW h^2
 * @param Tyear     Observation time in years
 * @return SNR value for BBO
 */
// double BBO_snr(
//     const std::vector<double>& freqVals,
//     const std::vector<double>& ampVals,
//     double Tyear = 4.0);

#endif // INCLUDE_SNR_HPP_H
