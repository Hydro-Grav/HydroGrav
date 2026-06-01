#ifndef INCLUDE_SNR_HPP_H
#define INCLUDE_SNR_HPP_H

#include <cmath>
#include <vector>
#include <string>

static constexpr double H_0 = 67.8 / 3.086e19;  // Hubble constant in SI units
static constexpr double h   = 0.678;            // Reduced Hubble constant

// Base detector class
class Detector {
public:
    virtual ~Detector() = default;

    virtual std::string name() const = 0;

    // Effective number of independent detector pairs used in the SNR integral.
    // Override in derived classes where ndet != 1 (e.g. DECIGO).
    virtual double ndet() const { return 1.0; }

    // Noise sensitivity curve in Omega_GW h^2 units.
    virtual double omega_noise(double f) const = 0;

protected:
    // Shared conversion from strain noise Sn to Omega_GW h^2.
    double sn_to_omegahsq(double f, double Sn) const {
        return (2.0 * M_PI * M_PI / (3.0 * H_0 * H_0)) * std::pow(f, 3) * Sn * h * h;
    }
};
    
// LISA
class LISA : public Detector {
public:
    std::string name() const override { return "LISA"; }

    double omega_noise(double f) const override {
        return sn_to_omegahsq(f, Sn(f));
    }

private:

    // Galactic binary confusion noise
    // double gb_omegahsq(double /*f*/) const { return 0.0; }

    // Extragalactic binary confusion noise
    // double eb_omegahsq(double /*f*/, const std::string& /*flag*/ = "central") const { return 0.0; }

    double P_oms(double f) const {
        return (1.5e-11) * (1.5e-11) * (1.0 + std::pow(2e-3 / f, 4));
    }

    double P_acc(double f) const {
        double term1 = (3e-15) * (3e-15) / std::pow(2.0 * M_PI * f, 4);
        double term2 = 1.0 + std::pow(0.4e-3 / f, 2);
        double term3 = 1.0 + std::pow(f / 8e-3, 4);
        return term1 * term2 * term3;
    }

    double Sn(double f) const {
        double L  = 2.5e9;                         // Arm length in metres
        double fs = 2.998e8 / (2.0 * M_PI * L);    // Transfer frequency

        double Sn = (10.0 / (3.0 * L * L))
                * (P_oms(f) + 4.0 * P_acc(f))
                * (1.0 + 0.54 * std::pow(f / fs, 2));
        return Sn;
    }
};

// BBO - AI SLOP DOUBLE CHECK THIS ONE
class BBO : public Detector {
public:
    std::string name() const override { return "BBO"; }

    double omega_noise(double f) const override {
        return sn_to_omegahsq(f, Sn(f));
    }

private:
    // Strain noise spectral density (Cutler & Harms 2006 approximation)
    double Sn(double f) const {
        double S_acc  = 9.0e-54 / std::pow(f, 4);
        double S_sn   = 3.6e-37 * std::pow(f, 2);
        double S_omni = 1.8e-44;
        return S_acc + S_sn + S_omni;
    }
};

// DECIGO - AI SLOP DOUBLE CHECK THIS ONE
class DECIGO : public Detector {
public:
    std::string name() const override { return "DECIGO"; }
    double ndet() const override { return 3.0; }  // 3 correlated detector pairs

    double omega_noise(double f) const override {
        return sn_to_omegahsq(f, Sn(f));
    }

private:
    // Strain noise spectral density (Seto et al. 2001 approximation)
    double Sn(double f) const {
        double f_p   = 7.36;
        double S_rp  = 1.05e-46 / (1.0 + std::pow(f / f_p, 2));  // Radiation pressure
        double S_sn  = 6.53e-51 * std::pow(f, 2);                // Shot noise
        double S_acc = 3.15e-60 / std::pow(f, 4);                // Acceleration
        return S_rp + S_sn + S_acc;
    }
};

// Einstein Telescope - AI SLOP DOUBLE CHECK THIS ONE
class EinsteinTelescope : public Detector {
public:
    std::string name() const override { return "ET"; }

    double omega_noise(double f) const override {
        return sn_to_omegahsq(f, Sn(f));
    }

private:
    // ET-D sensitivity curve approximation (Hild et al. 2011)
    double Sn(double f) const {
        if (f < 1.0) return 1.0;

        double a1   = 2.39e-27, a2 = -0.142, a3 = -3.0;
        double b1   = 0.349,    b2 = 1.76,   b3 = -2.0;
        double c1   = 1.76e-48, c2 = 5.0;

        double S_low  = a1 * std::pow(f, a2) + a3;
        double S_mid  = b1 * std::pow(f, b2) + b3;
        double S_high = c1 * std::pow(f, c2);

        return S_low * S_low + S_mid * S_mid + S_high;
    }
};

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
 * Convenience function to calculate SNR for BBO only.
 *
 * @param freqVals  Frequency values in Hz
 * @param ampVals   Signal amplitude in Omega_GW h^2
 * @param Tyear     Observation time in years
 * @return SNR value for BBO
 */
double BBO_snr(
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
 * Convenience function to calculate SNR for Einstein Telescope only.
 *
 * @param freqVals  Frequency values in Hz
 * @param ampVals   Signal amplitude in Omega_GW h^2
 * @param Tyear     Observation time in years
 * @return SNR value for Einstein Telescope
 */
double ET_snr(
    const std::vector<double>& freqVals,
    const std::vector<double>& ampVals,
    double Tyear = 4.0);

#endif // INCLUDE_SNR_HPP_H
