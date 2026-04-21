// physics.cpp
#include <cmath>
#include <vector>
#include <stdexcept>
#include <string>

#include "physics.hpp"
#include "maths_ops.hpp"

double gammaSq(double v) {
    return 1.0 / (1.0 - v*v);
}

// LISA SNR calculation - move this somewhere else?
// Optical Metrology System noise power spectral density
double P_oms(double f) {
    return 3.6e-41;
}

// Acceleration noise power spectral density
double P_acc(double f) {
    double term1 = 1.44e-48 / std::pow(2.0 * M_PI * f, 4);
    double term2 = 1.0 + std::pow(0.4e-3 / f, 2);
    return term1 * term2;
}

// LISA sensitivity curve in Omega_GW h^2
double get_LISA_omegahsq(double f) {
    double Poms = P_oms(f);
    double Pacc = P_acc(f);
    double L = 2.5e9;  // LISA arm length in meters (2.5 million km)
    double fs = 2.998e8 / (2.0 * M_PI * L);  // Transfer frequency
    
    // Strain noise spectral density
    double Sn = (40.0 / 3.0) * (Poms + 4.0 * Pacc) * (1.0 + std::pow(f / (4.0 * fs / 3.0), 2));
    
    // Convert to Omega_GW h^2
    return 4.0 * M_PI * M_PI / (3.0 * H_0 * H_0) * std::pow(f, 3) * Sn * h * h;
}

// Galactic binary confusion noise in strain
double gb_S(double f) {
    // Convert to millihertz
    double f_mHz = f * 1e3;
    
    double A = 8e-38;
    double f_ref = 1000.0;
    double a = 0.138;
    double b = -0.221;
    double c = 0.521;
    double d = 1.680;
    double fk = 1.13;
    
    double exponent = -(f_mHz / f_ref) * a - b * f_mHz * std::sin(c * f_mHz);
    double tanh_term = 1.0 + std::tanh(d * (fk - f_mHz));
    
    return A * std::pow(f_mHz, 7.0/3.0) * std::exp(exponent) * tanh_term;
}

// Galactic binary confusion noise in Omega_GW h^2
double gb_omegahsq(double f) {
    double Sn = gb_S(f);
    return 4.0 * M_PI * M_PI / (3.0 * H_0 * H_0) * std::pow(f, 3) * Sn * h * h;
}

/**
 * Extragalactic binary confusion noise in Omega_GW h^2
 * @param f Frequency in Hz
 * @param flag "central", "upper", or "lower" for different amplitude estimates
 */
double eb_omegahsq(double f, const std::string& flag) {
    double amp = 8.9e-10;
    
    if (flag == "upper") {
        amp += 12.6e-10;
    } else if (flag == "lower") {
        amp -= 5.6e-10;
    }
    // else flag == "central", use base amplitude
    
    double f_ref = 25.0;
    return amp * std::pow(f / f_ref, 2.0/3.0) * h * h;
}

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
                double Tyear) {
    
    if (freqVals.size() != ampVals.size()) {
        throw std::invalid_argument("freqVals and ampVals must have the same size");
    }
    
    if (freqVals.empty()) {
        throw std::invalid_argument("Input vectors cannot be empty");
    }
    
    size_t n = freqVals.size();
    std::vector<double> integrand(n);
    
    // Compute integrand: (Omega_signal / Omega_noise)^2
    for (size_t i = 0; i < n; ++i) {
        double f = freqVals[i];
        double omegat = ampVals[i];  // Signal
        
        // Total noise
        double omegan = get_LISA_omegahsq(f);
        
        // Add astrophysical foregrounds at low frequencies
        if (f < 1e-2) {
            omegan += gb_omegahsq(f);
            omegan += eb_omegahsq(f, "central");
        }
        
        // Avoid division by zero
        if (omegan <= 0.0) {
            throw std::runtime_error("Noise PSD is non-positive at f = " + std::to_string(f));
        }
        
        integrand[i] = std::pow(omegat / omegan, 2);
    }
    
    // Integrate over frequency
    double result = simpson_integrate(freqVals, integrand);
    
    // Convert observation time to seconds
    double Tsec = Tyear * 365.25 * 24.0 * 60.0 * 60.0;
    
    // SNR formula
    return std::sqrt(2.0 * Tsec * result);
}