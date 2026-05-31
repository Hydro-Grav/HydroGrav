#include <cmath>
#include <vector>
#include <stdexcept>
#include <string>

#include "snr.hpp"
#include "maths_ops.hpp"

// LISA SNR calculation - move this somewhere else?

// Optical Metrology System noise power spectral density
double P_oms(double f) {
    return (1.5e-11) * (1.5e-11) * (1.0 + std::pow(2e-3 / f, 4));
}

// Acceleration noise power spectral density
double P_acc(double f) {
    double term1 = (3e-15) * (3e-15) / std::pow(2.0 * M_PI * f, 4);
    double term2 = 1.0 + std::pow(0.4e-3 / f, 2);
    double term3 = 1.0 + std::pow(f / 8e-3, 4);
    return term1 * term2 * term3;
}

// LISA sensitivity curve in Omega_GW h^2
double get_LISA_omegahsq(double f) {
    double Poms = P_oms(f);
    double Pacc = P_acc(f);
    double L = 2.5e9;  // LISA arm length in meters (2.5 million km)
    double fs = 2.998e8 / (2.0 * M_PI * L);  // Transfer frequency
    
    // Strain noise spectral density
    double Sn = (10.0 / (3.0 * L * L)) * (Poms + 4.0 * Pacc) * (1.0 + 0.54 * std::pow(f / fs, 2));
    
    // Convert to Omega_GW h^2
    return 2.0 * M_PI * M_PI / (3.0 * H_0 * H_0) * std::pow(f, 3) * Sn * h * h;
}

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
        double omegan = get_LISA_omegahsq(f); // Noise spectrum
        
        if (omegan <= 0.0) {
            throw std::runtime_error("Noise PSD is non-positive at f = " + std::to_string(f));
        }
        
        integrand[i] = std::pow(omegat / omegan, 2);
    }
    
    // Integrate over frequency
    double result = simpson_integrate(freqVals, integrand);
    
    // Convert observation time to seconds
    double Tsec = Tyear * 365.25 * 24.0 * 60.0 * 60.0;

    // effective number of detectors (1 for LISA - autocorrelation measurement)
    double ndet = 1.0;
    
    return std::sqrt(ndet * Tsec * result);
}