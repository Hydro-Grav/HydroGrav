#include <cmath>
#include <vector>
#include <stdexcept>
#include <string>

#include "snr.hpp"
#include "maths.hpp"

double calculate_snr
(
    const std::vector<double>& freqVals,
    const std::vector<double>& ampVals, 
    const Detector& detector, 
    double Tyear) 
{
    if (freqVals.size() != ampVals.size()) 
    {
        throw std::invalid_argument("freqVals and ampVals must have the same size");
    }

    if (freqVals.empty()) 
    {
        throw std::invalid_argument("Input vectors cannot be empty");
    }

    size_t n = freqVals.size();
    std::vector<double> integrand(n);

    for (size_t i = 0; i < n; ++i) 
    {
        double f      = freqVals[i];
        double signal = ampVals[i];
        double noise  = detector.omega_noise(f);

        if (noise <= 0.0) 
        {
            throw std::runtime_error(detector.name() + ": noise PSD is non-positive at f = " + std::to_string(f));
        }

        integrand[i] = std::pow(signal / noise, 2);
    }

    double integral = simpson_integrate(freqVals, integrand);
    double Tsec = Tyear * 365.25 * 24.0 * 60.0 * 60.0;

    return std::sqrt(detector.ndet() * Tsec * integral);
}

std::vector<SNRResult> calculate_all_snrs
(
    const std::vector<double>& freqVals,
    const std::vector<double>& ampVals, 
    const std::vector<const Detector*>& detectors, 
    double Tyear) 
{

    std::vector<SNRResult> results;
    results.reserve(detectors.size());

    for (const Detector* det : detectors) 
    {
        results.push_back({ det->name(), calculate_snr(freqVals, ampVals, *det, Tyear) });
    }

    return results;
}


std::vector<SNRResult> get_SNR(const std::vector<double>& freqVals, const std::vector<double>& ampVals, double Tyear) 
{
    static const LISA lisa;
    static const DECIGO decigo;

    const std::vector<const Detector*> detectors = { &lisa, &decigo };

    return calculate_all_snrs(freqVals, ampVals, detectors, Tyear);
}

double LISA_snr(const std::vector<double>& freqVals, const std::vector<double>& ampVals, double Tyear) 
{
    static const LISA lisa;
    return calculate_snr(freqVals, ampVals, lisa, Tyear);
}

double DECIGO_snr(const std::vector<double>& freqVals, const std::vector<double>& ampVals, double Tyear) 
{
    static const DECIGO decigo;
    return calculate_snr(freqVals, ampVals, decigo, Tyear);
}

// double BBO_snr(const std::vector<double>& freqVals, const std::vector<double>& ampVals, double Tyear) 
// {
//     static const BBO bbo;
//     return calculate_snr(freqVals, ampVals, bbo, Tyear);
// }