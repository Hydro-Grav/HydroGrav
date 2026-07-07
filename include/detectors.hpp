#ifndef INCLUDE_DETECTORS_HPP_H
#define INCLUDE_DETECTORS_HPP_H

#include <cmath>
#include <vector>
#include <string>
#include "constants.hpp"

static constexpr double H_0 = 67.8 / 3.086e19; // Hubble constant in SI units
static constexpr double h   = 0.678; // Reduced Hubble constant

// Base detector class
class Detector {
public:
    virtual ~Detector() = default;

    virtual std::string name() const = 0;

    // Effective number of independent detector pairs used in the SNR integral.
    // Override in derived classes where ndet != 1 (e.g. DECIGO).
    virtual double ndet() const { return 1.0; }

    virtual double Pn(double f) const = 0;
    virtual double R(double f) const = 0;

    // Transfer frequency f_* = c / (2 pi L)
    double fs(double L) const {
        return c / (2.0 * M_PI * L);
    }

    double Sn(double f) const {
        return Pn(f) / R(f);
    }

    // Noise sensitivity curve in Omega_GW h^2 units.
    double omega_noise(double f) const {
        return (2.0 * M_PI * M_PI / (3.0 * H_0 * H_0)) * std::pow(f, 3) * Sn(f) * h * h;
    }    
};
    
// LISA
class LISA : public Detector {
public:
    std::string name() const override { return "LISA"; }

private:
    static constexpr double L  = 2.5e9;
    const double fs_LISA = fs(L);

    double P_oms(double f) const {
        return (1.5e-11) * (1.5e-11) * (1.0 + std::pow(2e-3 / f, 4));
    }

    double P_acc(double f) const {
        double term1 = (3e-15) * (3e-15);
        double term2 = 1.0 + std::pow(0.4e-3 / f, 2);
        double term3 = 1.0 + std::pow(f / 8e-3, 4);
        return term1 * term2 * term3;
    }

    double Pn(double f) const {
        return (P_oms(f) + 4.0 * P_acc(f) / std::pow(2.0 * M_PI * f, 4)) / (L * L);
    }

    double R(double f) const {
        return (3.0 / 10.0) / (1.0 + 0.54 * std::pow(f / fs_LISA, 2));
    }
};

// DECIGO
class DECIGO : public Detector {
public:
    std::string name() const override { return "DECIGO"; }
    double ndet() const override { return 2.0; }  // 2 correlated detector pairs

private:
    // Physical constants
    static constexpr double L       = 1.0e3;           // Arm length (m)
    static constexpr double la      = 532e-9;          // Laser wavelength (m)
    static constexpr double M_mir   = 100.0;           // Mirror mass (kg)
    static constexpr double P_las   = 10.0;            // Laser output power (W)
    static constexpr double P_eff   = 6.68;            // Effective laser power (W)
    static constexpr double F_cav   = 10.18;           // Cavity finesse

    const double fs_DECIGO = fs(L);

    // Shot noise PSD
    double P_shot(double f) const {
        return (hbar * c * M_PI * la / P_eff)
             * std::pow(1.0 / (4.0 * F_cav * L), 2)
             * (1.0 + std::pow(f / fs_DECIGO, 2));
    }

    // Radiation pressure noise PSD
    double P_rad(double f) const {
        return (hbar * P_las / (c * M_PI * la))
             * std::pow(16.0 * F_cav / (M_mir * L), 2)
             * std::pow(1.0 / (2.0 * M_PI * f), 4)
             / (1.0 + std::pow(f / fs_DECIGO, 2));
    }

    // Acceleration noise PSD
    double P_acc(double f) const {
        return (hbar * P_las / (c * M_PI * la))
             * std::pow(16.0 * F_cav / (3.0 * M_mir * L), 2)
             * std::pow(1.0 / (2.0 * M_PI * f), 4);
    }

    // Total detector noise PSD
    double Pn(double f) const {
        return P_shot(f) + P_rad(f) + P_acc(f);
    }

    // Signal response function
    double R(double f) const {
        return (3.0 / 10.0) / (1.0 + 0.54 * std::pow(f / fs_DECIGO, 2));
    }

    // Strain noise spectral density Sn = Pn / R
    double Sn(double f) const {
        return Pn(f) / R(f);
    }
};

// BBO
class BBO : public Detector {
public:
    std::string name() const override { return "BBO"; }
    double ndet() const override { return 2.0; }

private:
    static constexpr double L = 5e+7;
    static constexpr double P_oms = (1.4e-17) * (1.4e-17);
    static constexpr double P_acc = (3.0e-17) * (3.0e-17);
    
    double Pn(double f) const {
        return 4.0 * (P_oms + P_acc / std::pow(2.0 * M_PI * f, 4)) / (L * L);
    }

    // TO DO: add response function
    double R(double f) const {
        return 0.0 * f;
    }
};

#endif // INCLUDE_DETECTORS_HPP_H