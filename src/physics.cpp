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

// Galactic binary confusion noise in strain
// double gb_S(double f) {
//     // Convert to millihertz
//     double f_mHz = f * 1e3;
    
//     double A = 8e-38;
//     double f_ref = 1000.0;
//     double a = 0.138;
//     double b = -0.221;
//     double c = 0.521;
//     double d = 1.680;
//     double fk = 1.13;
    
//     double exponent = -(f_mHz / f_ref) * a - b * f_mHz * std::sin(c * f_mHz);
//     double tanh_term = 1.0 + std::tanh(d * (fk - f_mHz));
    
//     return A * std::pow(f_mHz, 7.0/3.0) * std::exp(exponent) * tanh_term;
// }

// // Galactic binary confusion noise in Omega_GW h^2
// double gb_omegahsq(double f) {
//     double Sn = gb_S(f);
//     return 4.0 * M_PI * M_PI / (3.0 * H_0 * H_0) * std::pow(f, 3) * Sn * h * h;
// }

// /**
//  * Extragalactic binary confusion noise in Omega_GW h^2
//  * @param f Frequency in Hz
//  * @param flag "central", "upper", or "lower" for different amplitude estimates
//  */
// double eb_omegahsq(double f, const std::string& flag) {
//     double amp = 8.9e-10;
    
//     if (flag == "upper") {
//         amp += 12.6e-10;
//     } else if (flag == "lower") {
//         amp -= 5.6e-10;
//     }
//     // else flag == "central", use base amplitude
    
//     double f_ref = 25.0;
//     return amp * std::pow(f / f_ref, 2.0/3.0) * h * h;
// }