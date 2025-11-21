// hydrodynamics.cpp

/*
TO DO:
- change short functions to inline functions and move to header
- check for domain errors when calling cublic spline in Ap_sq()
- could change f and l integrations to use profile.v_vals() rather than calling profile.v_prof(xi) every time (might be faster?)
    - i'm not sure if i would even need an interpolator if this is the case
- NOTE: on profile integrals f and l, xi is bound by 0<xi<1 so integration to infinity not strictly necessary (int from 0->1 and 0->inf agree to numerical precision)
    - might be safer to keep int to inf, but a bit quicker to just integrate to 1
- construct spline interpolation for f' and l integrands and use boost.math or integrate over a vector of integrand values using custom simpson/trapezoidal integrator?
    - for f' and l, they seem to be quite well behaved so both would be okay, but Ap_sq is highly oscillatory so i could reuse a custom integrator here (fitting spline to highly oscillatory function not good)
    - if adaptive integration (changing step-size for different parts of the func) is needed, calc vector of integrand values and construct a lambda from this somehow to use boost.math (avoids spline)
*/

#include <cmath>
#include <complex>
#include <functional>
#include <chrono>
#include <omp.h>
#include <fstream>

#include "maths_ops.hpp"
#include "profile.hpp"

namespace Hydrodynamics {

std::function<double(double)> lifetime_distribution_function(const std::string& nuc_type) {
    if (nuc_type == "exp") {
        return [](double Ttilde) -> double {
            return std::exp(-Ttilde);
        };
    } else if (nuc_type == "sim") {
        return [](double Ttilde) -> double {
            const auto exp_fac = - Ttilde * Ttilde * Ttilde / 6.0;
            return 0.5 * Ttilde * Ttilde * std::exp(exp_fac);
        };
    } else {
        throw std::invalid_argument("Invalid nucleation type: " + nuc_type);
    }
}


std::pair<std::vector<double>, std::vector<double>>
fluid_profile_integrals(const std::vector<double>& chi_vals, const FluidProfile& prof)
{
    const auto &xi_vals = prof.xi_vals();
    const auto &v_vals  = prof.v_vals();
    const auto &la_vals = prof.la_vals();

    const size_t N = xi_vals.size();
    if (N == 0) return {{},{}};

    const double fac = 4.0 * M_PI;
    const size_t M = chi_vals.size();
    
    std::vector<double> fd(M, 0.0), l(M, 0.0);

    #pragma omp parallel for
    for (size_t j = 0; j < M; ++j) {
        const double chi = chi_vals[j];
        const double inv_chi = 1.0 / chi;

        double f_sin_int = 0.0;
        double f_cos_int = 0.0;
        double l_sin_int = 0.0;

        for (size_t i = 0; i + 1 < N; ++i) {
            const double xi_i = xi_vals[i];
            const double xi_ip1 = xi_vals[i + 1];
            const double dx = xi_ip1 - xi_i;

            const double sin_chi_xi_i = std::sin(chi * xi_i);
            const double sin_chi_xi_ip1 = std::sin(chi * xi_ip1);

            const double y_fsin_i = -v_vals[i] * sin_chi_xi_i;
            const double y_fsin_ip1 = -v_vals[i + 1] * sin_chi_xi_ip1;
            f_sin_int += 0.5 * (y_fsin_i + y_fsin_ip1) * dx;

            const double y_fcos_i = v_vals[i] * xi_i * std::cos(chi * xi_i);
            const double y_fcos_ip1 = v_vals[i + 1] * xi_ip1 * std::cos(chi * xi_ip1);
            f_cos_int += 0.5 * (y_fcos_i + y_fcos_ip1) * dx;

            const double y_lsin_i = la_vals[i] * xi_i * sin_chi_xi_i;
            const double y_lsin_ip1 = la_vals[i + 1] * xi_ip1 * sin_chi_xi_ip1;
            l_sin_int += 0.5 * (y_lsin_i + y_lsin_ip1) * dx;
        }

        fd[j] = fac * inv_chi * (f_cos_int + inv_chi * f_sin_int);
        l[j] = fac * inv_chi * l_sin_int;
    }

    return {fd, l};
}

// |A_+|^2
std::vector<double> Ap_sq(const std::vector<double>& chi_vals, const FluidProfile& prof) {
    const auto csq = prof.params()->cpsq();
    const auto [fd_int, l_int] = fluid_profile_integrals(chi_vals, prof);
    const auto m = chi_vals.size();

    std::vector<double> Apsq(m);
    for (size_t j = 0; j < m; j++) {
        const auto f = fd_int[j];
        const auto l = l_int[j];

        Apsq[j] = 0.25 * (f*f + csq * l*l);
    }

    std::ofstream ofs("Ap_sq_debug.csv");
    ofs << "chi,Apsq,fd,l\n";
    for (size_t j = 0; j < m; j++) {
        ofs << chi_vals[j] << "," << Apsq[j] << "," << fd_int[j] << "," << l_int[j] << "\n";
    }
    ofs.close();

    return Apsq;
}

} // namespace Hydrodynamics