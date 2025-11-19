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

#include <boost/math/quadrature/gauss_kronrod.hpp>

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

std::pair<std::vector<double>, std::vector<double>> fluid_profile_integrals(const std::vector<double>& chi_vals, const FluidProfile& prof) {

    const auto xi_vals = prof.xi_vals();
    const auto v_vals = prof.v_vals();
    const auto la_vals = prof.la_vals();

    // the resolution is high enough that neighbouring points have the samee y_values which
    // affects spline fitting.
    std::vector<double> low_res_xi_vals, low_res_v_vals, low_res_la_vals;
    for ( size_t i = 0; i < xi_vals.size(); i += 2) {
        low_res_xi_vals.push_back(xi_vals[i]);
        low_res_v_vals.push_back(v_vals[i]);
        low_res_la_vals.push_back(la_vals[i]);
    }

    alglib::real_1d_array xi_arr, v_arr, la_arr;
    xi_arr.setcontent(low_res_xi_vals.size(), low_res_xi_vals.data());
    v_arr.setcontent(low_res_v_vals.size(), low_res_v_vals.data());
    la_arr.setcontent(low_res_la_vals.size(), low_res_la_vals.data());

    alglib::spline1dinterpolant v_spline, la_spline;
    try {
        alglib::spline1dbuildcubic(xi_arr, v_arr, v_spline);
        alglib::spline1dbuildcubic(xi_arr, la_arr, la_spline);
    } catch (const alglib::ap_error& e) {
        throw std::runtime_error(std::string("Error building spline in prof_ints_fl: ") + e.msg);
    } catch (...) {
        throw std::runtime_error("Unknown error building spline in prof_ints_fl");
    }

    std::ofstream ofs_debug("spline_debug.csv");
    ofs_debug << "xi,v,v_sp,la,la_sp\n";
    for (size_t i = 0; i < low_res_xi_vals.size(); ++i) {
        ofs_debug << low_res_xi_vals[i] << "," << low_res_v_vals[i] << "," << alglib::spline1dcalc(v_spline, low_res_xi_vals[i]) << "," << low_res_la_vals[i] << "," << alglib::spline1dcalc(la_spline, low_res_xi_vals[i]) << "\n";
    }
    ofs_debug.close();
    
    const auto n = xi_vals.size();
    const auto m = chi_vals.size();
    const auto fac = 4.0 * M_PI;

    std::vector<double> fd(m); 
    std::vector<double> l(m);

    #pragma omp parallel
    {
        std::vector<double> fd_integrand(n);
        std::vector<double> l_integrand(n);

        #pragma omp for
        for (size_t j = 0; j < m; j++) {

            const auto chi = chi_vals[j];
            const auto inv_chi = 1.0 / chi;

            auto integrand_f_dash = [&](double xi) -> double {

                const auto v_p = alglib::spline1dcalc(v_spline, xi);

                if(v_p == 0) { return 0.0;}

                const auto chi_xi = chi * xi;
                const auto sin_cx = std::sin(chi_xi);
                const auto cos_cx = std::cos(chi_xi);

                return v_p * (xi * cos_cx - sin_cx * inv_chi);
            };

            auto integrand_l = [&](double xi) -> double {

                const auto lambda_p = alglib::spline1dcalc(la_spline, xi);

                if(lambda_p == 0) { return 0.0;}

                const auto chi_xi = chi * xi;
                const auto sin_cx = std::sin(chi_xi);
                const auto cos_cx = std::cos(chi_xi);

                return  lambda_p * xi * sin_cx;
            };

            const double tol = 1e-8;

            if (chi < 1e-10) {
                const int max_depth = 5;
                auto integrand_fd_small_chi = [&](double xi) -> double {
                    const auto v_p = alglib::spline1dcalc(v_spline, xi);
                    return -0.5 * v_p * xi * xi * xi;
                };
                fd[j] = fac * chi * boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand_fd_small_chi, xi_vals.front(), xi_vals.back(), max_depth, tol);
                
                auto integrand_l_small_chi = [&](double xi) -> double {
                    const auto lambda_p = alglib::spline1dcalc(la_spline, xi);
                    return lambda_p * xi * xi;
                };
                l[j] = fac * chi * boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand_l_small_chi, xi_vals.front(), xi_vals.back(), max_depth, tol);
            } else {
                const int max_depth = 8;
                fd[j] = fac * inv_chi * boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand_f_dash, xi_vals.front(), xi_vals.back(), max_depth, tol);
                l[j] = fac * inv_chi * boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand_l, xi_vals.front(), xi_vals.back(), max_depth, tol);
            }
        }
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