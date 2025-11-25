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
#include <boost/math/quadrature/trapezoidal.hpp>

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

void create_fluid_integrand_splines(const FluidProfile& prof, alglib::spline1dinterpolant& f_sin_spline, alglib::spline1dinterpolant& f_cos_spline, alglib::spline1dinterpolant& l_sin_spline)
{
    const auto xi_vals = prof.xi_vals();
    const auto v_vals = prof.v_vals();
    const auto la_vals = prof.la_vals();

    std::vector<double> x_vals, f_sin_vals, f_cos_vals, l_sin_vals;
    for ( size_t i = 0; i < xi_vals.size(); i += 2) {
        x_vals.push_back(xi_vals[i]);
        f_sin_vals.push_back( - v_vals[i]);
        f_cos_vals.push_back(v_vals[i] * xi_vals[i]);
        l_sin_vals.push_back(la_vals[i] * xi_vals[i]);
    }

    alglib::real_1d_array x_arr, f_sin_arr, f_cos_arr, l_sin_arr;
    x_arr.setcontent(x_vals.size(), x_vals.data());
    f_sin_arr.setcontent(f_sin_vals.size(), f_sin_vals.data());
    f_cos_arr.setcontent(f_cos_vals.size(), f_cos_vals.data());
    l_sin_arr.setcontent(l_sin_vals.size(), l_sin_vals.data());

    try {
        alglib::spline1dbuildcubic(x_arr, f_sin_arr, f_sin_spline);
        alglib::spline1dbuildcubic(x_arr, f_cos_arr, f_cos_spline);
        alglib::spline1dbuildcubic(x_arr, l_sin_arr, l_sin_spline);
    } catch (const alglib::ap_error& e) {
        throw std::runtime_error(std::string("Error building spline in prof_ints_fl: ") + e.msg);
    } catch (...) {
        throw std::runtime_error("Unknown error building spline in prof_ints_fl");
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

    alglib::spline1dinterpolant f_sin_spline, f_cos_spline, l_sin_spline;
    create_fluid_integrand_splines(prof, f_sin_spline, f_cos_spline, l_sin_spline);

    LevinIntegrator levin(16);
    boost::math::quadrature::gauss<double, 15> integrator;

    #pragma omp parallel for
    for (size_t j = 0; j < M; ++j) {
        const double chi = chi_vals[j];
        const double inv_chi = 1.0 / chi;

        const double chi_threshold = 1e1;

        double f_sin_int = 0.0;
        double f_cos_int = 0.0;
        double l_sin_int = 0.0;

        const double xi_min = prof.xi_min();
        const double xi_max = prof.xi_max();
        const double vw = prof.params()->vw();
        const auto mode = prof.mode();

        if ( chi < chi_threshold ) {

            auto integrand_f_sin = [&](double xi) -> double {

                const auto sin_amp = alglib::spline1dcalc(f_sin_spline, xi);

                if(sin_amp == 0) {return 0.0;}

                const auto chi_xi = chi * xi;
                const auto sin_cx = std::sin(chi_xi);

                return sin_amp * sin_cx;
            };

            auto integrand_f_cos = [&](double xi) -> double {

                const auto cos_amp = alglib::spline1dcalc(f_cos_spline, xi);

                if(cos_amp == 0) {return 0.0;}

                const auto chi_xi = chi * xi;
                const auto cos_cx = std::cos(chi_xi);

                return cos_amp * cos_cx;
            };

            auto integrand_l_sin = [&](double xi) -> double {

                const auto sin_amp = alglib::spline1dcalc(l_sin_spline, xi);

                if(sin_amp == 0) {return 0.0;}

                double sin_cx = std::sin(chi * xi);
                return sin_amp * sin_cx;
            };

            if ( mode == 0 || mode == 2) { 
                // deflagration or detonation
                f_cos_int = integrator.integrate(integrand_f_cos, xi_min, xi_max);
                f_sin_int = integrator.integrate(integrand_f_sin, xi_min, xi_max);
                l_sin_int = integrator.integrate(integrand_l_sin, xi_min, xi_max);
            } else {
                // hybrid, split regions up
                f_cos_int += integrator.integrate(integrand_f_cos, xi_min, vw);
                f_sin_int += integrator.integrate(integrand_f_sin, xi_min, vw);
                l_sin_int += integrator.integrate(integrand_l_sin, xi_min, vw);

                f_sin_int += integrator.integrate(integrand_f_sin, vw, xi_max);
                f_cos_int += integrator.integrate(integrand_f_cos, vw, xi_max);
                l_sin_int += integrator.integrate(integrand_l_sin, vw, xi_max);
            }

        } else {
            auto f_sin_func = [&](double xi) -> double {
                return alglib::spline1dcalc(f_sin_spline, xi);
            };
            
            auto f_cos_func = [&](double xi) -> double {
                return alglib::spline1dcalc(f_cos_spline, xi);
            };
            
            auto l_sin_func = [&](double xi) -> double {
                return alglib::spline1dcalc(l_sin_spline, xi);
            };

            if ( mode == 0 || mode == 2) { 
                // deflagration or detonation
                f_sin_int = levin.integrate_sin(f_sin_func, chi, xi_min, xi_max);
                f_cos_int = levin.integrate_cos(f_cos_func, chi, xi_min, xi_max);
                l_sin_int = levin.integrate_sin(l_sin_func, chi, xi_min, xi_max);
            } else {
                // hybrid, split regions up
                f_sin_int += levin.integrate_sin(f_sin_func, chi, xi_min, vw);
                f_cos_int += levin.integrate_cos(f_cos_func, chi, xi_min, vw);
                l_sin_int += levin.integrate_sin(l_sin_func, chi, xi_min, vw);

                f_cos_int += levin.integrate_cos(f_cos_func, chi, vw, xi_max);
                f_sin_int += levin.integrate_sin(f_sin_func, chi, vw, xi_max);
                l_sin_int += levin.integrate_sin(l_sin_func, chi, vw, xi_max);
            }
        }

        double fd_j = fac * inv_chi * (f_cos_int + inv_chi * f_sin_int);
        double l_j = fac * inv_chi * l_sin_int;

        const double lambda_min = prof.la_vals().front();
        if(lambda_min != 0) {
            const double cos_term = - xi_min * std::cos(xi_min * chi);
            const double sin_term = std::sin(xi_min * chi)*inv_chi;
            const double analytic_l = cos_term + sin_term;
            l_j += fac * inv_chi * inv_chi * lambda_min * analytic_l;
        }

        fd[j] = fd_j;
        l[j] = l_j;
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