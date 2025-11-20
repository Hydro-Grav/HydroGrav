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
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

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

// struct IntegrandParams {
//     const alglib::spline1dinterpolant* spline;
// };

// double integrand_qawo(double xi, void *params)
// {
//     auto *P = static_cast<IntegrandParams*>(params);
//     return alglib::spline1dcalc(*(P->spline), xi) * xi;
// }

// const double compute_gsl_QAWO(const double& chi, const alglib::spline1dinterpolant& v_spline, const gsl_integration_qawo_enum type)
// {
//     const double xi_min = 0.01;
//     const double xi_max = 0.99;
//     const double L = xi_max - xi_min;

//     gsl_error_handler_t* old_handler = gsl_set_error_handler_off();

//     size_t workspace_size = 10000;
//     size_t n_levels = 150;

//     gsl_integration_workspace *w = gsl_integration_workspace_alloc(workspace_size);
//     gsl_integration_qawo_table *table = gsl_integration_qawo_table_alloc(chi, L, type, n_levels);

//     IntegrandParams params{ &v_spline };

//     gsl_function F;
//     F.function = &integrand_qawo;
//     F.params   = &params;

//     double result, error;

//     double abs_tol = 1e-4;
//     double rel_tol = 1e-4;

//     int status = gsl_integration_qawo(&F, xi_min, abs_tol, rel_tol, workspace_size, w, table, &result, &error);

//     if (status != GSL_SUCCESS) {
//         std::cerr << "GSL QAWO error: " << gsl_strerror(status) << " for chi = " << chi << ", error = " << error << "\n";
//         gsl_integration_qawo_table_free(table);
//         gsl_integration_workspace_free(w);
//         gsl_set_error_handler(old_handler);
//         return -1;
//     }

//     gsl_integration_qawo_table_free(table);
//     gsl_integration_workspace_free(w);

//     gsl_set_error_handler(old_handler);

//     return result;
// }

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

std::pair<std::vector<double>, std::vector<double>> fluid_profile_integrals(const std::vector<double>& chi_vals, const FluidProfile& prof) {

    alglib::spline1dinterpolant f_sin_spline, f_cos_spline, l_sin_spline;
    create_fluid_integrand_splines(prof, f_sin_spline, f_cos_spline, l_sin_spline);

    const double xi_min = prof.xi_vals().front();
    const double xi_max = prof.xi_vals().back();
    
    const auto m = chi_vals.size();
    const auto fac = 4.0 * M_PI;

    std::vector<double> fd(m); 
    std::vector<double> l(m);

    #pragma omp parallel
    {
        #pragma omp for
        for (size_t j = 0; j < m; j++) {

            const auto chi = chi_vals[j];
            const auto inv_chi = 1.0 / chi;

            auto integrand_f_dash = [&](double xi) -> double {

                const auto sin_term = alglib::spline1dcalc(f_sin_spline, xi);
                const auto cos_term = alglib::spline1dcalc(f_cos_spline, xi);

                const auto chi_xi = chi * xi;
                const auto sin_cx = std::sin(chi_xi);
                const auto cos_cx = std::cos(chi_xi);

                return sin_term * inv_chi * sin_cx + cos_term * cos_cx;
            };

            auto integrand_l = [&](double xi) -> double {

                const auto sin_term = alglib::spline1dcalc(l_sin_spline, xi);

                if(sin_term == 0) {return 0.0;}

                double sin_cx = std::sin(chi * xi);

                return sin_term * sin_cx;
            };

            double abs_error = 1e-8 / (1.0 + chi);
            double rel_error = 1e-8;
            if(chi < 5e2) {
                boost::math::quadrature::gauss_kronrod<double, 1000> integrator;
                l[j] = fac * inv_chi * integrator.integrate(integrand_l, xi_min, xi_max, abs_error, rel_error);
                fd[j] = fac * inv_chi * integrator.integrate(integrand_f_dash, xi_min, xi_max, abs_error, rel_error);
            } else {
                l[j] = fac * inv_chi * Levin::levin_integrate(l_sin_spline, chi, xi_min, xi_max, 64, false).first;
                fd[j] = fac * inv_chi * (inv_chi * Levin::levin_integrate(f_sin_spline, chi, xi_min, xi_max, 64, false).first + Levin::levin_integrate(f_cos_spline, chi, xi_min, xi_max, 64, false).second);
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