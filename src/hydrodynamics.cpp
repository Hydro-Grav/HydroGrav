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
#include <fftw3.h>
#include <finufft.h>

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

static void build_nodes_and_samples(const FluidProfile &prof,
    std::vector<double>& x, std::vector<double>& f_sin, std::vector<double>& f_cos, std::vector<double>& l_sin)
{
    const auto &xi_vals = prof.xi_vals();
    const auto &v_vals  = prof.v_vals();
    const auto &la_vals = prof.la_vals();

    if (xi_vals.size() != v_vals.size() || xi_vals.size() != la_vals.size())
        throw std::runtime_error("profile vectors size mismatch");

    x.clear(); f_sin.clear(); f_cos.clear(); l_sin.clear();

    for (size_t i = 0; i + 1 < xi_vals.size(); i += 2) {
        x.push_back(xi_vals[i]);
        f_sin.push_back(-v_vals[i]);
        f_cos.push_back(v_vals[i] * xi_vals[i]);
        l_sin.push_back(la_vals[i] * xi_vals[i]);
    }
}

static std::vector<double> trapezoid_weights(const std::vector<double>& x) {
    const size_t N = x.size();
    if (N == 0) return {};
    std::vector<double> w(N);
    if (N == 1) { w[0] = 0.0; return w; }

    w[0] = 0.5 * (x[1] - x[0]);
    for (size_t i = 1; i + 1 < N; ++i) {
        w[i] = 0.5 * (x[i+1] - x[i-1]);
    }
    w[N-1] = 0.5 * (x[N-1] - x[N-2]);
    return w;
}

/*
    Below is a modification of this routine that utilises FFT for 
    small chi values, and resorts to a levine integration method for 
    large chi. 

    This code was produced using GPT 5.0 and Claude Sonnet 4.5
*/
std::pair<std::vector<double>, std::vector<double>>
fluid_profile_integrals(const std::vector<double>& chi_vals, const FluidProfile& prof)
{
    const double chi_threshold = 500.0;
    
    std::vector<double> x, f_sin, f_cos, l_sin;
    build_nodes_and_samples(prof, x, f_sin, f_cos, l_sin);

    const size_t N = x.size();
    if (N == 0) return {{},{}};

    const double fac = 4.0 * M_PI;
    const size_t M = chi_vals.size();
    
    std::vector<double> fd(M, 0.0), l(M, 0.0);

    std::vector<size_t> low_chi_indices, high_chi_indices;
    std::vector<double> low_chi_vals, high_chi_vals;
    
    for (size_t j = 0; j < M; ++j) {
        if (chi_vals[j] < chi_threshold) {
            low_chi_indices.push_back(j);
            low_chi_vals.push_back(chi_vals[j]);
        } else {
            high_chi_indices.push_back(j);
            high_chi_vals.push_back(chi_vals[j]);
        }
    }

    if (!low_chi_vals.empty()) {
        std::vector<double> w = trapezoid_weights(x);

        std::vector<std::complex<double>> c_fsin(N), c_fcos(N), c_lsin(N);
        for (size_t n = 0; n < N; ++n) {
            c_fsin[n] = std::complex<double>( f_sin[n] * w[n], 0.0 );
            c_fcos[n] = std::complex<double>( f_cos[n] * w[n], 0.0 );
            c_lsin[n] = std::complex<double>( l_sin[n] * w[n], 0.0 );
        }

        const size_t M_low = low_chi_vals.size();

        std::vector<double> xbuf = x;
        std::vector<double> kbuf = low_chi_vals;
        std::vector<std::complex<double>> out_fsin(M_low), out_fcos(M_low), out_lsin(M_low);

        const int iflag = +1;
        const double eps = 1e-12;
        int ier;

        ier = finufft1d3((int64_t)N, xbuf.data(),
                         reinterpret_cast<std::complex<double>*>(c_fsin.data()),
                         iflag, eps, (int64_t)M_low, kbuf.data(),
                         reinterpret_cast<std::complex<double>*>(out_fsin.data()), nullptr);
        if (ier != 0) throw std::runtime_error("FINUFFT error in f_sin transform");

        ier = finufft1d3((int64_t)N, xbuf.data(),
                         reinterpret_cast<std::complex<double>*>(c_fcos.data()),
                         iflag, eps, (int64_t)M_low, kbuf.data(),
                         reinterpret_cast<std::complex<double>*>(out_fcos.data()), nullptr);
        if (ier != 0) throw std::runtime_error("FINUFFT error in f_cos transform");

        ier = finufft1d3((int64_t)N, xbuf.data(),
                         reinterpret_cast<std::complex<double>*>(c_lsin.data()),
                         iflag, eps, (int64_t)M_low, kbuf.data(),
                         reinterpret_cast<std::complex<double>*>(out_lsin.data()), nullptr);
        if (ier != 0) throw std::runtime_error("FINUFFT error in l_sin transform");

        for (size_t i = 0; i < M_low; ++i) {
            const size_t j = low_chi_indices[i];
            const double chi = chi_vals[j];
            const double inv_chi = 1.0 / chi;

            l[j] = fac * inv_chi * std::imag(out_lsin[i]);

            double term_cos = std::real(out_fcos[i]);
            double term_sin = std::imag(out_fsin[i]);
            fd[j] = fac * inv_chi * (term_cos + inv_chi * term_sin);
        }
    }

    if (!high_chi_vals.empty()) {
        const auto &xi_vals = prof.xi_vals();
        const auto &v_vals  = prof.v_vals();
        const auto &la_vals = prof.la_vals();

        std::vector<double> x_for_spline, v_for_spline, v_xi_for_spline, la_xi_for_spline;
        for (size_t i = 0; i + 1 < xi_vals.size(); i += 2) {
            x_for_spline.push_back(xi_vals[i]);
            v_for_spline.push_back(v_vals[i]);
            v_xi_for_spline.push_back(v_vals[i] * xi_vals[i]);
            la_xi_for_spline.push_back(la_vals[i] * xi_vals[i]);
        }

        alglib::real_1d_array x_arr, v_arr, v_xi_arr, la_xi_arr;
        x_arr.setcontent(x_for_spline.size(), x_for_spline.data());
        v_arr.setcontent(v_for_spline.size(), v_for_spline.data());
        v_xi_arr.setcontent(v_xi_for_spline.size(), v_xi_for_spline.data());
        la_xi_arr.setcontent(la_xi_for_spline.size(), la_xi_for_spline.data());

        alglib::spline1dinterpolant v_spline, v_xi_spline, la_xi_spline;
        alglib::spline1dbuildcubic(x_arr, v_arr, v_spline);
        alglib::spline1dbuildcubic(x_arr, v_xi_arr, v_xi_spline);
        alglib::spline1dbuildcubic(x_arr, la_xi_arr, la_xi_spline);

        const double xi_min = x_for_spline.front();
        const double xi_max = x_for_spline.back();
        const int N_colloc = 32;

        for (size_t i = 0; i < high_chi_vals.size(); ++i) {
            const size_t j = high_chi_indices[i];
            const double chi = chi_vals[j];
            const double inv_chi = 1.0 / chi;

            auto [l_sin_int, l_cos_int] = Levin::levin_integrate(la_xi_spline, chi, xi_min, xi_max, N_colloc, false);
            l[j] = fac * inv_chi * l_sin_int;

            auto [v_sin_int, v_cos_int] = Levin::levin_integrate(v_spline, chi, xi_min, xi_max, N_colloc, false);
            auto [vxi_sin_int, vxi_cos_int] = Levin::levin_integrate(v_xi_spline, chi, xi_min, xi_max, N_colloc, false);

            fd[j] = fac * inv_chi * (vxi_cos_int - inv_chi * v_sin_int);
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