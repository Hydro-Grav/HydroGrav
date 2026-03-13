/**
 * @file ssm.cpp
 * @brief Gravitational wave calculator using the sound shell model (see arXiv:1608.04735)
 */

#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <filesystem>

#include "ap.h"
#include "interpolation.h"

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "config.hpp"
#include "maths_ops.hpp"
#include "phasetransition.hpp"
#include "profile.hpp"
#include "ssm.hpp"
#include "sici.hpp"

#ifdef ENABLE_MATPLOTLIB
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

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

    const double xi_min = prof.xi_min();
    const double xi_max = prof.xi_max();
    const double vw = prof.params()->vw();
    const auto mode = prof.mode();

    const double lambda_min = prof.la_vals().front();

    const size_t N = xi_vals.size();
    if (N == 0) return {{},{}};

    const double fac = 4.0 * M_PI;
    const size_t M = chi_vals.size();
    
    std::vector<double> fd(M, 0.0), l(M, 0.0);

    alglib::spline1dinterpolant f_sin_spline, f_cos_spline, l_sin_spline;
    create_fluid_integrand_splines(prof, f_sin_spline, f_cos_spline, l_sin_spline);

    FilonQuadrature levin(config::filon_polynomial_order);
    boost::math::quadrature::gauss<double, config::fd_l_gauss_legendre_samples> integrator;

    #pragma omp parallel for
    for (size_t j = 0; j < M; ++j) {
        const double chi = chi_vals[j];
        const double inv_chi = 1.0 / chi;

        double f_sin_int = 0.0;
        double f_cos_int = 0.0;
        double l_sin_int = 0.0;

        // evaluate fluid integrals over xi_min < xi < xi_max for given chi
        if ( chi < config::chi_threshold ) {

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
        double l_j;

        // analytic evaluation of l(chi) for 0 < xi < xi_min generates a lot of noise
        // at small chi since the analytic result is proportional to 1/chi^3, so we do
        // this part of the integration numerically
        if (chi < 1.0) { // numeric evaluation for small chi
            double l_spline_part = l_sin_int;
            if (lambda_min != 0.0) {
                auto const_integrand = [&](double xi) -> double {
                    return lambda_min * xi * std::sin(chi * xi);
                };

                double l_extra = levin.integrate_sin(const_integrand, chi, 0.0, xi_min);
                l_spline_part += l_extra;
            }
            l_j = fac * inv_chi * l_spline_part;
        } else { // analytic evaluation for large chi
            l_j = fac * inv_chi * l_sin_int;
            if(lambda_min != 0) {
                const double cos_term = - xi_min * std::cos(xi_min * chi);
                const double sin_term = std::sin(xi_min * chi)*inv_chi;
                const double analytic_l = cos_term + sin_term;
                l_j += fac * inv_chi * inv_chi * lambda_min * analytic_l;
            }
        }

        fd[j] = fd_j;
        l[j] = l_j;
    }

    return {fd, l};
}

std::vector<double> Ap_sq(const std::vector<double>& chi_vals, const FluidProfile& prof) {
    // NOTE: csq here defined in terms of background fluid (i.e. csq in symmetric phase at T=TN)
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

namespace Spectrum {

/***** PowerSpec class *****/

// Define ctors
// PowerSpec::PowerSpec(const std::vector<double>& K_vals, std::vector<double>& P_vals, const PhaseTransition::PTParams& params)
PowerSpec::PowerSpec(const std::vector<double>& K_vals, std::vector<double>& P_vals, const Hydrodynamics::FluidProfile& profile, double dtau)
    : freq_vals_(), K_vals_(K_vals), P_vals_(P_vals),
      dtau_(dtau),
      profile_(profile), params_(profile.params()) {
        if (K_vals.size() != P_vals.size()) {
            throw std::invalid_argument("PowerSpec: k and P vectors must be the same size!");
        }

        // store frequency values
        for (const auto& K : K_vals_) {
            const auto HsRs = params_->un().Hs() * params_->Rs();
            const auto freq_val = (2.6e-5) * (K / 10.0) * (1.0 / HsRs) * (params_->un().Ts() / 100.0) * std::pow(params_->un().gs() / 100.0, 1.0 / 6.0);
            freq_vals_.push_back(freq_val);
        }
    }

std::pair<double, double> PowerSpec::peak_vals() const {
    const auto &Pv = P();
    const auto it = std::max_element(Pv.begin(), Pv.end());
    const auto idx = std::distance(Pv.begin(), it);
    return {freq_vals_[idx], Pv[idx]};
}

void PowerSpec::write(const std::string& filename) const {
    namespace fs = std::filesystem;

    std::cout << "Writing power spectrum to disk... ";

    // check directory exists
    fs::path filepath(filename);
    fs::path dir = filepath.parent_path();
    if (!dir.empty() && !fs::exists(dir)) {
        if (!fs::create_directories(dir)) {
            std::cerr << "Failed to create directory: " << dir << "\n";
            return;
        }
    }

    std::ofstream file(filename);

    // header
    file << "# Power spectrum generated by DeepPhase\n";
    file << "# Universe Parameters:\n";
    file << "# T0 = " << params_->un().T0() << "\n";
    file << "# Ts = " << params_->un().Ts() << "\n";
    file << "# H0 = " << params_->un().H0() << "\n";
    file << "# Hs = " << params_->un().Hs() << "\n";
    file << "# g0 = " << params_->un().g0() << "\n";
    file << "# gs = " << params_->un().gs() << "\n";
    file << "#\n";

    file << "# Phase Transition Parameters:\n";
    file << "# vw = " << params_->vw() << "\n";
    file << "# alpha_N = " << params_->alN() << "\n";
    file << "# TN = " << params_->TN() << "\n";
    file << "# beta = " << params_->beta() << "\n";
    file << "# betaHs = " << params_->betaHs() << "\n";
    file << "# Rs = " << params_->Rs() << "\n";
    file << "# nuc_type = " << params_->nuc_type() << "\n";
    file << "# dtau/Rs (sound wave duration) = " << dtau() / params_->Rs() << "\n";
    file << "# f_peak = " << peak_vals().first << "\n";
    file << "#\n";

    file << "# Fluid Parameters:\n";
    file << "# eos = " << params_->eos_to_string() << "\n";
    file << "# cp = " << std::sqrt(params_->cpsq()) << "\n";
    file << "# cm = " << std::sqrt(params_->cmsq()) << "\n";
    file << "# mode = " << profile_.mode_str() << "\n";
    file << "# -----------------------------------------\n";
    file << "#\n";

    // data
    file << "freq,K,P\n";

    for (size_t i = 0; i < K_vals_.size(); ++i) {
        file << freq_vals_[i] << "," << K_vals_[i] << "," << P_vals_[i] << "\n";
    }
    file.close();
    std::cout << "Saved to " << filename << "!\n";

    return;
}

#ifdef ENABLE_MATPLOTLIB
void PowerSpec::plot(const std::string& filename) const {
    plt::figure_size(800, 600);
    plt::loglog(K(), P(), "k-");
    plt::suptitle("vw = " + to_string_with_precision(params_->vw()) + ", alN = " + to_string_with_precision(params_->alN()));
    plt::xlabel("K=kRs");
    plt::ylabel("Omega_GW(K)");
    plt::xlim(K().front(), K().back());
    plt::grid(true);
    plt::save(filename);

    std::cout << "Saved to '" << filename << "'" << std::endl;

    return;
}
#endif

CubicSpline<double> PowerSpec::interpolate() const {
    return CubicSpline(K(), P());
}

// PowerSpec [op] Scalar arithmetic (can't use copy/move assignments if passing in PTParams to PowerSpec)
PowerSpec operator*(const PowerSpec& spec, double scalar) {
    std::vector<double> scaled_P;
    scaled_P.reserve(spec.P().size());

    for (double p : spec.P()) {
        scaled_P.push_back(p * scalar);
    }

    return PowerSpec(spec.K(), scaled_P, spec.profile());
}

PowerSpec operator*(double scalar, const PowerSpec& spec) {
    return spec * scalar;
}

PowerSpec& PowerSpec::operator*=(double scalar) {
    // for (auto& p : data_.second) {
    for (auto& p : P_vals_ ) {
        p *= scalar;
    }
    return *this;
}

PowerSpec operator/(const PowerSpec& spec, double scalar) {
    if (scalar == 0)
        throw std::invalid_argument("PowerSpec: Division by zero!");
    return spec * (1.0 / scalar);
}

PowerSpec& PowerSpec::operator/=(double scalar) {
    if (scalar == 0.0)
        throw std::invalid_argument("PowerSpec: Division by zero!");

    // for (auto& p : data_.second) {
    for (auto& p : P_vals_) {
        p /= scalar;
    }
    return *this;
}
/***************************/

// minimises pt(k,p,z=1)
double find_min_pt(const std::vector<double>& k_vals, const std::vector<double>& p_vals) {
    size_t i = 0, j = 0;
    double min_pt = std::numeric_limits<double>::infinity();

    while (i < k_vals.size() && j < p_vals.size()) {
        double k = k_vals[i];
        double p = p_vals[j];
        double diff = k - p;
        double pt = std::abs(diff);

        if (pt > 1e-10) {  // skip exact matches
            if (pt < min_pt) {
                min_pt = pt;
            }
        }

        // Move the pointer for the smaller value
        if (k < p)
            ++i;
        else
            ++j;
    }

    return min_pt;
}

void build_kinetic_spectrum_spline(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& profile, alglib::spline1dinterpolant& log_zk_spline) 
{
    const auto kinetic_spectrum = zetaKin(kRs_vals, profile);

    // kinetic_spectrum.write("zetakin_debug.csv");

    alglib::real_1d_array x_arr, y_arr;
    x_arr.setlength(kinetic_spectrum.K().size());
    y_arr.setlength(kinetic_spectrum.K().size());
    for (size_t i = 0; i < kinetic_spectrum.K().size(); i++) 
    {
        x_arr[i] = kinetic_spectrum.K()[i];

        const auto P_val = kinetic_spectrum.P()[i];
        if (P_val <= 0.0 || std::isnan(P_val) || std::isinf(P_val)) 
        {
            y_arr[i] = -700;
        } else 
        {
            y_arr[i] = std::log(P_val);
        }
    }

    try 
    {
        alglib::spline1dbuildcubic(x_arr, y_arr, log_zk_spline);
    } catch (const alglib::ap_error& e) 
    {
        std::cerr << "ALGLIB error building spline for zetaKin(ptRs): " << e.msg << std::endl;
        throw;
    } catch (const std::exception& e) 
    {
        std::cerr << "Error building spline for zetaKin(ptRs): " << e.what() << std::endl;
        throw;
    } catch (...) 
    {
        std::cerr << "Unknown error building spline for zetaKin(ptRs)" << std::endl;
        throw;
    }
}

/*** GW power spectrum ***/
PowerSpec GWSpec(const std::vector<double>& kRs_vals, const PhaseTransition::PTParams& params, const bool calc_dtau) {

    const auto ti = std::chrono::high_resolution_clock::now();

    // std::cout << "[DEBUG] Config values: "
    //           << "pRs_samples=" << config::pRs_samples 
    //           << ", z_samples=" << config::z_samples 
    //           << ", pRs_tolerance=" << config::pRs_tolerance 
    //           << ", z_tolerance=" << config::z_tolerance << "\n";

    const Hydrodynamics::FluidProfile profile(params);

    const auto cs = std::sqrt(params.cpsq());
    const auto Rs_inv = 1.0 / params.Rs();
    
    const auto tau_s = params.tau_s();
    const auto dtau = (calc_dtau) ? get_nl_timescale(profile) : dtau_approx(params);
    const auto tau_fin = tau_s + dtau;

    const auto nk = kRs_vals.size();

    const auto pRs_vals = logspace(log10(config::pRs_minimum), log10(config::pRs_maximum), config::n_pRs); // P = p*Rs

    const auto kinetic_spectrum_spline_lower_bound = (1.0 - config::kinetic_spectrum_spline_factor) * find_min_pt(kRs_vals, pRs_vals);
    const auto kinetic_spectrum_spline_upper_bound = (1.0 + config::kinetic_spectrum_spline_factor) * ptilde(kRs_vals.back(), pRs_vals.back(), -1.0);

    const auto kinetic_spectrum_K_values = logspace(log10(kinetic_spectrum_spline_lower_bound), log10(kinetic_spectrum_spline_upper_bound), config::kinetic_spectrum_spline_points);

    alglib::spline1dinterpolant log_zk_spline;
    build_kinetic_spectrum_spline(kinetic_spectrum_K_values, profile, log_zk_spline);

    std::cout << "Calculating gravitational wave power spectrum...\n";

    const auto prefac = gw_prefac(kRs_vals, profile);

    std::vector<double> GW_P_vals(nk);
    
    #pragma omp parallel 
    {
        #pragma omp for schedule(static)
        for (size_t kk = 0; kk < nk; kk++ ) 
        {
            const auto kRs = kRs_vals[kk];
            const auto k = kRs * Rs_inv;
            const auto kRs3 = kRs * kRs * kRs;

            auto pRs_integrand = [&](double log_pRs) -> double
            {
                const auto pRs = exp(log_pRs);
                const auto p = pRs * Rs_inv;
                const auto pRs_sq = pRs*pRs;
                const auto zk_pRs_val = std::exp(alglib::spline1dcalc(log_zk_spline, pRs));
                const auto zk_pRs_fac = zk_pRs_val * pRs_sq;

                auto z_integrand = [&](double z) -> double 
                {
                    const auto ptRs = ptilde(kRs, pRs, z);

                    if (ptRs == 0.0 
                        || ptRs < kinetic_spectrum_spline_lower_bound 
                        || ptRs > kinetic_spectrum_spline_upper_bound) 
                    {
                        return 0.0;
                    }
                
                    const auto dlt = dlt_SSM(k, p, ptRs * Rs_inv, cs, tau_s, tau_fin);      
                    double zk_ptRs_val = std::exp(alglib::spline1dcalc(log_zk_spline, ptRs));

                    const auto z_fac = 1.0 - z*z;
                    const auto z_fac2 = z_fac * z_fac;
                    const auto ptRs4_inv = 1.0 / (ptRs * ptRs * ptRs * ptRs);

                    return z_fac2 * ptRs4_inv * zk_ptRs_val * dlt;
                };

                const double z_result = boost::math::quadrature::gauss_kronrod<double, config::z_samples>::integrate(z_integrand, -1.0, 1.0, config::z_max_refinements, config::z_tolerance);

                return pRs * zk_pRs_fac * z_result;
            };

            double pRs_result = boost::math::quadrature::gauss_kronrod<double, config::pRs_samples>::integrate(pRs_integrand, log(config::pRs_minimum), log(config::pRs_maximum), config::pRs_max_refinements, config::pRs_tolerance);
            GW_P_vals[kk] = prefac * kRs3 * pRs_result;
        }
    }

    std::cout << "Gravitational power spectrum constructed!\n";

    /***************************** CLOCK ******************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "GWSpec Timer: " << duration.count() << " s" << std::endl;
    /******************************************************************/

    return PowerSpec(kRs_vals, GW_P_vals, profile, dtau);
}
/***************************/

double dlt_SSM(double k, double p, double pt, const double cs, const double tau_s, const double tau_fin) {
    const auto ptcs = pt * cs;
    const auto pcs = p * cs;

    auto dlt = 0.0;
    for (int m = -1; m < 2; m+=2) { // m = {-1,1}
        const auto pmn_1 = pcs + m * ptcs;
        for (int n = -1; n < 2; n+=2) { // n = {-1,1}
            const auto pmn = pmn_1 + n * k; // pmn = (p + m*pt)*cs + n*k

            const auto x1 = pmn * tau_fin;
            const auto x2 = pmn * tau_s;

            double Si_x1, Ci_x1, Si_x2, Ci_x2;

            alglib::sinecosineintegrals(x1, Si_x1, Ci_x1);
            alglib::sinecosineintegrals(x2, Si_x2, Ci_x2);

            const auto dSi = Si_x1 - Si_x2;
            const auto dCi = Ci_x1 - Ci_x2;

            dlt += 0.25 * (dCi * dCi + dSi * dSi);
        }
    }
    return dlt;
}

/*** dlt spectrum ***/
// inline this later?
double ff(double tau_m, double kcs) {
    // kcs = k*cs -> ff called this way to make dlt faster
    return std::cos(kcs * tau_m); // for SSM -> NEED TO UPDATE THIS
}

/*** Kinetic spectrum ***/
// avoids duplicating fluid profile in GWSpec
PowerSpec Ekin(const std::vector<double>& kRs_vals, const PhaseTransition::PTParams& params) {
    return Ekin(kRs_vals, Hydrodynamics::FluidProfile(params));
}

PowerSpec Ekin(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& prof) {
    const auto beta = prof.params()->beta();
    const auto Rs = prof.params()->Rs();
    const auto nuc_type = prof.params()->nuc_type();

    auto lt_dist = Hydrodynamics::lifetime_distribution_function(nuc_type);

    const auto nk = kRs_vals.size();
    std::vector<double> P_vals(nk);

    const auto chi_vals = logspace(log10(config::chi_min), log10(config::chi_max), config::chi_points);
    const auto n = chi_vals.size();

    const auto Apsq = Hydrodynamics::Ap_sq(chi_vals, prof);

    const auto fac1 = beta * Rs * Rs / (2.0 * M_PI * M_PI);

    alglib::real_1d_array chi_arr, Apsq_arr;
    chi_arr.setcontent(n, chi_vals.data());
    Apsq_arr.setcontent(n, Apsq.data());
    alglib::spline1dinterpolant Apsq_spline;
    alglib::spline1dbuildcubic(chi_arr, Apsq_arr, Apsq_spline);

    #pragma omp parallel for
    for (size_t kk = 0; kk < nk; kk++) {
        const auto kRs = kRs_vals[kk];
        const auto kRs_inv = 1.0 / kRs;

        const auto fac2 = fac1 * power(kRs_inv, 5);
        const auto fac3 = beta * Rs * kRs_inv;

        auto integrand = [&](double log_chi) -> double 
        {
            const double chi = std::exp(log_chi);
            const double Apsq_val = alglib::spline1dcalc(Apsq_spline, chi);
            const double T_tilde = fac3*chi;
            return lt_dist(T_tilde) * power(chi, 7) * Apsq_val;
        };
        
        const double log_chi_min = std::log(config::chi_min);
        const double log_chi_max = std::log(config::chi_max);

        boost::math::quadrature::gauss_kronrod<double, config::Ekin_samples> integrator;
        const auto y = integrator.integrate(integrand, log_chi_min, log_chi_max, config::Ekin_max_refinements, config::Ekin_tolerance);

        P_vals[kk] = fac2 * y;
    }

    return PowerSpec(kRs_vals, P_vals, prof);
}

PowerSpec norm_spec(const PowerSpec& spec) {
    const auto spec_max = spec.peak_vals().second;
    if (spec_max == 0.0) {
        throw std::runtime_error("Division by zero in norm_spec from spec.peak_vals().second = 0");
    } else if (isnan(spec_max)) {
        throw std::runtime_error("In norm_spec: spec.peak_vals().second = nan");
    }

    const auto zk = spec / spec_max;
    if (abs(zk.peak_vals().second - 1.0) > 1e-15) {
        throw std::runtime_error("In norm_spec: Power spectrum failed normalisation test");
    }

    return zk;
}

PowerSpec zetaKin(const std::vector<double>& kRs_vals, const PhaseTransition::PTParams& params) {
    const auto Ek = Ekin(kRs_vals, params);
    return norm_spec(Ek);
}

PowerSpec zetaKin(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& prof) {
    const auto Ek = Ekin(kRs_vals, prof);
    return norm_spec(Ek);
}
/***************************/

// computes timescale for non-linearities to appear in fluid (used as sound wave duration)
double get_nl_timescale(const Hydrodynamics::FluidProfile& prof) {
    // tau_nl ~ Rs / sqrt(Omega_K)
    // Omega_K = avg kinetic energy of sound waves
    const auto Rs = prof.params()->Rs();
    const auto kRs_vals = logspace(log10(config::kRs_minimum), log10(config::kRs_maximum), config::n_kRs);

    const auto Ek = Ekin(kRs_vals, prof);
    const auto Ek_int = simpson_integrate(Ek.K(), Ek.P());    

    return std::sqrt(Rs * Rs * Rs / Ek_int);
}

// approximation used for dtau in arXiv:2308.12943
double dtau_approx(const PhaseTransition::PTParams& params) {
    return 10.0 * params.Rs();
}

double gw_prefac(double Ekin_max, double Rs, double wNeN_rat, double T0, double Ts, double H0, double Hs, double g0, double gs) {
    // Transfer function (redshift of spectrum - eq 13 arXiv:2308.12943)
    const auto g0gs_rat = g0 / gs;
    const auto TH_rat = (T0 * T0 / H0) / (Ts * Ts / Hs);
    const auto TGW = std::pow(g0gs_rat, 4./3.) * TH_rat * TH_rat;

    // Normalised kinetic energy density OmegaK / KK (eq 42 arXiv:2308.12943)
    // OmegaK = total kinetic energy density, KK = critical energy density
    const auto OmegaK_KK = Ekin_max / Rs;

    // std::cout << "Ekin_max=" << Ekin_max << ", wNeN_rat=" << wNeN_rat << "\n";
    // std::cout << "prefac=" << 3.0 * wNeN_rat * wNeN_rat * TGW * OmegaK_KK * OmegaK_KK << "\n";
    
    return 3.0 * wNeN_rat * wNeN_rat * TGW * OmegaK_KK * OmegaK_KK;
}

double gw_prefac(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& profile) {
    const auto params = profile.params();
    const auto un = params->un();

    const auto Ek = Ekin(kRs_vals, profile);

    return gw_prefac(Ek.peak_vals().second, params->Rs(), params->wNeN_rat(), un.T0(), un.Ts(), un.H0(), un.Hs(), un.g0(), un.gs());
}

#ifdef ENABLE_MATPLOTLIB
void plot_spectra(const PowerSpec& gw_spec_bag, const PowerSpec& gw_spec_munu, const PowerSpec& gw_spec_veff, const std::string& filename, const double f_min, const double f_max) {
    const auto freq_bag = gw_spec_bag.freq();
    const auto P_bag = gw_spec_bag.P();

    const auto freq_munu = gw_spec_munu.freq();
    const auto P_munu = gw_spec_munu.P();

    const auto freq_veff = gw_spec_veff.freq();
    const auto P_veff = gw_spec_veff.P();
    
    std::map<std::string, std::string> opts_bag, opts_munu, opts_veff;
    opts_bag["label"] = "Bag";
    opts_bag["color"] = "red";
    opts_bag["linestyle"] = "--";

    opts_munu["label"] = "mu nu";
    opts_munu["color"] = "black";
    opts_munu["linestyle"] = "-.";

    opts_veff["label"] = "Veff";
    opts_veff["color"] = "blue";
    opts_veff["linestyle"] = "-";

    plt::figure_size(800, 600);
    // Set log scaling manually
    PyRun_SimpleString("import matplotlib.pyplot as plt\n"
                    "plt.xscale('log')\n"
                    "plt.yscale('log')");

    plt::plot(freq_bag, P_bag, opts_bag);
    plt::plot(freq_munu, P_munu, opts_munu);
    plt::plot(freq_veff, P_veff, opts_veff);
    plt::xlabel("K=kRs");
    plt::ylabel("Omega_GW(K)");
    plt::xlim(f_min, f_max);
    plt::grid(true);
    plt::legend();    

    plt::save(filename);
    std::cout << "GW spectrum saved to " << filename << "\n";
}
#endif

} // namespace Spectrum