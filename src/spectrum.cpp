// spectrum.cpp
#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <shared_mutex>

#include "ap.h"
#include "interpolation.h"
#include "specialfunctions.h"

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

#include "maths_ops.hpp"
#include "phasetransition.hpp"
#include "hydrodynamics.hpp"
#include "spectrum.hpp"
#include "physics.hpp"
#include "sici.hpp"

#ifdef ENABLE_MATPLOTLIB
#include "matplotlibcpp.h"
#endif

/*
TO DO:
- update prefac to allow for non-bag model
- update prefac to do actual calculation of TGW, OmegaK_KK
- remove instances of std::pow when possible - it is slow
- change throw exception for P() and K() so that it uses P() and K() when wrong one is called
- update Ekin to pass in Profile class (or maybe just PTParams?)
- implement adaptive step-size in Ekin integration (and dlt later too)
- change input of GWSpec2 to frequencies then convert to kRs internally
*/

namespace Spectrum {

/***** PowerSpec class *****/

// Define ctors
// PowerSpec::PowerSpec(const std::vector<double>& K_vals, std::vector<double>& P_vals, const PhaseTransition::PTParams& params)
PowerSpec::PowerSpec(const std::vector<double>& K_vals, std::vector<double>& P_vals, const Hydrodynamics::FluidProfile& profile)
    : freq_vals_(), K_vals_(K_vals), P_vals_(P_vals),
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

// Public functions
double PowerSpec::max() const {
    const auto &Pv = P();
    return *std::max_element(Pv.begin(), Pv.end());
}

void PowerSpec::write(const std::string& filename) const {
    std::cout << "Writing power spectrum to disk... ";
    std::ofstream file(filename);
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
    namespace plt = matplotlibcpp;

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

    // kinetic_spectrum.write("zetaKin_ptRs.csv");

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
PowerSpec GWSpec(const std::vector<double>& kRs_vals, const PhaseTransition::PTParams& params) {

    const auto ti = std::chrono::high_resolution_clock::now();

    const auto cs = std::sqrt(params.cpsq());
    const auto tau_s = params.tau_s();
    const auto tau_fin = params.tau_fin();
    const auto Rs_inv = 1.0 / params.Rs();

    const Hydrodynamics::FluidProfile profile(params);

    const auto nk = kRs_vals.size();

    const double pRs_minimum = 1e-3;
    const double pRs_maximum = 1e+3;
    const auto n_pRs = 500;

    const auto pRs_vals = logspace(pRs_minimum, pRs_maximum, n_pRs); // P = p*Rs

    const auto kinetic_spectrum_spline_lower_bound = 0.99 * find_min_pt(kRs_vals, pRs_vals);
    const auto kinetic_spectrum_spline_upper_bound = 1.01 * ptilde(kRs_vals.back(), pRs_vals.back(), -1.0);

    const auto kinetic_spectrum_K_values = logspace(kinetic_spectrum_spline_lower_bound, kinetic_spectrum_spline_upper_bound, 2*n_pRs);

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

                const double z_result = boost::math::quadrature::gauss_kronrod<double, 31>::integrate(z_integrand, -1.0, 1.0, 5, 1e-6);

                return pRs * zk_pRs_fac * z_result;
            };

            // double pRs_result = boost::math::quadrature::trapezoidal(pRs_integrand, log(pRs_minimum), log(pRs_maximum), 1e-6);
            double pRs_result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(pRs_integrand, log(pRs_minimum), log(pRs_maximum), 5, 1e-6);
            GW_P_vals[kk] = prefac * kRs3 * pRs_result;
        }
    }

    std::cout << "Gravitational power spectrum constructed!\n";

    /***************************** CLOCK ******************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer (GWSpec2): " << duration.count() << " s" << std::endl;
    /******************************************************************/

    return PowerSpec(kRs_vals, GW_P_vals, profile);
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

// tau_vals is the same for each call of dlt() -> redundance when calling dlt() in a loop since it recalculates tau_m each time
// best to store dlt as a nk x np x npt tensor to avoid this
// NOTE: dlt takes in k, NOT K=kRs
// generic dlt (not just SSM) - very slow!
std::vector<std::vector<std::vector<double>>> dlt(const int nt, const std::vector<double>& k_vals, const std::vector<double>& p_vals, const std::vector<double>& z_vals, const PhaseTransition::PTParams& params) {
    /***************************** CLOCK ******************************/
    const auto ti = std::chrono::high_resolution_clock::now();
    /******************************************************************/

    const auto cs = std::sqrt(params.cpsq());

    const auto tau_s = params.tau_s();
    const auto tau_fin = params.tau_fin();

    // integrand becomes very large for small tau -> use logspace for accuracy of integration
    // const auto nt = 50;
    const auto tau_vals = logspace(tau_s, tau_fin, nt);
    const auto ntsq = nt * nt;

    const auto nk = k_vals.size();
    const auto np = p_vals.size();
    const auto nz = z_vals.size();

    // store tau_m = tau2 - tau1 and tau_sq_inv = 1/(tau1 * tau2) values
    // avoids repeated calculation in loops over k, p, z (much quicker!)
    std::vector<double> tau_m(ntsq);
    std::vector<double> tau_sq_inv(ntsq);
    #pragma omp parallel for
    for (int idx = 0; idx < ntsq; idx++) {
        const int i = idx / nt; // idx = i * nt + j;
        const int j = idx % nt;

        const auto tau1 = tau_vals[i];
        const auto tau2 = tau_vals[j];

        tau_m[idx] = tau2 - tau1;
        tau_sq_inv[idx] = 1.0 / (tau1 * tau2);
    }

    // fill ff (reduces redundancy)
    std::vector<std::vector<double>> ff1_cache(ntsq, std::vector<double>(np));
    std::vector<std::vector<double>> ff3_cache(ntsq, std::vector<double>(nk));
    #pragma omp parallel for
    for (int idx = 0; idx < ntsq; idx++) {
        const auto tau_minus = tau_m[idx];
        // fill ff3
        for (size_t kk = 0; kk < nk; kk++) {
            const auto k = k_vals[kk];
            ff3_cache[idx][kk] = std::cos(k * tau_minus);
        }
        // fill ff1
        for (size_t pp = 0; pp < np; pp++) {
            const auto p = p_vals[pp];
            ff1_cache[idx][pp] = ff(tau_minus, p*cs); // redundancy computing p*cs here - change?
        }
    }

    // reserve memory for integration
    std::vector<std::vector<std::vector<double>>> result(nk, std::vector<std::vector<double>>(np, std::vector<double>(nz)));
    std::vector<std::vector<double>> integrands(omp_get_max_threads(), std::vector<double>(ntsq));

    // precompute weights for integration
    const auto weights = precompute_simpson_weights_2d(tau_vals, tau_vals);
    const auto Ax_weights = weights.Ax_weights;
    const auto Ay_weights = weights.Ay_weights;
    const auto dx = weights.dx;
    const auto dy = weights.dy;

    // integration routine
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::vector<double>& integrand = integrands[tid]; // local thread copy of integrand

        #pragma omp for
        for (size_t idx = 0; idx < nk * np * nz; idx++) {
            const int kk = idx / (np * nz);
            const int pp = (idx / nz) % np;
            const int zz = idx % nz;
        // #pragma omp for collapse(3) schedule(dynamic)
        // for (int kk = 0; kk < nk; kk++)
        // for (int pp = 0; pp < np; pp++)
        // for (int zz = 0; zz < nz; zz++) {
            const auto k = k_vals[kk];
            const auto p = p_vals[pp];
            const auto z = z_vals[zz];

            if (std::isnan(z)) {
                throw std::runtime_error("bad z");
            }

            const auto pt = ptilde(k, p, z); // collapsing loops a lot quicker than breaking up ptilde calc so redundancy here is ok!
            const auto ptcs = pt * cs;

            // integration routine
            for (int i = 0; i < ntsq; i++) {
                const auto tau_minus = tau_m[i];
                const auto ff3 = ff3_cache[i][kk];
                const auto ff1 = ff1_cache[i][pp];
                const auto ff2 = ff(tau_minus, ptcs);

                integrand[i] = ff1 * ff2 * ff3 * tau_sq_inv[i];
            }

            // result[kk][pp][zz] = simpson_2d_nonuniform_flat(tau_vals, tau_vals, integrand);
            result[kk][pp][zz] = simpson_2d_nonuniform_flat_weighted(tau_vals, tau_vals, integrand, Ax_weights, Ay_weights, dx, dy);
        }
    }

    /***************************** CLOCK ******************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer (dlt): " << duration.count() << " s" << std::endl;
    /******************************************************************/

    return result;
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

    const auto chi_vals = logspace(1e-3, 5e3, 1000);
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

        const double log_chi_min = std::log(1e-3);
        const double log_chi_max = std::log(5000.0);
        
        double error;
        const double tol = 1e-6;
        const int max_iter = 5;

        P_vals[kk] = fac2 * boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, log_chi_min, log_chi_max, max_iter, tol, &error);
    }

    return PowerSpec(kRs_vals, P_vals, prof);
}

PowerSpec norm_spec(const PowerSpec& spec) {
    const auto spec_max = spec.max();
    if (spec_max == 0.0) {
        throw std::runtime_error("Division by zero in norm_spec from spec.max() = 0");
    } else if (isnan(spec_max)) {
        throw std::runtime_error("In norm_spec: spec.max() = nan");
    }

    const auto zk = spec / spec_max;
    if (abs(zk.max() - 1.0) > 1e-15) {
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

double gw_prefac(double Ekin_max, double Rs, double wNeN_rat, double T0, double Ts, double H0, double Hs, double g0, double gs) {
    // Transfer function (redshift of spectrum - eq 13 arXiv:2308.12943)
    const auto g0gs_rat = g0 / gs;
    const auto TH_rat = (T0 * T0 / H0) / (Ts * Ts / Hs);
    const auto TGW = std::pow(g0gs_rat, 4./3.) * TH_rat * TH_rat;

    // Normalised kinetic energy density OmegaK / KK (eq 42 arXiv:2308.12943)
    // OmegaK = total kinetic energy density, KK = critical energy density
    const auto OmegaK_KK = Ekin_max / Rs;
    
    return 3.0 * wNeN_rat * wNeN_rat * TGW * OmegaK_KK * OmegaK_KK;
}

double gw_prefac(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& profile) {
    const auto params = profile.params();
    const auto un = params->un();

    const auto Ek = Ekin(kRs_vals, profile);

    return gw_prefac(Ek.max(), params->Rs(), params->wNeN_rat(), un.T0(), un.Ts(), un.H0(), un.Hs(), un.g0(), un.gs());
}

} // namespace Spectrum