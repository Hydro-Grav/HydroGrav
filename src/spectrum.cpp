// spectrum.cpp
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

#include "ap.h"
#include "interpolation.h"
#include "specialfunctions.h"

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
*/

namespace Spectrum {

/***** PowerSpec class *****/

// Define ctors
PowerSpec::PowerSpec(const std::vector<double>& k_vals, std::vector<double>& P_vals, const PhaseTransition::PTParams& params)
    : data_(Spectrum{k_vals, P_vals}),
      params_(params) {
        if (k_vals.size() != P_vals.size()) {
            throw std::invalid_argument("PowerSpec: k and P vectors must be the same size!");
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
    file << "k,P\n";

    const auto k_vals = data_.first;
    const auto P_vals = data_.second;
    for (size_t i = 0; i < k_vals.size(); ++i) {
        file << k_vals[i] << "," << P_vals[i] << "\n";
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
    plt::suptitle("vw = " + to_string_with_precision(params_.vw()) + ", alN = " + to_string_with_precision(params_.alN()));
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

    return PowerSpec(spec.K(), scaled_P, spec.params());
}

PowerSpec operator*(double scalar, const PowerSpec& spec) {
    return spec * scalar;
}

PowerSpec& PowerSpec::operator*=(double scalar) {
    for (auto& p : data_.second) {
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

    for (auto& p : data_.second) {
        p /= scalar;
    }
    return *this;
}
/***************************/

/*** GW power spectrum ***/
// add option for inputing pRs_vals, Ttilde_vals and z_vals?
PowerSpec GWSpec(const std::vector<double>& kRs_vals, const PhaseTransition::PTParams& params) {
    const Hydrodynamics::FluidProfile profile(params); // generate fluid profile
    const auto prefac = gw_prefac(kRs_vals, profile); // prefactor

    const auto nk = kRs_vals.size();

    const auto pRs_vals = logspace(1e-2, 1e+3, 1000); // P = p*Rs
    const auto np = pRs_vals.size();

    std::vector<double> pRs2_vals(np), p_vals(np); // keep here otherwise have to calculate for each k
    for (size_t i = 0; i < np; i++) {
        const auto pRs = pRs_vals[i];
        pRs2_vals[i] = pRs * pRs;
    }

    const auto z_vals = linspace(-1.0, 1.0, 1000); // logspace gives nan over this domain
    const auto nz = z_vals.size();    

    /********** precompute normalised kinetic spectrum **********/
    // NOTE: this is currently a bit buggy and domain of interpolating function requires fine tuning depending on range of k,p,z
    // zetaKin(ptRs) can't be precomputed since ptRs = ptRs(k,p,z) -> use interpolator (much faster than constructing PowerSpec objects inside loops)

    // calc temp ptRs vals for interpolating func
    const auto pRs_max = pRs_vals.back();
    const auto kRs_max = kRs_vals.back();
    const auto ptRs_max = 10.0 * (kRs_max + pRs_max); // max of pt=sqrt(k^2-2kpz+p^2)
    const auto ptRs_min = 1e-7; // not sure how to choose best min val - update later

    const auto ptRs_vals_tmp = logspace(ptRs_min, ptRs_max, 2.0*np);

    const auto zk_pRs_spec = zetaKin(pRs_vals, profile);
    const auto zk_pRs_vals = zk_pRs_spec.P(); // store zetaKin(pRs) vals (quicker than calling interpolator)

    const auto zk_ptRs_spec = zetaKin(ptRs_vals_tmp, profile);    
    const auto zk_ptRs_interp = zk_ptRs_spec.interpolate(); // interpolating function for zetaKin(ptRs)
    /************************************************************/

    std::cout << "Calculating gravitational wave power spectrum...\n";

    // precompute dlt
    // const auto delta = dlt_SSM(k_vals, p_vals, z_vals, params);
    const auto delta = dlt_SSM2(kRs_vals, pRs_vals, z_vals, params);
    std::vector<double> GW_P_vals(nk);

    #pragma omp parallel for schedule(static)
    for (size_t kk = 0; kk < nk; kk++ ) {
        const auto kRs = kRs_vals[kk];
        const auto kRs3 = kRs * kRs * kRs;

        std::vector<std::vector<double>> integrand(np, std::vector<double>(nz));
        for (size_t pp = 0; pp < np; pp++) {
            const auto pRs = pRs_vals[pp];
            const auto zk_pRs_fac = kRs3 * zk_pRs_vals[pp] * pRs2_vals[pp]; // kRs^3 * zetaKin(pRs) * pRs^2

            for (size_t zz = 0; zz < nz; zz++) {
                const auto z = z_vals[zz];
                const auto ptRs = ptilde(kRs, pRs, z);

                if (ptRs == 0.0) { // careful! need to check this converges properly for pt=0!
                    integrand[pp][zz] = 0.0;
                    continue;
                }

                const auto z_fac = 1.0 - z;
                const auto z_fac2 = z_fac * z_fac;
                const auto ptRs4_inv = 1.0 / (ptRs * ptRs * ptRs * ptRs);
                
                // const auto dlt = delta[kk][pp][zz];
                const auto dlt = delta[kk * np * nz + pp * nz + zz];

                // integrand[pp][zz] = (ptRs != 0.0) ? z_fac2 * ptRs4_inv * zk_pRs_fac * zk_ptRs_interp(ptRs) * dlt : 0.0;
                integrand[pp][zz] = z_fac2 * ptRs4_inv * zk_pRs_fac * zk_ptRs_interp(ptRs) * dlt;
            }
        }

        GW_P_vals[kk] = prefac * simpson_2d_integrate(pRs_vals, z_vals, integrand);
        // GW_P_vals[kk] = simpson_2d_integrate(pRs_vals, z_vals, integrand);
    }

    std::cout << "Gravitational power spectrum constructed!\n";

    return PowerSpec(kRs_vals, GW_P_vals, params);
}
/***************************/

/*** dlt spectrum ***/
double ptilde(double k, double p, double z) {
    const auto arg = k*k - 2.0 * k * p * z + p*p;

    if (std::abs(arg) < 1e-10)
        return 0.0; // avoids numerical precision issues giving arg < 0

    return std::sqrt(arg);
}

double ff(double tau_m, double kcs) {
    // kcs = k*cs -> ff called this way to make dlt faster
    return std::cos(kcs * tau_m); // for SSM -> NEED TO UPDATE THIS
}

double dtau_fin(double tau_fin, double tau_s) {
    return tau_fin - tau_s;
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

// same as dlt_SSM but with flattened index
std::vector<double> dlt_SSM2(const std::vector<double>& kRs_vals, const std::vector<double>& pRs_vals, const std::vector<double>& z_vals, const PhaseTransition::PTParams& params) {
    /***************************** CLOCK ******************************/
    const auto ti = std::chrono::high_resolution_clock::now();
    /******************************************************************/

    const auto cs = std::sqrt(params.cpsq());
    const auto tau_s = params.tau_s();
    const auto tau_fin = params.tau_fin();
    const auto Rs_inv = 1.0 / params.Rs();

    const auto nk = kRs_vals.size();
    const auto np = pRs_vals.size();
    const auto nz = z_vals.size();

    // reserve memory for integration
    std::vector<double> result(nk * np * nz);
    constexpr std::array<double,2> sum_vals = {-1.0, 1.0};

    #pragma omp parallel
    {
        #pragma omp for collapse(3) schedule(static)
        for (size_t kk = 0; kk < nk; kk++)
        for (size_t pp = 0; pp < np; pp++)
        for (size_t zz = 0; zz < nz; zz++) {
            const auto k = kRs_vals[kk] * Rs_inv;
            const auto p = pRs_vals[pp] * Rs_inv;
            const auto z = z_vals[zz];

            const auto pt = ptilde(k, p, z);
            // const auto pt = std::sqrt(k*k - 2.0 * k * p * z + p*p);
            auto dlt_temp = 0.0;

            // regular for loop
            for (int i = 0; i < 2; i++) { // loop over m
                const auto m = sum_vals[i];
                const auto pmn_1 = (p + m * pt) * cs;
                for (int j = 0; j < 2; j++) { // loop over n
                    const auto n = sum_vals[j];
                    const auto pmn = pmn_1 + n * k;

                    const auto x1 = pmn * tau_fin;
                    const auto x2 = pmn * tau_s;

                    double Si_x1, Ci_x1, Si_x2, Ci_x2;
                    sici(x1, Si_x1, Ci_x1);
                    sici(x2, Si_x2, Ci_x2);
                    // alglib::sinecosineintegrals(x1, Si_x1, Ci_x1);
                    // alglib::sinecosineintegrals(x2, Si_x2, Ci_x2);

                    // Im(Si(x))=0 for real x
                    // Im(Ci(x))=pi (x<0), 0 (x>0)
                    // Taking difference dCi -> imaginary part cancels since sign of x1, x2 always the same
                    const auto dSi = Si_x1 - Si_x2;
                    const auto dCi = Ci_x1 - Ci_x2;
                    // const auto dSi_dCi = dSiCi(x1, x2, 1000);
                    // const auto dSi = dSi_dCi[0];
                    // const auto dCi = dSi_dCi[1];

                    dlt_temp += 0.25 * (dCi * dCi + dSi * dSi);
                }
            }

            result[kk * np * nz + pp * nz + zz] = dlt_temp;
        }
    }

    /***************************** CLOCK ******************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer (dlt_SSM): " << duration.count() << " s" << std::endl;
    /******************************************************************/

    return result;
}

std::vector<std::vector<std::vector<double>>> dlt_SSM(const std::vector<double>& k_vals, const std::vector<double>& p_vals, const std::vector<double>& z_vals, const PhaseTransition::PTParams& params) {
    /***************************** CLOCK ******************************/
    const auto ti = std::chrono::high_resolution_clock::now();
    /******************************************************************/

    const auto cs = std::sqrt(params.cpsq());
    const auto tau_s = params.tau_s();
    const auto tau_fin = params.tau_fin();

    const auto nk = k_vals.size();
    const auto np = p_vals.size();
    const auto nz = z_vals.size();

    // reserve memory for integration
    std::vector<std::vector<std::vector<double>>> result(nk, std::vector<std::vector<double>>(np, std::vector<double>(nz)));
    constexpr std::array<double,2> sum_vals = {-1.0, 1.0};

    #pragma omp parallel
    {
        #pragma omp for collapse(3) schedule(static)
        for (size_t kk = 0; kk < nk; kk++)
        for (size_t pp = 0; pp < np; pp++)
        for (size_t zz = 0; zz < nz; zz++) {
            const auto k = k_vals[kk];
            const auto p = p_vals[pp];
            const auto z = z_vals[zz];

            const auto pt = ptilde(k, p, z);
            // const auto pt = std::sqrt(k*k - 2.0 * k * p * z + p*p);
            auto dlt_temp = 0.0;

            // regular for loop
            for (int i = 0; i < 2; i++) { // loop over m
                const auto m = sum_vals[i];
                const auto pmn_1 = (p + m * pt) * cs;
                for (int j = 0; j < 2; j++) { // loop over n
                    const auto n = sum_vals[j];
                    const auto pmn = pmn_1 + n * k;

                    const auto x1 = pmn * tau_fin;
                    const auto x2 = pmn * tau_s;

                    double Si_x1, Ci_x1, Si_x2, Ci_x2;
                    alglib::sinecosineintegrals(x1, Si_x1, Ci_x1);
                    alglib::sinecosineintegrals(x2, Si_x2, Ci_x2);

                    // Im(Si(x))=0 for real x
                    // Im(Ci(x))=pi (x<0), 0 (x>0)
                    // Taking difference dCi -> imaginary part cancels since sign of x1, x2 always the same
                    const auto dSi = Si_x1 - Si_x2;
                    const auto dCi = Ci_x1 - Ci_x2;

                    dlt_temp += 0.25 * (dCi * dCi + dSi * dSi);
                }
            }

            // result[kk][pp][zz] = dlt_temp[thread_id];
            result[kk][pp][zz] = dlt_temp;
        }
    }

    /***************************** CLOCK ******************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer (dlt_SSM): " << duration.count() << " s" << std::endl;
    /******************************************************************/

    return result;
}
/***************************/

/*** Kinetic spectrum ***/
// avoids duplicating fluid profile in GWSpec
PowerSpec Ekin(const std::vector<double>& kRs_vals, const PhaseTransition::PTParams& params) {
    return Ekin(kRs_vals, Hydrodynamics::FluidProfile(params));
}

PowerSpec Ekin(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& prof) {
    const auto beta = prof.params().beta();
    const auto Rs = prof.params().Rs();
    const auto nuc_type = prof.params().nuc_type();

    auto lt_dist = Hydrodynamics::lifetime_dist_func(nuc_type);

    const auto nk = kRs_vals.size();
    std::vector<double> P_vals(nk);
    // std::vector<double> P_vals;


    // define Ttilde from chi = Ttilde * k / beta (makes calling Apsq simpler)
    // using K = k * Rs below
    /*
    NOTE:
    - Ap_sq = inf at 0
    - chi_vals = logspace(1e-3, 3000, 5000) gives good convergence
    */
    const auto chi_vals = logspace(1e-3, 3000, 5000); // bad to hard code?
    const auto n = chi_vals.size();

    const auto Apsq = Hydrodynamics::Ap_sq(chi_vals, prof);

    const auto fac1 = beta * Rs * Rs / (2.0 * M_PI * M_PI);

    std::vector<std::vector<double>> integrands(omp_get_max_threads(), std::vector<double>(n));
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::vector<double>& integrand = integrands[tid]; // local thread copy of integrand

        #pragma omp for
        for (size_t kk = 0; kk < nk; kk++) {
            const auto kRs = kRs_vals[kk];
            const auto kRs_inv = 1.0 / kRs;

            const auto fac2 = fac1 * power(kRs_inv, 5);
            const auto fac3 = beta * Rs * kRs_inv;

            for (size_t i = 0; i < n; i++) {
                const auto chi = chi_vals[i];
                integrand[i] = fac2 * lt_dist(fac3 * chi) * power(chi, 6) * Apsq[i];
            }

            P_vals[kk] = simpson_integrate(chi_vals, integrand);
        }
    }

    return PowerSpec(kRs_vals, P_vals, prof.params());
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

// not finished
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
    const auto un = params.un();

    const auto Ek = Ekin(kRs_vals, profile);

    return gw_prefac(Ek.max(), params.Rs(), params.wNeN_rat(), un.T0(), un.Ts(), un.H0(), un.Hs(), un.g0(), un.gs());
}

} // namespace Spectrum