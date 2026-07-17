// profile.cpp
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <optional>
#include <chrono>

#include "ap.h"
#include "interpolation.h"

#include "profile.hpp"
#include "phasetransition.hpp"
#include "physics.hpp"
#include "maths.hpp"
#include "config.hpp"

#ifdef ENABLE_MATPLOTLIB
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

/*
TO DO:
- update generate_streamplot_data() to include points of interest (fixed pts, detonation/deflag/hybrid regions)
- fix calc of w profile in generate_streamplot_data()
- write summary of detonation/deflag/hybrid at top of ctor
- allocate memory to y_sol_tmp etc (inside solve) -> get rid of .push_back() and change to indexing
- update get_mode for veff
- update "continue" to "break" in find_shock_idx? (if (xi_sh <= vw_) continue;) -> i think this will be faster since xi runs backwards and after the first 'continue' the rest will be 'continue' too
*/

namespace Hydrodynamics { // calculate bubble profile

double mu(double xi, double v) {
    return (xi - v) / (1.0 - xi * v);
}

/***************************** Fluid equations of motion *****************************/
/* EoM for a perfect fluid comes from \partial_{\mu} T^{mu nu} = 0:                  */
/*    2v/xi = gamma^2 (1 - v*xi)(mu^2/cs^2 - 1) v'                                   */
/*    w'/w = gamma^2 * mu (1 + 1/cs^2) v'                                            */
/* where xi=r/t, gamma = 1/sqrt(1-v^2), v'=dv/dxi, w'=dw/dxi and mu=(xi-v)/(1-xi*v). */

// these should be vUF - change for consistentsy
double dxidv(double xi, double v, const double csq) {
    const auto mu_val = mu(xi, abs(v));
    if (v < 1e-15) return 0.0; // prevents divergence at (xi,v)=(csq,0)
    return gammaSq(v) * (1.0 - v * xi) * (mu_val * mu_val / csq - 1.0) * xi / (2.0 * v);
}

double dwdv(double xi, double v, double w, const double csq) {
    return w * gammaSq(v) * mu(xi, v) * (1.0 + 1.0 / csq);
}

double dTdv(double xi, double v, double T) {
    return T * gammaSq(v) * mu(xi, v);
}

state_type dydv_vec(double v, const state_type& y, double vw, double cmsq, double cpsq) {
    const auto xi = y[0];
    const auto w = y[1];
    const auto T = y[2];
    const auto csq = (xi <= vw) ? cmsq : cpsq;

    return { dxidv(xi, v, csq), dwdv(xi, v, w, csq), dTdv(xi, v, T) };
}
/*************************************************************************************/

/****************************** FluidProfile class ******************************/

    /********************** Notation for fluid parameters **********************/
    /* xi_sh = position of shock wave                                          */
    /*                                                                         */
    /* bubble wall frame:                                                      */
    /*   vm, vp (velocity of fluid behind (m) and in front (p) of bubble wall) */
    /*   wm, wp (enthalpy of fluid behind (m) and in front (p) of bubble wall) */
    /*                                                                         */
    /* shock frame:                                                            */
    /*   v1, v2 (velocity of fluid behind (1) and in front of (2) shock)       */
    /*                                                                         */
    /* centre of bubble/centre of shock frame (universe frame):                */
    /*   As above, but ending with 'UF'                                        */
    /***************************************************************************/

FluidProfile::FluidProfile(const PhaseTransition::PTParams& params, const size_t n, const bool dev_log)
    : params_(&params),
      cpsq_(params.cpsq()), cmsq_(params.cmsq()),
      vw_(params.vw()), alN_(params.alN()),
      alp_min_(std::numeric_limits<double>::quiet_NaN()), 
      alp_max_(std::numeric_limits<double>::quiet_NaN()),
      mode_(),
      xi_vals_(), v_vals_(), w_vals_(), T_vals_(), la_vals_(),
      shock_flag_(true),
      dev_log_(dev_log)
    {
        std::vector<prof_type> profiles;

        // define hydrodynamic mode for bag model
        const auto mode_bag = get_mode_bag(vw_, cmsq_, alN_);

        // check alN large enough for shock (deflag/hybrid only)
        if (mode_ == 0 || mode_ == 1) {
            const auto alp_minmax = get_alp_minmax(vw_);
            alp_min_ = alp_minmax[0];
            alp_max_ = alp_minmax[1];

            // alN > alp > alp_min (can't properly constrain from above since we need alp)
            if (alN_ <= alp_min_) throw std::runtime_error("alN too small for shock!");
        }

        // calculate fluid profiles v(xi), w(xi), la(xi)
        switch (params.eos()) {
            case PhaseTransition::PTParams::ModelType::Bag:
                if (cpsq_ == 1.0 / 3.0 && cmsq_ == cpsq_) { // bag
                    std::cout << "Constructing fluid profile using bag equation of state (cp=cm=1/sqrt(3))...\n";
                } else { // mu nu
                    std::cout << "Constructing fluid profile using mu-nu (modified bag) equation of state (cp=" << std::sqrt(cpsq_) << ", cm=" << std::sqrt(cmsq_) << ")...\n";
                }

                bag_params_ = &dynamic_cast<const PhaseTransition::PTParams_Bag&>(params);
                mode_ = mode_bag;
                profiles = solve_profile(n);
                break;
            case PhaseTransition::PTParams::ModelType::Veff:
                std::cout << "Constructing fluid profile using generic equation of state from Veff...\n";

                veff_params_ = &dynamic_cast<const PhaseTransition::PTParams_Veff&>(params);
                // mode_ = get_mode_veff(vw_, cmsq_);
                std::cout << "Determining hydrodynamic mode from EoS... ";
                try {
                    mode_ = get_mode_veff(vw_, cmsq_);
                    std::cout << "Hydrodynamic mode determined successfully!\n";
                } catch (std::exception& e) {
                    std::cerr << "WARNING: " << e.what() << " Using get_mode_bag() to estimate hydrodynamic mode instead!\n";
                    mode_ = mode_bag;
                }

                // check if hydrodynamic modes agree between bag and veff
                if (mode_ != mode_bag) {
                    std::cerr << "WARNING: hydrodynamic modes do not agree between mu-nu (mode=" << mode_bag << ") and Veff (mode=" << mode_ << ") EoS! Using hydrodynamic mode from Veff.\n";
                }

                profiles = solve_profile_veff(n);
                break;
            default:
                throw std::runtime_error("Unknown EOS type");
        }

        xi_vals_ = profiles[0];
        v_vals_ = profiles[1];
        w_vals_ = profiles[2];
        T_vals_ = profiles[3];
        la_vals_ = profiles[4];

        std::cout << "Fluid profile constructed!\n";
    }

// Public functions
std::string FluidProfile::mode_str() const {
    if (mode_ == 0) return "deflagration";
    else if (mode_ == 1) return "hybrid";
    else if (mode_ == 2) return "detonation";
    else return "unknown";
}

void FluidProfile::write(const std::string& filename) const {
    std::cout << "Writing fluid profile to disk... ";

    std::ofstream file(filename);
    file << "xi,v,w,T,la\n";

    for (size_t i = 0; i < xi_vals_.size(); ++i) {
        file << xi_vals_[i] << "," << v_vals_[i] << "," << w_vals_[i] << "," << T_vals_[i] << "," << la_vals_[i] << "\n";
    }
    file.close();

    std::cout << "Saved to " << filename << "!\n";

    return;
}

#ifdef ENABLE_MATPLOTLIB
void FluidProfile::plot(const std::string& filename, double xi_min, double xi_max) const {
    std::cout << "Generating fluid profile plot... ";

    plt::figure_size(2400, 800);

    // v(xi)
    plt::subplot2grid(2, 2, 0, 0);
    plt::plot(xi_vals_, v_vals_);
    plt::xlabel("xi");
    plt::ylabel("v(xi)");
    plt::xlim(xi_min, xi_max);
    plt::grid(true);

    // w(xi)
    plt::subplot2grid(2, 2, 0, 1);
    plt::plot(xi_vals_, w_vals_);
    plt::xlabel("xi");
    plt::ylabel("w(xi)");
    plt::xlim(xi_min, xi_max);
    plt::grid(true);

    // T(xi)
    plt::subplot2grid(2, 2, 1, 0);
    plt::plot(xi_vals_, T_vals_);
    plt::xlabel("xi");
    plt::ylabel("T(xi)");
    plt::xlim(xi_min, xi_max);
    plt::grid(true);

    // la(xi)
    plt::subplot2grid(2, 2, 1, 1);
    plt::plot(xi_vals_, la_vals_);
    plt::xlabel("xi");
    plt::ylabel("la(xi)");
    plt::xlim(xi_min, xi_max);
    plt::grid(true);

    plt::suptitle("vw = " + to_string_with_precision(vw_) + ", alpha = " + to_string_with_precision(alN_));
    plt::save(filename);

    std::cout << "Saved to '" << filename << "'." << std::endl;

    return;
}
#endif

// Private functions

/**************************************** Bag EoS ****************************************/
int FluidProfile::get_mode_bag(double vw, double cmsq, double alN) const {
    const auto vwsq = vw * vw;

    if (vwsq < cmsq) return 0; // deflagration
    if (vw < vJ_det(alN)) return 1; // hybrid
    return 2; // detonation
}

double FluidProfile::vJ_det(double alp) const {
    return vp_from_matching(std::sqrt(cmsq_), alp); // vJ(alp) = vp(|vm|=cm, alp)
}

std::array<double, 2> FluidProfile::get_alp_minmax(double vw) const {
    const auto cm = std::sqrt(cmsq_);

    const auto vm = std::min(cm, vw); // vw for deflag (vw < cm), cm for hybrid (cp < vw)
    const auto vp_min = 0.0;
    const auto vp_max = vm; // |v+| < |v-|
    
    const auto al_max = get_alp_wall(vp_min, vm);
    const auto al_min = get_alp_wall(vp_max, vm);

    if (al_min < 0.0) return {0.0, al_max}; // numerical precision issue

    return {al_min, al_max};
}

double FluidProfile::get_alp_wall(double vp, double vm) const {
    return gammaSq(vp) * (vp * vp / cmsq_ - vp * vm / cmsq_ - vp / vm + 1.0) / 3.0;
}

double FluidProfile::vp_from_matching(double vm, double alp) const { // vp(vm,alp) from matching eqs
    const auto sgn = 1.0; // for detonations
    // const auto sgn = (abs(vm) > std::sqrt(cmsq_)) ? 1.0 : -1.0;

    const auto fac1 = 1.0 / (2.0 * (1.0 / (3.0 * cmsq_) + alp));
    const auto fac2 = 1.0 / (3.0 * vm);
    const auto fac3 = fac2 - vm / (3.0 * cmsq_);

    return fac1 * (fac2 + vm / (3.0 * cmsq_) + sgn * std::sqrt(fac3 * fac3 + 4.0 * alp * alp + 4.0 * (1.0/cmsq_ - 1.0) * alp / 3.0));
}

double FluidProfile::vm_from_matching(double vp, double alp) const { // inverse of vp(vm,alp)
    // const auto sgn = 1.0;
    const auto sgn = (abs(vp) < std::sqrt(cmsq_)) ? 1.0 : -1.0;
    const auto fac = vp + cmsq_ * (1.0 - 3.0 * alp * (1.0 - vp * vp)) / vp;

    return 0.5 * (fac + sgn * std::sqrt(fac * fac - 4.0 * cmsq_)); // mu nu
}

double FluidProfile::get_TmTN(double wmwN) const {
    // const auto mu = 1.0 + 1.0 / cpsq_;
    // const auto nu = 1.0 + 1.0 / cmsq_;
    // const auto r = 1.0; // ap/am ratio

    // if (cpsq_ == cmsq_) return std::pow(r * wmwN, 1.0 / mu);

    // const auto fac = (mu / nu) * r * wmwN;
    // return std::pow(fac, 1.0 / nu) * std::pow(bag_params_->TN(), mu / nu - 1.0);

    const auto r = 1.0; // ap/am ratio

    if (cpsq_ == cmsq_) return std::pow(r * wmwN, 0.25);
    return std::pow(r * (1.0 + cpsq_) / (1.0 + cmsq_) * wmwN, 0.25);
}

std::pair<double, double> FluidProfile::wT_from_shock(double xi_sh) const { // w1/wN
    // w1*gammaSq(v1)*v1 = wN*gammaSq(v2)*v2, v2=xi_sh, v1*v2=cpsq
    const auto xi_sh_sq = xi_sh * xi_sh;
    const auto w1wN = (xi_sh_sq - cpsq_ * cpsq_) / (cpsq_ * (1.0 - xi_sh_sq));
    const auto T1TN = std::pow(w1wN, 0.25);

    return {w1wN, T1TN};
}

// deflagrations
// change bracket so vpUF_max = vpUF that gives alp=alN (any larger and alp > alN which is not allowed)
double FluidProfile::find_vpUF(const deriv_func& dydv, const size_t n) const {
    const auto vm = (mode_ == 0) ? vw_ : std::sqrt(cmsq_);

    auto safe_residual = [this, &dydv, vm, n] (double vpUF) {
        try {
            const auto resi = alN_residual(dydv, vpUF, vm, n);
            if (!std::isfinite(resi)) return 1e+5;
            return resi;
        } catch (std::exception& e) {
            return 1e+5;
        }
    };

    // test_alN_residual(dydv, vm, 500);

    // find bracket where residual is defined and minimum lies
    // const auto vpUF_min = 0.0;
    // const auto vpUF_max = 1.0;
    // auto bracket = find_bracket(safe_residual, vpUF_min, vpUF_max);

    // // fallback if find_bracket fails
    // if (!std::isfinite(bracket[0]) || !std::isfinite(bracket[1])) {
    //     bracket = find_bracket(safe_residual, 1e-4, 0.999);
    // }

    // if (!std::isfinite(bracket[0]) || !std::isfinite(bracket[1])) {
    //     throw std::runtime_error("find_vpUF failed: bracket not found!");
    // }

    // find bracket where residual is defined and minimum lies
    const std::array<std::pair<double,double>, 2> vpUF_ranges = {{
        {0.0,  1.0},
        {1e-4, 0.999} // fallback if find_bracket fails across vpUF in [0,1]
    }};

    std::array<double, 2> bracket = {std::numeric_limits<double>::quiet_NaN(),
                                    std::numeric_limits<double>::quiet_NaN()};

    for (const auto& [vpUF_min, vpUF_max] : vpUF_ranges) {
        bracket = find_bracket(safe_residual, vpUF_min, vpUF_max);
        if (std::isfinite(bracket[0]) && std::isfinite(bracket[1])) break;
    }

    if (!std::isfinite(bracket[0]) || !std::isfinite(bracket[1])) {
        throw std::runtime_error("find_vpUF failed: bracket not found!");
    }

    // find vp that minimise residual
    const auto vpUF = golden_section_minimize(safe_residual, bracket[0], bracket[1]);
    if (vpUF < 0.0) throw std::runtime_error("find_vpUF failed: vpUF < 0!");
    if (vpUF >= 1.0) throw std::runtime_error("find_vpUF failed: vpUF > 1!");

    return vpUF;
}

std::tuple<std::vector<double>, std::vector<state_type>, bool> FluidProfile::deflagration_profile(const deriv_func& dydv, double vpUF, const bool final_prof, const size_t n) const {
    // construct profile using IC {vpUF,wpwN,TpTN}={vpUF,1,1}
    auto [v_sol, y_sol, shock_flag] = deflagration_profile_internal(dydv, vpUF, 1.0, 1.0, final_prof, n);
    const auto w_end_val = y_sol.back()[1];
    const auto T_end_val = y_sol.back()[2];

    // read off shock and calculate fluid behind shock
    // const auto xi_sh = std::max(y_sol.back()[0], std::sqrt(cpsq_));
    const auto xi_sh = (shock_flag) ? std::sqrt(cpsq_) : std::max(y_sol.back()[0], std::sqrt(cpsq_));

    const auto [w1wN, T1TN] = wT_from_shock(xi_sh);

    const auto w_fac = w1wN / w_end_val;
    const auto T_fac = T1TN / T_end_val;

    // correctly normalise w/wN, T/TN profiles
    for (size_t i = 0; i < v_sol.size(); i++) {
        y_sol[i][1] *= w_fac; // normalize w profile so w1/wN=1 & scale to correct w1/wN
        y_sol[i][2] *= T_fac;
    }

    // check shock matching condition for final profile
    if (final_prof) {
        const auto r1 = abs(mu(xi_sh, v_sol.back()) * xi_sh - cpsq_);
        const auto r2 = abs(xi_sh/mu(xi_sh, v_sol.back()) - (T1TN*T1TN*T1TN*T1TN + cpsq_) / (cpsq_ * T1TN*T1TN*T1TN*T1TN + 1.0));

        if (r1 > config::sh_resi_tol || r2 > config::sh_resi_tol) {
            std::cerr << "Warning in deflagration_profile: Shock residual above tolerance (" << config::sh_resi_tol << "). R1=" << r1 << ", R2=" << r2 << "\n";
        }
    }

    // return {v_sol, y_sol};
    return std::make_tuple(v_sol, y_sol, shock_flag);
}

double FluidProfile::alN_residual(const deriv_func& dydv, double vpUF, double vm, const size_t n) const {
    const auto y_prof = std::get<1>(deflagration_profile(dydv, vpUF, false, n));
    const auto wpwN = y_prof.front()[1];

    const auto fac = (1.0 / cmsq_ - 1.0 / cpsq_) / (1.0 / cpsq_ + 1) / 3.0;
    const auto alN_calc = wpwN * get_alp_wall(mu(vw_, vpUF), vm) + fac * (wpwN - 1.0);
    // const auto alN_calc = wpwN * get_alp_wall(mu(vw_, vpUF), vm);

    return abs(alN_ - alN_calc);
}

std::tuple<std::vector<double>, std::vector<state_type>, bool> FluidProfile::deflagration_profile_internal(const deriv_func& dydv, double vpUF, double wpwN, double TpTN, const bool final_prof, const size_t n) const {
    // solve EoM - integrate from v(xi)=vpUF -> v(xi)=0
    const state_type yp = {vw_, wpwN, TpTN}; // {xi_w, wpwN, TpTN}
    const auto [v_sol_tmp, y_sol_tmp] = rk4_solver(dydv, vpUF, 1e-10, yp, n);

    // find shock & truncate solution
    const auto find_sh_idx = find_shock_idx(v_sol_tmp, y_sol_tmp, final_prof);
    const auto sh_idx = std::min(find_sh_idx, v_sol_tmp.size()-1);

    // handles edge case where xi_sh = cp (i.e. shock discontinuity dissappears)
    // uses flag to handle shock value rather than integrating from vpUF -> 0.0 since eom diverges at v=0
    const bool shock_flag = (find_sh_idx == v_sol_tmp.size() - 1) ? true : false;

    // re-integrate for final profile
    // avoids insufficient no. points for integrating profile when shock is close to vpUF
    if (final_prof) {
        const auto [v_sol, y_sol] = rk4_solver(dydv, vpUF, v_sol_tmp[sh_idx], yp, n); // integrates from vpUF->v1UF
        return std::make_tuple(v_sol, y_sol, shock_flag);
    }

    const std::vector<double> v_sol(v_sol_tmp.begin(), v_sol_tmp.begin() + sh_idx + 1);
    const std::vector<state_type> y_sol(y_sol_tmp.begin(), y_sol_tmp.begin() + sh_idx + 1);

    // return {v_sol, y_sol};
    return std::make_tuple(v_sol, y_sol, shock_flag);
}

size_t FluidProfile::find_shock_idx(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol, const bool final_prof) const {
    if (final_prof) {
        // test_shock_bag(v_sol, y_sol);
    }

    std::vector<double> resi_vals(v_sol.size());
    int pass_count = 0;

    for (size_t i = 0; i < v_sol.size(); i++) {
        const auto xi_sh = y_sol[i][0];
        if (xi_sh <= vw_ || xi_sh > 1.0) {
            resi_vals[i] = 1.0;
            continue;
        }
        
        const auto v1UF = v_sol[i];
        const auto v1 = mu(xi_sh, abs(v1UF)); // v1=mu(xi_sh, v1UF)
        const auto v2 = xi_sh; // v2=xi_sh

        // shock condition mu(xi_sh, v(xi_sh)) xi_sh = cp^2
        resi_vals[i] = abs(v1 * v2 - cpsq_);
        pass_count++;
    }

    if (pass_count == 0) {
        throw std::runtime_error("find_shock_idx_veff failed (no shock found in fluid profile)!");
    }

    // index where residual is minimised
    const auto it = std::min_element(resi_vals.begin(), resi_vals.end());
    const auto idx = std::distance(resi_vals.begin(), it);

    return idx;
}

// dev
void FluidProfile::test_alN_residual(const deriv_func& dydv, double vm, const size_t n) const {
    std::cout << "Running test for alN residual...\n";

    const auto vpUF_vals = linspace(1e-4, 0.2, n);
    std::vector<double> resi_vals(n);

    for (size_t i = 0; i < n; i++) {
        double resi = std::numeric_limits<double>::quiet_NaN();
        try {
            resi = alN_residual(dydv, vpUF_vals[i], vm);
        } catch (std::exception& e) {
            // std::cout << "Failed for vpUF=" << vpUF_vals[i] << ":" << e.what() << "\n";
        }

        resi_vals[i] = resi;
    }

    const auto it = std::min_element(resi_vals.begin(), resi_vals.end());
    const auto idx = std::distance(resi_vals.begin(), it);
    std::cout << "test_alN_residual: vpUF=" << vpUF_vals[idx] << ", min_resi=" << resi_vals[idx] << "\n";
    
    #ifdef ENABLE_MATPLOTLIB
    plt::figure_size(800, 800);
    plt::plot(vpUF_vals, resi_vals);
    plt::xlabel("vpUF");
    plt::ylabel("residual");
    // plt::xlim(0.6, 0.002);
    // plt::ylim(0.0, 0.1);
    plt::grid(true);
    plt::save("alN_resi.png");
    #endif

    std::cout << "Test complete. alN residual saved to 'alN_resi.png'\n";
}

void FluidProfile::test_shock_bag(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol) const {
    std::cout << "Running test for shock residual...\n";

    std::vector<double> xi_vals, resi_vals, r1_vals, r2_vals;
    for (size_t i = 0; i < v_sol.size(); i++) {
        const auto xi_sh = y_sol[i][0];
        if (xi_sh <= vw_ || xi_sh > 1.0) continue;
        
        const auto v1UF = v_sol[i];
        const auto v1 = mu(xi_sh, abs(v1UF)); // v1=mu(xi_sh, v1UF)
        const auto v2 = xi_sh; // v2=xi_sh
        const auto T1TN = y_sol[i][2];

        // shock condition mu(xi_sh, v(xi_sh)) xi_sh = (p1-pN)/(e1-eN)
        xi_vals.push_back(xi_sh);
        const auto r1 = abs(v1 * v2 - cpsq_);
        const auto r2 = abs(v2/v1 - (T1TN*T1TN*T1TN*T1TN + cpsq_) / (cpsq_ * T1TN*T1TN*T1TN*T1TN + 1.0));
        r1_vals.push_back(r1);
        r2_vals.push_back(r2);

        resi_vals.push_back(r1);
    }

    const auto it = std::min_element(resi_vals.begin(), resi_vals.end());
    const auto idx = std::distance(resi_vals.begin(), it);
    std::cout << "test_shock: xi_sh=" << xi_vals[idx] << ", R1=" << r1_vals[idx] << ", R2=" << r2_vals[idx] << ", min_resi=" << resi_vals[idx] << "\n";

    #ifdef ENABLE_MATPLOTLIB
    plt::figure_size(800, 800);
    plt::plot(xi_vals, resi_vals);
    plt::xlabel("xi_sh");
    plt::ylabel("residual");
    // plt::xlim(1.0, 1.1);
    plt::grid(true);
    plt::save("shock_resi_bag.png");
    #endif

    std::cout << "Test complete. Shock residual saved to 'shock_resi_bag.png'\n";
}

// detonations
std::pair<double, state_type> FluidProfile::get_IC_detonation() const {
    // initial conditions:
    //     v0 = v(xi_w) = vm(UF)
    //     w0 = w(xi_w) = wm/wN
    //     T0 = T(xi_w) = Tm/TN
    
    const auto xi0 = vw_;
    const auto vp = -vw_;
    const auto wpwN = 1.0; // w+ = wN
    const auto TpTN = 1.0;

    const auto alp = alN_; // alpha_+ = alpha_N

    const auto vm = vm_from_matching(vp, alp);
    if (abs(vm) >= vw_) throw std::runtime_error("Detonation IC failed: vm<vw required!");

    const auto vmUF = mu(vw_, abs(vm));
    if (abs(vmUF) >= 1.0) throw std::runtime_error("Detonation IC failed: vmUF must be <1!");
    
    const auto wmwN = w_from_matching(wpwN, vp, vm);
    if (wmwN <= wpwN) throw std::runtime_error("Detonation IC failed: wm>wp required!");

    const auto TmTN = get_TmTN(wmwN);
    if (TmTN <= TpTN) throw std::runtime_error("Detonation IC failed: Tm>Tp required!");

    const state_type y0 = {xi0, wmwN, TmTN};
    
    return {vmUF, y0};
}

// lambda profiles
double FluidProfile::lambda_b(double wowN) const {
    // la(xi) behind bubble wall (detonations)
    return (wowN - (1.0 + 3.0 * cmsq_ * alN_)) / (1.0 + cmsq_);
}

double FluidProfile::lambda_s(double wowN) const {
    // la(xi) in front of bubble wall (deflagrations)
    return (wowN - 1.0) / (1.0 + cpsq_);
}

/*****************************************************************************************/
/**************************************** Veff EoS ***************************************/
int FluidProfile::get_mode_veff(double vw, double cmsq) const {
    const auto cm = std::sqrt(cmsq);

    if (vw < cm) return 0; // deflagration

    // calculate Jouguet detonation velocity from matching eqs
    const auto vm = cm; // vJ_det=vp when |vm|=cm
    const auto TpTN = 1.0;

    const auto vp_min = vm;
    const auto vp_max = 1.0;
    const auto TmTN_min = TpTN;
    const auto TmTN_max = veff_params_->TTN_max();

    if (veff_params_->TTN_min() > TmTN_min)
        throw std::runtime_error("TmTN_min lies outside of bounds of EoS, check for bad EoS data!");

    const std::array<double, 2> bounds_min = {vp_min, TmTN_min};
    const std::array<double, 2> bounds_max = {vp_max, TmTN_max};

    auto matching_helper = [this, vm, TpTN](const std::array<double, 2>& vp_TmTN) {
        return matching_eqs_wall(vp_TmTN[0], TpTN, vm, vp_TmTN[1]);
    };

    // rough grid search to find initial guess for {vp, TmTN}
    const auto vp_TmTN_guess = grid_search_2d(matching_helper, bounds_min, bounds_max);
    // const auto vp_TmTN_guess = grid_search_2d(matching_helper, bounds_min, bounds_max, 50, 50, true, "get_mode_grid_search.csv");

    // solve matching eqs
    std::array<double, 2> sol;
    try { // solve using bounded newton's method
        sol = newton_solve_2d_bounded(matching_helper, vp_TmTN_guess, bounds_min, bounds_max, 1e-12, 100, 1e-12);
    } catch (std::exception& e) { // nelder-mead minimisation method fallback
        std::cerr << "WARNING: " << e.what()
                  << " Solving matching eqs by minimising residuals instead!\n";
        
        auto matching_helper2 = [this, vm, TpTN](const std::array<double, 2>& vp_TmTN) {
            const auto resi = matching_eqs_wall(vp_TmTN[0], TpTN, vm, vp_TmTN[1]);
            return std::array<double, 2>{resi[0] * resi[0], resi[1] * resi[1]};
        };

        sol = nelder_mead_minimise_2d(matching_helper2, vp_TmTN_guess[0], vp_TmTN_guess[1], 0.1, 0.01, 1e-12);

        const auto resi = matching_helper({sol[0], sol[1]});
        if (resi[0] > config::minimiser_tol || resi[1] > config::minimiser_tol) {
            throw std::runtime_error("Solving matching equations failed in get_mode_veff (residual too large)!");
        }
    }

    const auto vJ_det = sol[0]; // vp
    if (abs(vJ_det) <= cm) throw std::runtime_error("Invalid Jouguet detonation velocity calculated (|vJ_det| <= cm)!");
    
    const auto TmTN = sol[1];
    if (TmTN < TpTN) throw std::runtime_error("Invalid temperature for Jouguet detonation solution (TmTN < TpTN)!");

    if (vw < vJ_det) return 1; // hybrid
    return 2; // detonation
}

// change these to lambda functions in solve?
double FluidProfile::lambda_s_veff(double ToTN, const double eN, const double wN_inv) const {
    // la(xi) = (es(T(xi)) - eN) / wN
    const auto es_T = veff_params_->es_val(ToTN); // es(T/TN)
    return (es_T - eN) * wN_inv;
}

double FluidProfile::lambda_b_veff(double ToTN, const double eN, const double wN_inv) const {
    // la(xi) = (eb(T(xi)) - eN) / wN
    const auto eb_T = veff_params_->eb_val(ToTN); // eb(T/TN)
    return (eb_T - eN) * wN_inv;
}

std::array<double, 2> FluidProfile::matching_eqs_shock(double v1, double T1TN, double v2, double T2TN) const {    
    if (T1TN < veff_params_->TTN_min() || T1TN > veff_params_->TTN_max()) {
        throw std::runtime_error("T1/TN is called out of bounds for spline!");
    }
    if (T2TN < veff_params_->TTN_min() || T2TN > veff_params_->TTN_max()) {
        throw std::runtime_error("T2/TN is called out of bounds for spline!");
    }

    const auto p1 = veff_params_->ps_val(T1TN); // p_1, e_1
    const auto e1 = veff_params_->es_val(T1TN);
    const auto p2 = veff_params_->ps_val(T2TN); // p_2, e_2
    const auto e2 = veff_params_->es_val(T2TN);

    // scale residuals to O(1) (prevents issues with newton solver)
    // NOTE: p(T), e(T) ~ 1e+8 so residual must converge to at least 1e-10!!
    const auto scale = std::max({std::abs(p1), std::abs(p2), std::abs(e1), std::abs(e2)});
    const auto eq1 = (v2 * v1 * (e1 - e2) - (p1 - p2)) / scale;
    const auto eq2 = (v1 * (e1 + p2) - v2 * (e2 + p1)) / scale;

    // const auto eq1 = (v2 * v1 * (e1 - e2) - (p1 - p2));
    // const auto eq2 = (v1 * (e1 + p2) - v2 * (e2 + p1));

    return {eq1, eq2};
}

std::array<double, 2> FluidProfile::matching_eqs_wall(double vp, double TpTN, double vm, double TmTN) const {        
    if (TmTN < veff_params_->TTN_min() || TmTN > veff_params_->TTN_max()) {
        throw std::runtime_error("Tm/TN is called out of bounds for spline!");
    }
    if (TpTN < veff_params_->TTN_min() || TpTN > veff_params_->TTN_max()) {
        throw std::runtime_error("Tp/TN is called out of bounds for spline!");
    }

    const auto pp = veff_params_->ps_val(TpTN); // p_+, e_+
    const auto ep = veff_params_->es_val(TpTN);
    const auto pm = veff_params_->pb_val(TmTN); // p_-, e_-
    const auto em = veff_params_->eb_val(TmTN);

    // scale residuals to O(1) (prevents issues with newton solver)
    // NOTE: p(T), e(T) ~ 1e+8 so residual must converge to at least 1e-10!!
    const auto scale = std::max({std::abs(pp), std::abs(pm), std::abs(ep), std::abs(em)});
    const auto eq1 = (vp * vm * (em - ep) - (pm - pp)) / scale;
    const auto eq2 = (vm * (em + pp) - vp * (ep + pm)) / scale;

    return {eq1, eq2};
}

// dev
void FluidProfile::test_residual_veff(const deriv_func& dydv, const size_t n) const {
    std::cout << "Running test for T2/TN residual...\n";

    const auto TmTN_vals = linspace(veff_params_->TTN_min(), veff_params_->TTN_max(), n);
    // const auto TmTN_vals = linspace(1.03597, 1.03615, n);

    std::vector<double> resi_vals, TmTN_vals_pass;
    for (size_t i = 0; i < n; i++) {  
        double resi = std::numeric_limits<double>::quiet_NaN();
        try {
            resi = v1_residual_veff(dydv, TmTN_vals[i]);
            resi_vals.push_back(resi);
            TmTN_vals_pass.push_back(TmTN_vals[i]);
        } catch (std::exception& e) {
            // std::cout << "Failed for Tm/TN=" << TmTN_vals[i] << ": " << e.what() << "\n";
        }
    }

    if (resi_vals.size() == 0) {
        std::cout << "Warning: test_residual veff failed (resi_vals.size()=0)\n";
    }

    // index where residual is minimised
    const auto it = std::min_element(resi_vals.begin(), resi_vals.end());
    const auto idx = std::distance(resi_vals.begin(), it);
    std::cout << "test_residual_veff: TmTN=" << TmTN_vals_pass[idx] << ", min_resi=" << resi_vals[idx] << "\n";
    
    #ifdef ENABLE_MATPLOTLIB
    plt::figure_size(800, 800);
    plt::plot(TmTN_vals_pass, resi_vals);
    plt::xlabel("TmTN");
    plt::ylabel("residual");
    // plt::xlim(TmTN_vals.front(), TmTN_vals.back());
    plt::grid(true);
    plt::save("T2_resi.png");
    #endif

    std::cout << "Test complete. T2/TN residual saved to 'T2_resi.png'\n";
}

void FluidProfile::test_shock_veff(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol) const {
    std::cout << "Running test for shock residual...\n";

    std::vector<double> xi_vals, resi_vals;
    for (size_t i = 0; i < v_sol.size(); i++) {
        const auto xi_sh = y_sol[i][0];
        if (xi_sh <= vw_ || xi_sh > 1.0) continue;
        
        const auto v1UF = v_sol[i];
        const auto v1 = mu(xi_sh, abs(v1UF)); // v1=mu(xi_sh, v1UF)
        const auto T1TN = y_sol[i][2];

        const auto v2 = xi_sh; // v2=xi_sh

        // shock condition mu(xi_sh, v(xi_sh)) xi_sh = (p1-pN)/(e1-eN)
        const auto resi = matching_eqs_shock(v1, T1TN, v2, 1.0); // T2TN=1
        xi_vals.push_back(xi_sh);
        resi_vals.push_back(abs(resi[0]));
    }

    // index where residual is minimised
    const auto it = std::min_element(resi_vals.begin(), resi_vals.end());
    const auto idx = std::distance(resi_vals.begin(), it);
    std::cout << "test_shock_veff: xi_sh=" << xi_vals[idx] << ", min_resi=" << resi_vals[idx] << "\n";

    #ifdef ENABLE_MATPLOTLIB
    plt::figure_size(800, 800);
    plt::plot(xi_vals, resi_vals);
    plt::xlabel("xi_sh");
    plt::ylabel("residual");
    plt::grid(true);
    plt::save("shock_resi_veff.png");
    #endif

    std::cout << "Test complete. Shock residual saved to 'shock_resi_veff.png'\n";
}

double FluidProfile::find_TmTN_veff(const deriv_func& dydv, const bool fallback) const {
    auto safe_residual = [this, fallback, &dydv] (double TmTN) {
        try {
            const auto resi = v1_residual_veff(dydv, TmTN, fallback);
            if (!std::isfinite(resi)) return 1e+5;
            return resi;
        } catch (std::exception& e) {
            return 1e+5;
        }
    };

    // testing residual
    // test_residual_veff(dydv, 1000);

    // find bracket where residual is defined and minimum lies
    const auto TmTN_min = veff_params_->TTN_min();
    const auto TmTN_max = veff_params_->TTN_max();
    const auto bracket = find_bracket(safe_residual, TmTN_min, TmTN_max);
    if (!std::isfinite(bracket[0]) || !std::isfinite(bracket[1]))
        throw std::runtime_error("find_TmTN_veff failed (bracket not found)!");

    // find TmTN that minimise residual
    const auto TmTN = golden_section_minimize(safe_residual, bracket[0], bracket[1]);
    if (TmTN < 0.0) throw std::runtime_error("find_TmTN_veff failed (TmTN < 0)!");

    if (dev_log_) {
        std::cout << "find_TmTN_veff: TmTN=" << TmTN << ", resi=" << safe_residual(TmTN) << "\n";
    }

    return TmTN;
}

double FluidProfile::v1_residual_veff(const deriv_func& dydv, double TmTN, const bool fallback, const size_t n) const {
    const auto vm = (mode_ == 0) ? vw_ : std::sqrt(veff_params_->csq_b(TmTN)); // vm=vw (deflagration), vm=cm (hybrid)
    const auto wmwN = veff_params_->wb_val(TmTN) / veff_params_->wN();

    const auto [v_sol, y_sol] = deflagration_profile_veff(dydv, vm, wmwN, TmTN, fallback, false, n);

    // fluid behind shock
    const auto xi_sh = y_sol.back()[0];
    const auto v1UF = v_sol.back();
    const auto v1 = mu(xi_sh, abs(v1UF));
    // const auto w1wN = y_sol.back()[1];
    const auto T1TN = y_sol.back()[2];

    // residual for 2nd shock matching condition
    double resi;
    if (fallback) { // fallback routine if shock matching fails
        const auto TpTN = y_sol.front()[2];
        const auto cpsq = veff_params_->csq_s(TpTN);
        const auto T1TN4 = T1TN * T1TN * T1TN * T1TN;
        resi = v1 - xi_sh * (cpsq * T1TN4 + 1.0) / (T1TN4 + cpsq);
    } else {
        resi = matching_eqs_shock(v1, T1TN, xi_sh, 1.0)[1]; // T2/TN=1
    }

    return abs(resi);
}

std::tuple<double, double, double> FluidProfile::wall_matching_veff(const double vm, const double wmwN, const double TmTN) const {
    // matching at wall: vm,Tm -> vp,Tp
    const auto vp_min = 0.0;
    const auto vp_max = vm;
    // const auto vp_max = 1.0;
    const auto TpTN_min = veff_params_->TTN_min();
    const auto TpTN_max = veff_params_->TTN_max();

    const std::array<double, 2> bounds_min = {vp_min, TpTN_min};
    const std::array<double, 2> bounds_max = {vp_max, TpTN_max};

    auto matching_helper = [this, vm, TmTN](const std::array<double, 2>& vp_TpTN) {
        return matching_eqs_wall(vp_TpTN[0], vp_TpTN[1], vm, TmTN);
    };
    
    // rough grid search to find initial guess for {vp, TpTN}
    const auto vp_TpTN_guess = grid_search_2d(matching_helper, bounds_min, bounds_max);
    // const auto vp_TpTN_guess = grid_search_2d(matching_helper, bounds_min, bounds_max, 50, 50, true, "deflag_grid_search.csv");

    // solve matching eqs
    std::array<double, 2> sol;
    try { // solve using bounded newton's method
        sol = newton_solve_2d_bounded(matching_helper, vp_TpTN_guess, bounds_min, bounds_max);
    } catch (std::exception& e) { // nelder-mead minimisation method fallback  
        auto matching_helper2 = [this, vm, TmTN](const std::array<double, 2>& vp_TpTN) {
            const auto resi = matching_eqs_wall(vp_TpTN[0], vp_TpTN[1], vm, TmTN);
            return std::array<double, 2>{resi[0] * resi[0], resi[1] * resi[1]};
        };

        sol = nelder_mead_minimise_2d(matching_helper2, vp_TpTN_guess[0], vp_TpTN_guess[1]);

        const auto resi = matching_helper({sol[0], sol[1]});
        if (resi[0] > config::minimiser_tol || resi[1] > config::minimiser_tol) {
            throw std::runtime_error("Solving matching equations failed in deflagration_profile_veff (residual too large)!");
        }
    }

    // fluid in front of wall
    const auto vp = sol[0];
    const auto wpwN = w_from_matching(wmwN, vm, vp);
    const auto TpTN = sol[1];

    return {vp, wpwN, TpTN};
}

std::pair<std::vector<double>, std::vector<state_type>> FluidProfile::deflagration_profile_veff(const deriv_func& dydv, double vm, double wmwN, double TmTN, const bool fallback, const bool final_prof, const size_t n) const {
    const auto [vp, wpwN, TpTN] = wall_matching_veff(vm, wmwN, TmTN);
    const auto vpUF = mu(vw_, abs(vp));

    const state_type yp = {vw_, wpwN, TpTN}; // {xi_w, wpwN, TpTN}

    // solve EoM - integrate from v(xi)=vpUF -> v(xi)=0
    auto [v_sol_tmp, y_sol_tmp] = rk4_solver(dydv, vpUF, 1e-10, yp, n);

    // truncate solution to where dxi/dv=0 (gives better behaved shock residual)
    std::vector<double> resi_vals(n);
    for (size_t i = 0; i < v_sol_tmp.size(); i++) {
        resi_vals[i] = abs(dxidv(y_sol_tmp[i][0], v_sol_tmp[i], veff_params_->csq_s(y_sol_tmp[i][2])));
    }
    const auto it = std::min_element(resi_vals.begin(), resi_vals.end());
    const auto dxidv_root_idx = std::distance(resi_vals.begin(), it); // idx where dxi/dv=0

    // only erases part of solution in front of dxidv=0
    // if dxidv /= 0 across whole solution, ignores this
    if (resi_vals[dxidv_root_idx] < 1e-2) {
        v_sol_tmp.erase(v_sol_tmp.begin() + dxidv_root_idx + 1, v_sol_tmp.end());
        y_sol_tmp.erase(y_sol_tmp.begin() + dxidv_root_idx + 1, y_sol_tmp.end());
    }

    // find shock & truncate solution
    const auto find_sh = find_shock_idx_veff(v_sol_tmp, y_sol_tmp, fallback); // uses first matching cond. at shock
    const auto sh_idx = std::min(find_sh, v_sol_tmp.size() - 1);

    // re-integrate for final profile
    // avoids insufficient no. points for integrating profile when shock is close to vpUF
    if (final_prof) {
        const auto [v_sol, y_sol] = rk4_solver(dydv, vpUF, v_sol_tmp[sh_idx], yp, n); // integrates from vpUF->v1UF
        return {v_sol, y_sol};
    }

    const std::vector<double> v_sol(v_sol_tmp.begin(), v_sol_tmp.begin() + sh_idx + 1);
    const std::vector<state_type> y_sol(y_sol_tmp.begin(), y_sol_tmp.begin() + sh_idx + 1);

    return {v_sol, y_sol};
}

bool FluidProfile::check_shock_convergence(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol) {
    const auto resi = matching_eqs_shock(mu(y_sol.back()[0], v_sol.back()), y_sol.back()[2], y_sol.back()[0], 1.0); 

    // check shock matching condition for final profile
    if (abs(resi[0]) > config::minimiser_tol || abs(resi[1]) > config::minimiser_tol) {
        shock_flag_ = false;
        return shock_flag_;
    }

    return true;
}

bool FluidProfile::check_shock_convergence_fallback(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol) const {
    const auto xi_sh = y_sol.back()[0];
    const auto v1 = mu(y_sol.back()[0], v_sol.back());
    const auto T1TN = y_sol.back()[2];
    const auto T1TN4 = T1TN * T1TN * T1TN * T1TN;

    const auto TpTN = y_sol.front()[2];
    const auto cpsq = veff_params_->csq_s(TpTN);

    const std::array<double, 2> resi = {v1 * xi_sh - cpsq, v1 - xi_sh * (cpsq * T1TN4 + 1.0) / (T1TN4 + cpsq)};

    // check shock matching condition for final profile
    if (abs(resi[0]) > config::sh_fallback_tol || abs(resi[1]) > config::sh_fallback_tol) return false;
    return true;
}

size_t FluidProfile::find_shock_idx_veff(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol, const bool fallback) const {
    const auto TpTN = y_sol[0][2];
    const auto cpsq = veff_params_->csq_s(TpTN);

    std::vector<double> resi_vals(v_sol.size());
    int pass_count = 0;

    for (size_t i = 0; i < v_sol.size(); i++) {
        const auto xi_sh = y_sol[i][0];

        if (xi_sh <= vw_ || xi_sh > 1.0) {
            resi_vals[i] = 1.0;
            continue;
        }
        
        const auto v1UF = v_sol[i];
        const auto v1 = mu(xi_sh, abs(v1UF)); // v1=mu(xi_sh, v1UF)
        const auto T1TN = y_sol[i][2];
        const auto v2 = xi_sh; // v2=xi_sh

        // shock condition mu(xi_sh, v(xi_sh)) xi_sh = (p1-pN)/(e1-eN)
        if (fallback) {
            resi_vals[i] = abs(v1 * v2 - cpsq);
        } else {
            const auto resi = matching_eqs_shock(v1, T1TN, v2, 1.0); // T2/TN=1
            resi_vals[i] = abs(resi[0]);
        }
        
        pass_count++;
    }

    if (pass_count == 0)
        throw std::runtime_error("find_shock_idx_veff failed (no shock found in fluid profile)!");

    // index where residual is minimised
    const auto it = std::min_element(resi_vals.begin(), resi_vals.end());
    auto idx = std::distance(resi_vals.begin(), it);

    return idx;
}

std::pair<double, state_type> FluidProfile::get_IC_detonation_veff(const double vm_min, const double vm_max) const {
    // uses matching conditions to get vm, Tm/TN from vp, Tp/TN
    const auto xi0 = vw_;
    const auto vp = vw_;
    const auto TpTN = 1.0;
    const auto wpwN = 1.0;

    // newton solver seems to fail sometimes when tighter bounds on vm and Tm are used
    // i.e. 0 < vm < vp and TpTN < TmTN < TTN_max
    // leave as is and have consistency check to test vm < vp and Tm > Tp after solver
    const auto TmTN_min = veff_params_->TTN_min();
    const auto TmTN_max = veff_params_->TTN_max();

    const std::array<double, 2> bounds_min = {vm_min, TmTN_min};
    const std::array<double, 2> bounds_max = {vm_max, TmTN_max};

    auto matching_helper = [this, vp, TpTN] (std::array<double, 2> vm_TmTN) {
        return matching_eqs_wall(vp, TpTN, vm_TmTN[0], vm_TmTN[1]); // ym_guess = {vm_guess, TmTN_guess}
    };

    // rough grid search to find initial guess for {vm, TmTN}
    auto vm_TmTN_guess = grid_search_2d(matching_helper, bounds_min, bounds_max);
    // auto vm_TmTN_guess = grid_search_2d(matching_helper, bounds_min, bounds_max, 100, 100, true, "det_grid_search.csv");

    // prevents solution from landing in wrong minima if grid search fails
    if (dxidv(vw_, mu(vw_, vm_TmTN_guess[0]), veff_params_->csq_b(vm_TmTN_guess[1])) < 0.0) {
        vm_TmTN_guess = {vp, TpTN};
    }

    // solve matching eqs
    std::array<double, 2> sol;
    try { // solve using bounded newton's method
        sol = newton_solve_2d_bounded(matching_helper, vm_TmTN_guess, bounds_min, bounds_max, 1e-12, 100, 1e-12);
    } catch (std::exception& e) { // nelder-mead minimisation method fallback
        std::cerr << "WARNING: " << e.what()
                  << " Solving matching eqs by minimising residuals instead!\n";
        
        auto matching_helper2 = [this, vp, TpTN] (std::array<double, 2> vm_TmTN) {
            const auto resi = matching_eqs_wall(vp, TpTN, vm_TmTN[0], vm_TmTN[1]);
            return std::array<double, 2>{resi[0] * resi[0], resi[1] * resi[1]};
        };

        sol = nelder_mead_minimise_2d(matching_helper2, vm_TmTN_guess[0], vm_TmTN_guess[1], 0.1, 0.01, 1e-12);

        const auto resi = matching_helper({sol[0], sol[1]});
        if (resi[0] > config::minimiser_tol || resi[1] > config::minimiser_tol) {
            throw std::runtime_error("Solving matching equations failed in get_IC_detonation_veff (residual too large)!");
        }
    }

    // fluid behind wall
    const auto vm = sol[0];
    if (abs(vm) >= abs(vp)) throw std::runtime_error("|vm| must be < |vp|=vw for detonation");

    const auto vmUF = mu(vw_, abs(vm));
    if (vmUF < 0.0) throw std::runtime_error("vmUF must be > 0 for detonation");

    const auto TmTN = sol[1];
    if (TmTN <= TpTN) throw std::runtime_error("Tm/TN must be > Tp/TN for detonation");

    const auto wmwN = w_from_matching(wpwN, vp, vm);
    if (wmwN <= wpwN) throw std::runtime_error("wm/wN must be > wp/wN for detonation");

    const auto resi = matching_eqs_wall(vp, TpTN, vm, TmTN);
    std::cout << "get_IC_detonation_veff: vm=" << vm << ", TmTN=" << TmTN << ", resi[0]=" << resi[0] << ", resi[1]=" << resi[1] << "\n";

    const state_type y0 = {xi0, wmwN, TmTN};
    return {vmUF, y0};
}

/*****************************************************************************************/
/************************************* Both Bag/Veff *************************************/
// from matching condition wp * vp * gamma_p^2 = wm * vm * gamma_m^2
double FluidProfile::w_from_matching(double wp, double vp, double vm) const {
    return wp * abs(vp) * (1.0 - vm * vm) / (abs(vm) * (1.0 - vp * vp));
}
/*****************************************************************************************/

std::vector<prof_type> FluidProfile::solve_profile(int n) {
    // check valid hydrodynamic mode
    if (!(mode_ == 0 || mode_ == 1 || mode_ == 2)) {
            throw std::runtime_error("Hydrodynamic mode must be: 0 (deflagration), 1 (hybrid) or 2 (detonation)");
    }

    std::cout << "Solving fluid profile for hydrodynamic mode=";
    if (mode_ == 0) {
        std::cout << "deflagration...";
    } else if (mode_ == 1) {
        std::cout << "hybrid...";
    } else {
        std::cout << "detonation...";
    }
    std::cout << "\n";

    auto dydv = [this] (double vUF, const state_type& y) -> state_type {
        return dydv_vec(vUF, y, vw_, cmsq_, cpsq_);
    };

    double xi0, xif;
    std::vector<state_type> y_sol_tmp;
    prof_type xi_sol_tmp, v_sol_tmp, w_sol_tmp, T_sol_tmp, la_sol_tmp;

    const auto w_start_val = 1.0; // w+/wN (det), w2/wN (deflag/hybrid)
    const auto T_start_val = 1.0; // T+/TN (det), T2/TN (deflag/hybrid);
    const auto la_start_val = 0.0;
    
    double w_end_val, T_end_val, la_end_val;

    if (mode_ < 2) { // deflagration & hybrid
        // hybrid and deflagration ICs the same for xi_w < xi < xi_sh
        // v(xi_sh) = v1UF, w(xi_sh) = w1wN, T(xi_sh) = T1TN

        const auto vpUF = find_vpUF(dydv);
        const auto vp = mu(vw_, vpUF);

        const auto prof_tmp = deflagration_profile(dydv, vpUF, true, n);
        v_sol_tmp = std::get<0>(prof_tmp);
        y_sol_tmp = std::get<1>(prof_tmp);
        const bool shock_flag = std::get<2>(prof_tmp);


        std::reverse(v_sol_tmp.begin(), v_sol_tmp.end());
        std::reverse(y_sol_tmp.begin(), y_sol_tmp.end());

        for (size_t i = 0; i < v_sol_tmp.size(); i++) {
            xi_sol_tmp.push_back(y_sol_tmp[i][0]);
            w_sol_tmp.push_back(y_sol_tmp[i][1]);
            T_sol_tmp.push_back(y_sol_tmp[i][2]);

            la_sol_tmp.push_back(lambda_s(w_sol_tmp[i]));
        }

        // cp <= xi_sh < 1
        // solver will find xi_sh < cp if shock is exactly xi_sh=cp (numerical error)
        // ensures difference isn't too large!
        const auto cp = std::sqrt(cpsq_);

        const auto xi_sh = (shock_flag) ? cp : std::max(xi_sol_tmp.front(), cp);
        if (xi_sh > 1.0) {
            throw std::runtime_error("Shock must be less than speed of light (cp <= xi_sh < 1)");
        }

        const auto wpwN = w_sol_tmp.back();
        const auto TpTN = T_sol_tmp.back();

        const auto v1UF = v_sol_tmp.front();
        const auto w1wN = w_sol_tmp.front();
        const auto T1TN = T_sol_tmp.front();

        xi0 = xi_sh;
        xif = vw_;

        // check alp okay
        auto alp_wall = [this] (double vp, double vm) {
            const auto alp = get_alp_wall(vp, vm);

            if (alp >= alN_) throw std::runtime_error("alpha_+ must be < alpha_N");
            if (alp < alp_min_) throw std::runtime_error("alpha_+ too small for shock");
            if (alp > alp_max_) throw std::runtime_error("alpha_+ too large for shock");

            return alp;
        };

        if (mode_ == 0) {
            // fix end-value for enthalpy
            const auto vm = -vw_;
            const auto alp = alp_wall(vp, abs(vm));

            const auto wpwN = w_sol_tmp.back(); // w(xi_w + dlt) = w+/wN

            w_end_val = w_from_matching(wpwN, vp, vm); // wm/wN, from matching condition at wall
            T_end_val = get_TmTN(w_end_val);
            la_end_val = lambda_b(w_end_val);

            if (dev_log_) {
                std::clog << "Deflagration profile:\n"
                          << "  vm = " << vm << ", vmUF = " << mu(vw_, abs(vm)) << "\n"
                          << "  wmwN = " << w_end_val << ", TmTN = " << T_end_val << "\n"
                          << "  vp = " << vp << ", vpUF = " << vpUF << "\n"
                          << "  wpwN = " << wpwN << ", TpTN = " << T_sol_tmp.back() << "\n"
                          << "  alp = " << alp << "\n"
                          << "  v1 = " << mu(xi_sh, abs(v1UF)) << ", v1UF = " << v1UF << "\n"
                          << "  w1wN = " << w1wN << ", T1TN = " << T1TN << "\n"
                          << "  xi_sh = " << xi_sh << "\n";
            }
                      

        } else { // hybrid
            const auto vm = -std::sqrt(cmsq_);
            const auto alp = alp_wall(vp, abs(vm));

            // consistency check - I'm not sure why this is violated for some hybrids...
            // if (vw_ >= vJ_det(alp)) {
            //     std::cout << "Warning: Hybrid condition violated (vw >= vJ_det). Hydrodynamic mode may be incorrect!\n";
            // }
            
            // initial conditions for rarefaction wave
            const auto xi0_rf = vw_;
            
            const auto vmUF = mu(vw_, abs(vm));

            const auto vp = mu(vw_, abs(vpUF));
            const auto wmwN = w_from_matching(wpwN, vp, vm);

            const auto TmTN = get_TmTN(wmwN);

            const state_type y0_rf = {xi0_rf, wmwN, TmTN};
            auto [v_sol_rf_tmp, y_sol_rf_tmp] = rk4_solver(dydv, vmUF, 1e-10, y0_rf, n);

            // remove any numerical errors in final point (xi<0)
            if (!y_sol_rf_tmp.empty() && y_sol_rf_tmp.back()[0] < 0.0) {
                v_sol_rf_tmp.pop_back();
                y_sol_rf_tmp.pop_back();
            }

            // combine rarefaction wave with shockwave part of solution
            for (size_t i = 0; i < v_sol_rf_tmp.size(); i++) {
                v_sol_tmp.push_back(v_sol_rf_tmp[i]);
                xi_sol_tmp.push_back(y_sol_rf_tmp[i][0]);
                w_sol_tmp.push_back(y_sol_rf_tmp[i][1]);
                T_sol_tmp.push_back(y_sol_rf_tmp[i][2]);

                la_sol_tmp.push_back(lambda_b(y_sol_rf_tmp[i][1]));
            }

            // xif = xif_rf; // update xif value to behind rarefaction wave
            w_end_val = w_sol_tmp.back();
            T_end_val = T_sol_tmp.back();
            la_end_val = la_sol_tmp.back();

            if (dev_log_) {
                std::clog << "Hybrid profile:\n"
                          << "  vm = " << vm << ", vmUF=" << mu(vw_, abs(vm)) << "\n"
                          << "  wmwN = " << w_end_val << ", TmTN = " << TmTN << "\n"
                          << "  vp = " << vp << ", vpUF = " << vpUF << "\n"
                          << "  wpwN = " << wpwN << ", TpTN = " << TpTN << "\n"
                          << "  alp = " << alp << "\n"
                          << "  v1 = " << mu(xi_sh, abs(v1UF)) << ", v1UF = " << v1UF << "\n"
                          << "  w1wN = " << w1wN << ", T1TN = " << T1TN << "\n"
                          << "  xi_sh = " << xi_sh << "\n";
            }
        }

    } else { // detonation
        // cm < xi < xi_w
        xi0 = vw_;
        // xif = std::sqrt(cmsq_);

        // solver
        const auto ics = get_IC_detonation();
        const auto vmUF = ics.first;
        const auto y0 = ics.second; // {xi0, w0, T0}

        const auto sol = rk4_solver(dydv, vmUF, 1e-10, y0, n);
        v_sol_tmp = sol.first;
        y_sol_tmp = sol.second;

        // remove any numerical errors in final point (xi<0)
        if (!y_sol_tmp.empty() && y_sol_tmp.back()[0] < 0.0) {
            v_sol_tmp.pop_back();
            y_sol_tmp.pop_back();
        }

        // fill v, w vectors
        // Optimisation Note: pre-allocate memory for these vectors to speed up?
        for (size_t i = 0; i < v_sol_tmp.size(); i++) {
            xi_sol_tmp.push_back(y_sol_tmp[i][0]);
            w_sol_tmp.push_back(y_sol_tmp[i][1]);
            T_sol_tmp.push_back(y_sol_tmp[i][2]);

            la_sol_tmp.push_back(lambda_b(w_sol_tmp[i]));
        }

        w_end_val = w_sol_tmp.back();
        T_end_val = T_sol_tmp.back();
        la_end_val = la_sol_tmp.back();

        if (dev_log_) {
            std::clog << "Detonation profile:\n"
                      << "  vm = " << mu(vw_, abs(vmUF)) << ", vmUF=" << vmUF << "\n"
                      << "  wmwN = " << y0[1] << ", TmTN = " << y0[2] << "\n"
                      << "  w_end = " << w_end_val << ", T_end = " << T_end_val << "\n";
        }
    }

    // update final point manually where dxidv is singular
    if (mode_ != 0) {
        w_end_val = w_sol_tmp.back() - dwdv(std::sqrt(cmsq_), 0.0, w_end_val, cmsq_) * v_sol_tmp.back();
        T_end_val = T_sol_tmp.back() - dTdv(std::sqrt(cmsq_), 0.0, T_end_val) * v_sol_tmp.back();
        la_end_val = lambda_b(w_end_val);

        const size_t idx = v_sol_tmp.size() - 1; // index of last point
        v_sol_tmp[idx] = 0.0;
        xi_sol_tmp[idx] = std::sqrt(cmsq_);
        w_sol_tmp[idx] = w_end_val;
        T_sol_tmp[idx] = T_end_val;
        la_sol_tmp[idx] = la_end_val;

        xif = xi_sol_tmp.back();
    }

    // store start/endpoints of profile for integration
    xi_min_integrate_ = xif;
    xi_max_integrate_ = xi0;
    
    // define start & end points where profile=const (outside integration)
    const auto xi_start = linspace(0.99, xi0, n); // backwards integration
    const auto xi_end = linspace(xif, 0.01, n);

    const prof_type v_start(n, 0.0);
    const prof_type v_end = v_start;

    const prof_type w_start(n, w_start_val);
    const prof_type w_end(n, w_end_val);

    const prof_type T_start(n, T_start_val);
    const prof_type T_end(n, T_end_val);

    const prof_type la_start(n, la_start_val);
    const prof_type la_end(n, la_end_val);

    prof_type xi_sol, v_sol, w_sol, T_sol, la_sol;

    // concatenate xi vals
    xi_sol.insert(xi_sol.end(), xi_start.begin(), xi_start.end());
    xi_sol.insert(xi_sol.end(), xi_sol_tmp.begin(), xi_sol_tmp.end());
    xi_sol.insert(xi_sol.end(), xi_end.begin(), xi_end.end());

    // concatenate v(xi) vals
    v_sol.insert(v_sol.end(), v_start.begin(), v_start.end());
    v_sol.insert(v_sol.end(), v_sol_tmp.begin(), v_sol_tmp.end());
    v_sol.insert(v_sol.end(), v_end.begin(), v_end.end());

    // concatenate w(xi) vals
    w_sol.insert(w_sol.end(), w_start.begin(), w_start.end());
    w_sol.insert(w_sol.end(), w_sol_tmp.begin(), w_sol_tmp.end());
    w_sol.insert(w_sol.end(), w_end.begin(), w_end.end());

    // concatenate T(xi) vals
    T_sol.insert(T_sol.end(), T_start.begin(), T_start.end());
    T_sol.insert(T_sol.end(), T_sol_tmp.begin(), T_sol_tmp.end());
    T_sol.insert(T_sol.end(), T_end.begin(), T_end.end());

    // concatenate la(xi) vals
    la_sol.insert(la_sol.end(), la_start.begin(), la_start.end());
    la_sol.insert(la_sol.end(), la_sol_tmp.begin(), la_sol_tmp.end());
    la_sol.insert(la_sol.end(), la_end.begin(), la_end.end()); 

    // reformat from backwards integration
    std::reverse(xi_sol.begin(), xi_sol.end());
    std::reverse(v_sol.begin(), v_sol.end());
    std::reverse(w_sol.begin(), w_sol.end());
    std::reverse(T_sol.begin(), T_sol.end());
    std::reverse(la_sol.begin(), la_sol.end());
    
    return {xi_sol, v_sol, w_sol, T_sol, la_sol};
}

std::vector<prof_type> FluidProfile::solve_profile_veff(int n) {
    std::cout << "Solving fluid profile for hydrodynamic mode=";
    if (mode_ == 0) {
        std::cout << "deflagration...";
    } else if (mode_ == 1) {
        std::cout << "hybrid...";
    } else if (mode_ == 2){
        std::cout << "detonation...";
    } else {
        throw std::runtime_error("Hydrodynamic mode must be: 0 (deflagration), 1 (hybrid) or 2 (detonation)");
    }
    std::cout << "\n";

    // wrapper for hydrodynamic equation of motion
    auto dydv = [this] (double vUF, const state_type& y) -> state_type {
        return dydv_vec(vUF, y, vw_, veff_params_->csq_b(y[2]), veff_params_->csq_s(y[2]));
    };

    const auto eN = veff_params_->eN();
    const auto wN_inv = 1.0 / veff_params_->wN();

    double xi0, xif;
    std::vector<state_type> y_sol_tmp;
    prof_type xi_sol_tmp, v_sol_tmp, w_sol_tmp, T_sol_tmp, la_sol_tmp;

    const auto w_start_val = 1.0; // w+/wN (det), w2/wN (deflag/hybrid)
    const auto T_start_val = 1.0; // T+/TN (det), T2/TN (deflag/hybrid);
    const auto la_start_val = 0.0;
    
    double w_end_val, T_end_val, la_end_val;

    if (mode_ < 2) { // deflagration & hybrid
        // hybrid and deflagration ICs the same for xi_w < xi < xi_sh
        // v(xi_sh) = v1UF, w(xi_sh) = w1wN, T(xi_sh) = T1TN

        auto TmTN = find_TmTN_veff(dydv);
        auto vm = (mode_ == 0) ? vw_ : std::sqrt(veff_params_->csq_b(TmTN)); // vm=vw (deflagration), vm=cm (hybrid)
        auto wmwN = veff_params_->wb_val(TmTN) / veff_params_->wN();
        auto prof_tmp = deflagration_profile_veff(dydv, vm, wmwN, TmTN, false, true, n);
        
        // check convergence of profile
        const bool pass = check_shock_convergence(prof_tmp.first, prof_tmp.second);
        if (!pass) { // fallback method
            std::cerr << "Warning: Shock residual above tolerance. Using mu-nu fallback method for shock finding!\n";
            TmTN = find_TmTN_veff(dydv, true);
            vm = (mode_ == 0) ? vw_ : std::sqrt(veff_params_->csq_b(TmTN));
            wmwN = veff_params_->wb_val(TmTN) / veff_params_->wN();
            prof_tmp = deflagration_profile_veff(dydv, vm, wmwN, TmTN, true, true, n);

            if (!check_shock_convergence_fallback(prof_tmp.first, prof_tmp.second)) {
                throw std::runtime_error("Shock convergence failed!");
            }
        }

        v_sol_tmp = prof_tmp.first;
        y_sol_tmp = prof_tmp.second;

        std::reverse(v_sol_tmp.begin(), v_sol_tmp.end());
        std::reverse(y_sol_tmp.begin(), y_sol_tmp.end());

        // fill v, w vectors
        for (size_t i = 0; i < v_sol_tmp.size(); i++) {
            xi_sol_tmp.push_back(y_sol_tmp[i][0]);
            w_sol_tmp.push_back(y_sol_tmp[i][1]);
            
            const auto T_sol = y_sol_tmp[i][2];
            T_sol_tmp.push_back(T_sol);
            la_sol_tmp.push_back(lambda_s_veff(T_sol, eN, wN_inv)); // lambda in symmetric phase
        }

        xi0 = xi_sol_tmp.front(); // xi_sh
        if (xi0 <= vw_) throw std::runtime_error("solve_profile failed (xi_sh must be > vw)");

        const auto vpUF = v_sol_tmp.back();
        if (vpUF <= 0.0) throw std::runtime_error("solve_profile failed (vpUF < 0)!");

        const auto wpwN = w_sol_tmp.back();
        if (wpwN <= 0.0) throw std::runtime_error("solve_profile failed (wpwN < 0)!");
        if (wpwN <= 1.0) throw std::runtime_error("solve_profile failed (wp/wN must be > 1)!");
        if (wpwN <= w_sol_tmp.front()) throw std::runtime_error("solve_profile failed (wp/wN must be > w1/wN)!");

        const auto TpTN = T_sol_tmp.back();
        if (TpTN <= 0.0) throw std::runtime_error("solve_profile failed (TpTN < 0)!");
        if (TpTN <= 1.0) throw std::runtime_error("solve_profile failed (Tp/TN must be > 1)!");
        if (TpTN <= T_sol_tmp.front()) throw std::runtime_error("solve_profile failed (Tp/TN must be > T1/TN)!");

        // update sound speed values
        cpsq_ = veff_params_->csq_s(TpTN);
        cmsq_ = veff_params_->csq_b(TmTN);

        if (mode_ == 0) { // deflagration
            xif = vw_;

            // fix end values 
            w_end_val = wmwN;
            T_end_val = TmTN;
            la_end_val = lambda_b_veff(T_end_val, eN, wN_inv); // lambda just behind wall (broken phase)

            if (dev_log_) {
                std::clog << "Deflagration profile:\n"
                          << "  vm = " << vm << ", vmUF=" << mu(vw_, abs(vm)) << "\n"
                          << "  wmwN = " << wmwN << ", TmTN = " << TmTN << "\n"
                          << "  vp = " << mu(vw_, abs(vpUF)) << ", vpUF = " << vpUF << "\n"
                          << "  wpwN = " << wpwN << ", TpTN = " << T_sol_tmp.back() << "\n"
                          << "  v1 = " << mu(xi0, abs(v_sol_tmp.front())) << ", v1UF = " << v_sol_tmp.front() << "\n"
                          << "  w1wN = " << w_sol_tmp.front() << ", T1TN = " << T_sol_tmp.front() << "\n"
                          << "  xi_sh = " << xi0 << "\n"
                          << "  cp = " << std::sqrt(veff_params_->csq_s(T_sol_tmp.back())) << ", cm = " << std::sqrt(veff_params_->csq_b(TmTN)) << "\n";
            }
        } else { // hybrid
            // initial conditions for rarefaction wave
            const auto xi0_rf = vw_;

            const auto vmUF = mu(vw_, abs(vm));

            // TO DO: fix this issue
            if (vmUF < 0.0) throw std::runtime_error("solve_profile_veff failed (vmUF < 0)!"); // fails if detonation part of solution is too small

            const state_type y0 = {xi0_rf, wmwN, TmTN}; // {xi0, wmwN, TmTN}
            auto [v_sol_rf_tmp, y_sol_rf_tmp] = rk4_solver(dydv, vmUF, 1e-10, y0, n);

            // remove any numerical errors in final point (xi<0)
            if (!y_sol_rf_tmp.empty() && y_sol_rf_tmp.back()[0] < 0.0) {
                v_sol_rf_tmp.pop_back();
                y_sol_rf_tmp.pop_back();
            }

            // combine rarefaction wave with shockwave part of solution
            for (size_t i = 0; i < v_sol_rf_tmp.size(); i++) {
                v_sol_tmp.push_back(v_sol_rf_tmp[i]);
                xi_sol_tmp.push_back(y_sol_rf_tmp[i][0]);
                w_sol_tmp.push_back(y_sol_rf_tmp[i][1]);
            
                const auto T_sol_rf = y_sol_rf_tmp[i][2];
                T_sol_tmp.push_back(T_sol_rf);
                la_sol_tmp.push_back(lambda_b_veff(T_sol_rf, eN, wN_inv));
            }

            w_end_val = w_sol_tmp.back();
            T_end_val = T_sol_tmp.back();
            la_end_val = la_sol_tmp.back();

            if (dev_log_) {
                std::clog << "Hybrid profile (veff):\n"
                          << "  vm = " << vm << ", vmUF=" << vmUF << "\n"
                          << "  wmwN = " << wmwN << ", TmTN = " << TmTN << "\n"
                          << "  vp = " << mu(vw_, abs(vpUF)) << ", vpUF = " << vpUF << "\n"
                          << "  wpwN = " << wpwN << ", TpTN = " << TpTN << "\n"
                          << "  v1 = " << mu(xi0, abs(v_sol_tmp.front())) << ", v1UF = " << v_sol_tmp.front() << "\n"
                          << "  w1wN = " << w_sol_tmp.front() << ", T1TN = " << T_sol_tmp.front() << "\n"
                          << "  xi_sh = " << xi0 << "\n"
                          << "  cp = " << std::sqrt(veff_params_->csq_s(TpTN)) << ", cm = " << std::sqrt(veff_params_->csq_b(TmTN)) << "\n";
            }
        }

    } else { // detonation
        // cm < xi < xi_w
        // v0 = v(xi_w) = vm(UF)
        // w0 = w(xi_w) = wm/wN
        // T0 = T(xi_w) = Tm/TN

        auto ics = get_IC_detonation_veff();
        auto vmUF = ics.first;
        auto y0 = ics.second; // {xi_w, wmwN, TmTN}

        // failsafe for if vmUF < cm
        if (vmUF < std::sqrt(veff_params_->csq_b(y0[2]))) {
            ics = get_IC_detonation_veff(std::sqrt(veff_params_->csq_b(y0[2])), 1.0);
            vmUF = ics.first;
            y0 = ics.second;
        }

        // solver
        const auto sol = rk4_solver(dydv, vmUF, 1e-10, y0, n);
        v_sol_tmp = sol.first;
        y_sol_tmp = sol.second;

        // remove any numerical errors in final point (xi<0)
        if (!y_sol_tmp.empty() && y_sol_tmp.back()[0] < 0.0) {
            v_sol_tmp.pop_back();
            y_sol_tmp.pop_back();
        }

        // fill v, w vectors
        for (size_t i = 0; i < v_sol_tmp.size(); i++) {
            xi_sol_tmp.push_back(y_sol_tmp[i][0]);
            w_sol_tmp.push_back(y_sol_tmp[i][1]);

            const auto T_sol = y_sol_tmp[i][2];
            T_sol_tmp.push_back(T_sol);
            la_sol_tmp.push_back(lambda_b_veff(T_sol, eN, wN_inv));
        }

        xi0 = vw_;

        w_end_val = w_sol_tmp.back();
        T_end_val = T_sol_tmp.back();
        la_end_val = la_sol_tmp.back();

        // update sound speed values (cp doesn't change since Tp=TN for detonations)
        cmsq_ = veff_params_->csq_b(y0[2]);

        if (dev_log_) {
            std::cout << "Detonation profile:\n"
                      << "  vm = " << mu(vw_, abs(vmUF)) << ", vmUF=" << vmUF << "\n"
                      << "  wmwN = " << y0[1] << ", TmTN = " << y0[2] << "\n";
        }
    }

    // update final point manually where dxidv is singular
    if (mode_ != 0) {
        w_end_val = w_sol_tmp.back() - dwdv(std::sqrt(veff_params_->csq_b(T_end_val)), 0.0, w_end_val, veff_params_->csq_b(T_end_val)) * v_sol_tmp.back();
        T_end_val = T_sol_tmp.back() - dTdv(std::sqrt(veff_params_->csq_b(T_end_val)), 0.0, T_end_val) * v_sol_tmp.back();
        la_end_val = lambda_b_veff(T_end_val, eN, wN_inv);

        const size_t idx = v_sol_tmp.size() - 1; // index of last point
        v_sol_tmp[idx] = 0.0;
        xi_sol_tmp[idx] = std::sqrt(veff_params_->csq_b(T_end_val));
        w_sol_tmp[idx] = w_end_val;
        T_sol_tmp[idx] = T_end_val;
        la_sol_tmp[idx] = la_end_val;

        xif = xi_sol_tmp.back();
    }

    // store start/endpoints of profile for integration
    xi_min_integrate_ = xif;
    xi_max_integrate_ = xi0;

    // define start & end points where profile=const (outside integration)
    const auto xi_start = linspace(0.99, xi0, n); // backwards integration
    const auto xi_end = linspace(xif, 0.01, n);

    const prof_type v_start(n, 0.0);
    const prof_type v_end = v_start;

    const prof_type w_start(n, w_start_val);
    const prof_type w_end(n, w_end_val);

    const prof_type T_start(n, T_start_val);
    const prof_type T_end(n, T_end_val);

    const prof_type la_start(n, la_start_val);
    const prof_type la_end(n, la_end_val);

    prof_type xi_sol, v_sol, w_sol, T_sol, la_sol;

    // concatenate xi vals
    xi_sol.insert(xi_sol.end(), xi_start.begin(), xi_start.end());
    xi_sol.insert(xi_sol.end(), xi_sol_tmp.begin(), xi_sol_tmp.end());
    xi_sol.insert(xi_sol.end(), xi_end.begin(), xi_end.end());

    // concatenate v(xi) vals
    v_sol.insert(v_sol.end(), v_start.begin(), v_start.end());
    v_sol.insert(v_sol.end(), v_sol_tmp.begin(), v_sol_tmp.end());
    v_sol.insert(v_sol.end(), v_end.begin(), v_end.end());

    // concatenate w(xi) vals
    w_sol.insert(w_sol.end(), w_start.begin(), w_start.end());
    w_sol.insert(w_sol.end(), w_sol_tmp.begin(), w_sol_tmp.end());
    w_sol.insert(w_sol.end(), w_end.begin(), w_end.end());

    // concatenate T(xi) vals
    T_sol.insert(T_sol.end(), T_start.begin(), T_start.end());
    T_sol.insert(T_sol.end(), T_sol_tmp.begin(), T_sol_tmp.end());
    T_sol.insert(T_sol.end(), T_end.begin(), T_end.end());

    // concatenate la(xi) vals
    la_sol.insert(la_sol.end(), la_start.begin(), la_start.end());
    la_sol.insert(la_sol.end(), la_sol_tmp.begin(), la_sol_tmp.end());
    la_sol.insert(la_sol.end(), la_end.begin(), la_end.end());    

    // reformat from backwards integration
    std::reverse(xi_sol.begin(), xi_sol.end());
    std::reverse(v_sol.begin(), v_sol.end());
    std::reverse(w_sol.begin(), w_sol.end());
    std::reverse(T_sol.begin(), T_sol.end());
    std::reverse(la_sol.begin(), la_sol.end());
    
    return {xi_sol, v_sol, w_sol, T_sol, la_sol};
}

/*******************************************************************************/

#ifdef ENABLE_MATPLOTLIB
void plot_profiles(const FluidProfile& fp_bag, const FluidProfile& fp_munu, const FluidProfile& fp_veff, const std::string& filename, const double xi_min, const double xi_max) {
    const auto xi_bag = fp_bag.xi_vals();
    const auto v_bag = fp_bag.v_vals();
    const auto w_bag = fp_bag.w_vals();
    const auto T_bag = fp_bag.T_vals();
    const auto la_bag = fp_bag.la_vals();

    const auto xi_munu = fp_munu.xi_vals();
    const auto v_munu = fp_munu.v_vals();
    const auto w_munu = fp_munu.w_vals();
    const auto T_munu = fp_munu.T_vals();
    const auto la_munu = fp_munu.la_vals();

    const auto xi_veff = fp_veff.xi_vals();
    const auto v_veff = fp_veff.v_vals();
    const auto w_veff = fp_veff.w_vals();
    const auto T_veff = fp_veff.T_vals();
    const auto la_veff = fp_veff.la_vals();
    
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

    plt::figure_size(2400, 800);

    // v(xi)
    plt::subplot2grid(2, 2, 0, 0);
    plt::plot(xi_bag, v_bag, opts_bag);
    plt::plot(xi_munu, v_munu, opts_munu);
    plt::plot(xi_veff, v_veff, opts_veff);
    plt::xlabel("xi");
    plt::ylabel("v(xi)");
    plt::xlim(xi_min, xi_max);
    plt::grid(true);
    plt::legend();

    // w(xi)
    plt::subplot2grid(2, 2, 0, 1);
    plt::plot(xi_bag, w_bag, opts_bag);
    plt::plot(xi_munu, w_munu, opts_munu);
    plt::plot(xi_veff, w_veff, opts_veff);
    plt::xlabel("xi");
    plt::ylabel("w(xi)");
    plt::xlim(xi_min, xi_max);
    plt::grid(true);

    // T(xi)
    plt::subplot2grid(2, 2, 1, 0);
    plt::plot(xi_bag, T_bag, opts_bag);
    plt::plot(xi_munu, T_munu, opts_munu);
    plt::plot(xi_veff, T_veff, opts_veff);
    plt::xlabel("xi");
    plt::ylabel("T(xi)");
    plt::xlim(xi_min, xi_max);
    plt::grid(true);

    // la(xi)
    plt::subplot2grid(2, 2, 1, 1);
    plt::plot(xi_bag, la_bag, opts_bag);
    plt::plot(xi_munu, la_munu, opts_munu);
    plt::plot(xi_veff, la_veff, opts_veff);
    plt::xlabel("xi");
    plt::ylabel("la(xi)");
    plt::xlim(xi_min, xi_max);
    plt::grid(true);

    // plt::suptitle("vw = " + to_string_with_precision(vw) + ", alpha = " + to_string_with_precision(alN));

    plt::save(filename);
    std::cout << "Fluid profile saved to " << filename << "\n";
}
#endif

} // namespace Hydrodynamics

// William Searle is a goofy goober and he will never find this message