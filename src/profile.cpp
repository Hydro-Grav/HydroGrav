// profile.cpp
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
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
#include "hydrodynamics.hpp"
#include "physics.hpp"
#include "maths_ops.hpp"

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
- remove dydxi stuff (used only in deflag part of calculation) - not critical, just cleans up code
- update get_mode for veff
- fix deflag solver for veff
- get_TmTN(wmwN) uses bag model currently (mu nu giving bad results)
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

double dvdxi(double xi, double v, const double csq) {
    const auto mu_val = mu(xi, abs(v));
    const auto denom = gammaSq(v) * (1.0 - v * xi) * (mu_val * mu_val / csq - 1.0);
    return (2.0 * v / xi) / denom;
}

double dwdxi(double xi, double v, double w, const double csq) {
    return w * gammaSq(v) * mu(xi, v) * (1.0 + 1.0 / csq) * dvdxi(xi, v, csq);
}

double dTdxi(double xi, double v, double T, const double csq) {
    return T * gammaSq(v) * mu(xi, v) * dvdxi(xi, v, csq);
}

state_type dydxi_vec(double xi, const state_type& y, double vw, double cpsq, double cmsq) {
    const auto v = y[0];
    const auto w = y[1];
    const auto T = y[2];
    const auto csq = (xi < vw) ? cmsq : cpsq;

    return { dvdxi(xi, v, csq), dwdxi(xi, v, w, csq), dTdxi(xi, v, T, csq) };
}

double dxidv(double xi, double v, const double csq) {
    const auto mu_val = mu(xi, abs(v));
    if (v < 1e-9) return 0.0; // dxidv -> 0 for (xi,v)->(csq,0)

    return gammaSq(v) * (1.0 - v * xi) * (mu_val * mu_val / csq - 1.0) * xi / (2.0 * v);
}

double dwdv(double xi, double v, double w, const double csq) {
    return w * gammaSq(v) * mu(xi, v) * (1.0 + 1.0 / csq);
}

double dTdv(double xi, double v, double T, const double csq) {
    return T * gammaSq(v) * mu(xi, v);
}

state_type dydv_vec(double v, const state_type& y, double vw, double cmsq, double cpsq) {
    const auto xi = y[0];
    const auto w = y[1];
    const auto T = y[2];
    const auto csq = (xi < vw) ? cmsq : cpsq;

    return { dxidv(xi, v, csq), dwdv(xi, v, w, csq), dTdv(xi, v, T, csq) };
}
/*************************************************************************************/

// Warning: doesn't work for w profile yet!
void generate_streamplot_data(const PhaseTransition::PTParams& params, int xi_pts, int y_pts, const std::string& filename) {
    // EoM parametrised by time coord tau
    auto dxi_dtau = [] (double xi, double v, const double csq) {
        return xi * ((xi-v)*(xi-v) - csq * (1-xi*v)*(1-xi*v));
    };

    auto dv_dtau = [] (double xi, double v, const double csq) {
        return 2.0 * v * csq * (1-v*v) * (1 - xi*v);
    };

    auto dw_dtau = [&dv_dtau] (double xi, double v, double w, const double csq) {
        return w * (1 + 1/csq) * gammaSq(v) * mu(xi, v) * dv_dtau(xi, v, csq);
    };
    
    std::cout << "Generating streamplot data for fluid profile... ";
    std::cout << "(warning: does not work for w(xi) profile yet!) ";
    
    std::ofstream file(filename);
    file << "xi,v,w,dxidtau,dvdtau,dwdtau\n";
    file << std::fixed << std::setprecision(8); // needed for compatibility with python streamplot

    // Define grid ranges (avoid's singularity at xi=0)
    const double xi_min = 0.01;
    const double xi_max = 0.99;
    const double y_min = 0.01; // bounds for v, w the same
    const double y_max = 0.99;

    // create grid for streamplot
    const auto xi_vals = linspace(xi_min, xi_max, xi_pts);
    const auto y_vals = linspace(y_min, y_max, y_pts);

    for (double xi : xi_vals) {
        const auto csq = (xi < params.vw()) ? params.cmsq() : params.cpsq();
        for (double y : y_vals) {
            // Avoid division by zero
            if (std::abs(1 - xi * y) < 1e-6) continue;

            const auto dxi = dxi_dtau(xi, y, csq);
            const auto dv = dv_dtau(xi, y, csq);
            const auto dw = dv; // UPDATE THIS -> need dw at fixed v, but what v?

            file << xi << "," << y << "," << y << "," << dxi << "," << dv << "," << dw << "\n";
        }
    }

    file.close();
    std::cout << "Streamplot data saved to " << filename << "\n";
    
    return;
}

void generate_streamplot_data(const PhaseTransition::PTParams& params) {
    generate_streamplot_data(params, 30, 30, "streamplot_data.csv");
    return;
}

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

FluidProfile::FluidProfile(const PhaseTransition::PTParams& params, const size_t n)
    : params_(&params),
      cpsq_(params.cpsq()), cmsq_(params.cmsq()),
      vw_(params.vw()), alN_(params.alN()),
      alp_min_(std::numeric_limits<double>::quiet_NaN()), 
      alp_max_(std::numeric_limits<double>::quiet_NaN()),
      mode_(),
      xi_vals_(), v_vals_(), w_vals_(), T_vals_(), la_vals_()
    {
        std::vector<prof_type> profiles;

        // define hydrodynamic mode
        mode_ = get_mode_bag(vw_, cmsq_, alN_);

        // check alN large enough for shock (deflag/hybrid only)
        if (mode_ == 0 || mode_ == 1) {
            const auto alp_minmax = get_alp_minmax(vw_);
            alp_min_ = alp_minmax[0];
            alp_max_ = alp_minmax[1];

            // alN > alp > alp_min (can't properly constrain from above since we need alp)
            if (alN_ <= alp_min_) throw std::invalid_argument("alN too small for shock!");
        }

        // calculate fluid profiles v(xi), w(xi), la(xi)
        switch (params.eos()) {
            case PhaseTransition::PTParams::ModelType::Bag:
                if (cpsq_ == cmsq_) { // bag
                    std::cout << "Calculating fluid profile using Bag equation of state\n";
                } else { // mu nu
                    std::cout << "Calculating fluid profile using modified Bag equation of state\n";
                }
                bag_params_ = &dynamic_cast<const PhaseTransition::PTParams_Bag&>(params);
                profiles = solve_profile(n);
                break;
            case PhaseTransition::PTParams::ModelType::Veff:
                std::cout << "Calculating fluid profile using generic equation of state from Veff\n";
                std::cout << "Warning: alN stored in PTParams is not used for Veff EoS!\n";
                veff_params_ = &dynamic_cast<const PhaseTransition::PTParams_Veff&>(params);
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

    std::cout << "Fluid profile saved to " << filename << "!\n";

    return;
}

#ifdef ENABLE_MATPLOTLIB
void FluidProfile::plot(const std::string& filename) const {
    plt::figure_size(2400, 800);

    // v(xi)
    plt::subplot2grid(2, 2, 0, 0);
    plt::plot(xi_vals_, v_vals_);
    plt::xlabel("xi");
    plt::ylabel("v(xi)");
    plt::xlim(0.0, 1.0);
    plt::grid(true);

    // w(xi)
    plt::subplot2grid(2, 2, 0, 1);
    plt::plot(xi_vals_, w_vals_);
    plt::xlabel("xi");
    plt::ylabel("w(xi)");
    plt::xlim(0.0, 1.0);
    plt::grid(true);

    // T(xi)
    plt::subplot2grid(2, 2, 1, 0);
    plt::plot(xi_vals_, T_vals_);
    plt::xlabel("xi");
    plt::ylabel("T(xi)");
    plt::xlim(0.0, 1.0);
    plt::grid(true);

    // la(xi)
    plt::subplot2grid(2, 2, 1, 1);
    plt::plot(xi_vals_, la_vals_);
    plt::xlabel("xi");
    plt::ylabel("la(xi)");
    plt::xlim(0.0, 1.0);
    plt::grid(true);

    plt::suptitle("vw = " + to_string_with_precision(vw_) + ", alpha = " + to_string_with_precision(alN_));
    plt::save(filename);

    std::cout << "Fluid profile plot saved to '" << filename << "'." << std::endl;

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

// need to check sign convention - see arXiv:1909.10040 eq B.6
double FluidProfile::vp_from_matching(double vm, double alp) const { // vp(vm,alp) from matching eqs
    const auto sgn = 1.0;

    const auto fac1 = 1.0 / (2.0 * (1.0 / (3.0 * cmsq_) + alp));
    const auto fac2 = 1.0 / (3.0 * vm);
    const auto fac3 = fac2 - vm / (3.0 * cmsq_);

    return fac1 * (fac2 + vm / (3.0 * cmsq_) + sgn * std::sqrt(fac3 * fac3 + 4.0 * alp * alp + 4.0 * (1.0/cmsq_ - 1.0) * alp / 3.0));
}

// need to check sign convention - see arXiv:1909.10040 eq B.7
double FluidProfile::vm_from_matching(double vp, double alp) const { // inverse of vp(vm,alp)
    const auto vp_abs = abs(vp);
    const auto sgn = 1.0; // not sure when to use which sign
    const auto fac = vp_abs + cmsq_ * (1.0 - 3.0 * alp * (1.0 - vp_abs * vp_abs)) / vp_abs;

    return 0.5 * (fac + sgn * std::sqrt(fac * fac - 4.0 * cmsq_)); // mu nu
}

double FluidProfile::w1wN_from_matching(double xi_sh) const { // w1/wN
    // alpha_1 w1 = alpha_N wN
    const auto xi_sh_sq = xi_sh * xi_sh;
    return (xi_sh_sq - cpsq_ * cpsq_) / (cpsq_ * (1.0 - xi_sh_sq)); // mu nu
}

double FluidProfile::get_T1TN(double w1wN) const {
    // const auto mu = 1.0 + 1.0 / cpsq_;
    // return std::pow(w1wN, 1.0 / mu);
    return std::pow(w1wN, 0.25);
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

double FluidProfile::v1UF_from_shock(double xi_sh) const {
    // shock condition (cp < xi_sh < 1) relaxed here since for some vw & alN, xi_shock VERY close to bounds
    // so root-finder takes xi_sh_min = cp, xi_sh_max = 1
    if (xi_sh < std::sqrt(cpsq_) || xi_sh > 1.0)
        throw std::invalid_argument("shock must be supersonic and less than speed of light (cp < xi_sh < 1)");

    if (xi_sh == std::sqrt(cpsq_)) return 1e-10; // avoid numerical precision errors
    return (xi_sh * xi_sh - cpsq_) / ((1.0 - cpsq_) * xi_sh);
}

std::array<double, 2> FluidProfile::get_alp_minmax(double vw) const {
    const auto cp = std::sqrt(cpsq_);

    const auto vm = std::min(cp, vw); // vw for deflag (vw < cm), cm for hybrid (cp < vw)
    const auto vp_min = 0.0;
    const auto vp_max = vm; // |v+| < |v-|
    
    // same as get_alp_wall but using vp, vm
    auto get_alp = [this] (double vp, double vm) {
        return gammaSq(vp) * (vp * vp / cmsq_ - vp * vm / cmsq_ - vp / vm + 1.0) / 3.0; // mu nu 
    };
    
    const auto al_max = get_alp(vp_min, vm);
    const auto al_min = get_alp(vp_max, vm);

    if (al_min < 0.0) return {0.0, al_max}; // numerical precision issue

    return {al_min, al_max};
}

// alpha_+ from wall condition
double FluidProfile::get_alp_wall(double vpUF, double vw) const {
    return gammaSq(vpUF) * vpUF * (1.0 + (1.0 / cmsq_ - 1.0) * vw * vpUF - vw * vw / cmsq_) / (3.0 * vw); // mu nu
}

double FluidProfile::alN_residual_func(double xi_sh, const deriv_func& dydxi, const int n) const {
    // initial conditions
    const auto xi0 = xi_sh - 0.001;
    const auto xif = vw_ + 0.001;

    const auto v1UF = v1UF_from_shock(xi_sh);
    const auto w1wN = w1wN_from_matching(xi_sh);

    const std::array<double, 3> y0 = {v1UF, w1wN, 0.0}; // v0 = v(xi_sh) = v1UF

    // solve fluid EoM to get vpUF
    // WARNING: choosing num steps too small gives bad result!
    const auto [xi_sol, y_sol] = rk4_solver(dydxi, xi0, xif, y0, n);
    const auto vpUF = y_sol.back()[0]; // vpUF = v(xi_w) (endpoint of integration)
    const auto wpwN = y_sol.back()[1];
    
    // calc alpha_N from wall constraint
    const auto alN_wall = wpwN * get_alp_wall(vpUF, vw_);

    // return alN_wall - alN_;
    return std::log(std::abs(alN_wall / alN_)); // doesn't always work for some vw, alN
}

double FluidProfile::lambda_b(double wowN) const {
    // la(xi) behind bubble wall (detonations)
    return (wowN - (1.0 + 3.0 * cmsq_ * alN_)) / (1.0 + cmsq_);
}

double FluidProfile::lambda_s(double wowN) const {
    // la(xi) in front of bubble wall (deflagrations)
    return (wowN - 1.0) / (1.0 + cpsq_);
}

double FluidProfile::find_shock(const deriv_func& dydxi) const {
    // Root-finding algorithm for initial condition v0 = v(xi_sh) = v1UF
    
    // f(xi_sh) = alN_calc - alN_actual
    auto residual = [this, &dydxi] (double xi_sh) {
        return alN_residual_func(xi_sh, dydxi);
    };

    // cp < xi_sh < 1 (shock must be supersonic and less than speed of light)
    const double xi_sh_min = std::max(std::sqrt(cpsq_), vw_); // xi_sh > cp > xi_w (deflag), xi_sh > xi_w > cp (hybrid)
    const double xi_sh_max = 1.0 - 1e-15; // xi_sh_max = 1 will cause residual to diverge (note: for extreme cases of hybrids, xi-sh can be very close to 1)

    // const size_t n = 500;
    // const auto xi_sh_vals = linspace(xi_sh_min, 0.6, n);
    // std::vector<double> res_vals;
    // for (const auto xi_sh : xi_sh_vals) {
    //     res_vals.push_back(residual(xi_sh));
    // }

    // plt::figure_size(800, 800);
    // plt::plot(xi_sh_vals, res_vals);
    // plt::xlabel("xi_sh");
    // plt::ylabel("residual");
    // plt::xlim(0.4, 0.6);
    // plt::grid(true);
    // plt::save("shock_residual.png");

    return root_finder(residual, xi_sh_min, xi_sh_max, 1e-7, 100);
}

std::pair<double, state_type> FluidProfile::get_IC_deflagration(const deriv_func& dydxi) const {
    const auto xi_sh = find_shock(dydxi);

    const auto v1UF = v1UF_from_shock(xi_sh);
    if (abs(v1UF) >= 1.0) throw std::invalid_argument("Deflagration IC failed: v1UF must be <1!");

    const auto w1wN = w1wN_from_matching(xi_sh);
    if (w1wN <= 0.0) throw std::invalid_argument("Deflagration IC failed: w1wN must be >0!");

    const auto T1TN = get_T1TN(w1wN);
    if (T1TN <= 0.0) throw std::invalid_argument("Deflagration IC failed: T1TN must be >0!");

    const std::array<double, 3> y1 = {v1UF, w1wN, T1TN};
    return {xi_sh, y1};
}

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
    if (abs(vm) >= vw_) throw std::invalid_argument("Detonation IC failed: vm<vw required!");

    const auto vmUF = mu(vw_, abs(vm));
    if (abs(vmUF) >= 1.0) throw std::invalid_argument("Detonation IC failed: vmUF must be <1!");
    
    const auto wmwN = w_from_matching(wpwN, vp, vm);
    if (wmwN <= wpwN) throw std::invalid_argument("Detonation IC failed: wm>wp required!");

    const auto TmTN = get_TmTN(wmwN);
    // if (TmTN <= TpTN) throw std::invalid_argument("Detonation IC failed: Tm>Tp required!");

    std::cout << "IC det: vm=" << vm << ", vmUF=" << vmUF << ", wmwN=" << wmwN << ", TmTN=" << TmTN << "\n";

    const state_type y0 = {xi0, wmwN, TmTN};
    
    return {vmUF, y0};
}

/*****************************************************************************************/
/**************************************** Veff EoS ***************************************/

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

std::array<double, 2> FluidProfile::matching_eqs_shock(double pN, double eN, double v2, double v1, double T1TN) const {    
    if (T1TN < veff_params_->TTN_min() || T1TN > veff_params_->TTN_max()) {
        throw std::out_of_range("T1/TN is called out of bounds for spline!");
    }

    // const auto p2 = veff_params_->ps_val(T2TN); // p_2, e_2
    // const auto e2 = veff_params_->es_val(T2TN);
    const auto p2 = pN;
    const auto e2 = eN;

    const auto p1 = veff_params_->ps_val(T1TN); // p_1, e_1
    const auto e1 = veff_params_->es_val(T1TN);

    // const auto eq1 = v2 * v1 * (e1 - e2) - (p1 - p2);
    // const auto eq2 = v1 * (e1 + p2) - v2 * (e2 + p1);
    const auto eq1 = v1 * v2 - (p1 - p2) / (e1 - e2);
    const auto eq2 = v1 / v2 - (e2 + p1) / (e1 + p2);

    return {eq1, eq2};
}

std::array<double, 2> FluidProfile::matching_eqs_shock2(double v1, double T1TN, double v2, double T2TN) const {    
    if (T1TN < veff_params_->TTN_min() || T1TN > veff_params_->TTN_max()) {
        throw std::out_of_range("T1/TN is called out of bounds for spline!");
    }
    if (T2TN < veff_params_->TTN_min() || T2TN > veff_params_->TTN_max()) {
        throw std::out_of_range("T2/TN is called out of bounds for spline!");
    }

    // clamping values leads to derivative=0 in newton_solve_1d, causing problems
    // const auto TTN_min = veff_params_->TTN_min();
    // const auto TTN_max = veff_params_->TTN_max();

    // const auto T1TN_clamped = std::clamp(T1TN, TTN_min, TTN_max);
    // const auto T2TN_clamped = std::clamp(T2TN, TTN_min, TTN_max);

    const auto p1 = veff_params_->ps_val(T1TN); // p_1, e_1
    const auto e1 = veff_params_->es_val(T1TN);

    const auto p2 = veff_params_->ps_val(T2TN); // p_2, e_2
    const auto e2 = veff_params_->es_val(T2TN);

    const auto eq1 = v1 * v2 - (p1 - p2) / (e1 - e2);
    const auto eq2 = v1 / v2 - (e2 + p1) / (e1 + p2);

    return {eq1, eq2};
}

std::array<double, 2> FluidProfile::matching_eqs_wall(double vp, double TpTN, double vm, double TmTN) const {        
    if (TmTN < veff_params_->TTN_min() || TmTN > veff_params_->TTN_max()) {
        throw std::out_of_range("Tm/TN is called out of bounds for spline!");
    }

    const auto pp = veff_params_->ps_val(TpTN); // p_+, e_+
    const auto ep = veff_params_->es_val(TpTN);

    const auto pm = veff_params_->pb_val(TmTN); // p_-, e_-
    const auto em = veff_params_->eb_val(TmTN);

    // using regular form seems to make detonation IC calculation go out of bounds for Tm
    const auto eq1 = vp * vm * (em - ep) - (pm - pp);
    const auto eq2 = vm * (em + pp) - vp * (ep + pm);
    // const auto eq1 = vp * vm - (pm - pp) / (em - ep);
    // const auto eq2 = vm / vp - (ep + pm) / (em + pp);

    return {eq1, eq2};
}

// testing purposes
void FluidProfile::test_residual_veff(const deriv_func& dydv, const size_t n) const {
    std::cout << "Running test for T2/TN residual... ";

    const auto TmTN_vals = linspace(veff_params_->TTN_min(), veff_params_->TTN_max(), n);
    // const auto TmTN_vals = linspace(1.0, 1.1, n);

    std::vector<double> resi_vals(n);
    for (int i = 0; i < n; i++) {
        double resi = std::numeric_limits<double>::quiet_NaN();
        try {
            resi = T2TN_residual(dydv, TmTN_vals[i]);
            // std::cout << "resi=" << resi << " for Tm/TN=" << TmTN_vals[i] << "\n";
        } catch (std::exception& e) {
            // std::cout << "failed for Tm/TN=" << TmTN_vals[i] << "\n";
            // std::cout << e.what() << " for Tm/TN=" << TmTN_vals[i] << "\n";
        }

        resi_vals[i] = resi;
    }
    
    #ifdef ENABLE_MATPLOTLIB
    plt::figure_size(800, 800);
    plt::plot(TmTN_vals, resi_vals);
    plt::xlabel("TmTN");
    plt::ylabel("residual");
    // plt::xlim(1.0, 1.1);
    plt::grid(true);
    plt::save("T2_resi.png");
    #endif

    std::cout << "Test complete. T2/TN residual saved to 'T2_resi.png'\n";
}

void FluidProfile::test_shock_veff(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol) const {
    std::cout << "Running test for shock residual... ";
    const auto pN = veff_params_->pN();
    const auto eN = veff_params_->eN();

    std::vector<double> xi_vals, resi_vals;
    for (int i = 0; i < v_sol.size(); i++) {
        const auto xi_sh = y_sol[i][0];
        if (xi_sh <= vw_) continue;

        xi_vals.push_back(xi_sh);
        
        const auto v1UF = v_sol[i];
        const auto v1 = mu(xi_sh, abs(v1UF)); // v1=mu(xi_sh, v1UF)
        const auto T1TN = y_sol[i][2];

        const auto v2 = xi_sh; // v2=xi_sh

        // shock condition mu(xi_sh, v(xi_sh)) xi_sh = (p1-pN)/(e1-eN)
        resi_vals.push_back(abs(matching_eqs_shock(pN, eN, v2, v1, T1TN)[0]));

        // std::cout << "xi_sh=" << xi_sh << ", resi=" << resi_vals[i] << "\n";
    }

    #ifdef ENABLE_MATPLOTLIB
    plt::figure_size(800, 800);
    plt::plot(xi_vals, resi_vals);
    plt::xlabel("xi_sh");
    plt::ylabel("residual");
    // plt::xlim(1.0, 1.1);
    plt::grid(true);
    plt::save("shock_resi.png");
    #endif

    std::cout << "Test complete. Shock residual saved to 'shock_resi.png'\n";
}

double FluidProfile::find_TmTN_veff(const deriv_func& dydv) const {
    auto safe_residual = [this, &dydv] (double TmTN) {
        try {
            const auto resi = T2TN_residual(dydv, TmTN);
            if (!std::isfinite(resi)) return 1e+5;
            return resi;
        } catch (std::exception& e) {
            return 1e+5;
        }
    };

    // testing residual
    // test_residual_veff(dydv, 500);

    // find bracket where residual is defined and minimum lies
    const auto TmTN_min = veff_params_->TTN_min();
    const auto TmTN_max = veff_params_->TTN_max();
    const auto bracket = find_bracket(safe_residual, TmTN_min, TmTN_max);
    if (!std::isfinite(bracket[0]) || !std::isfinite(bracket[1]))
        throw std::runtime_error("find_TmTN_veff failed (bracket not found)!");

    // find TmTN that minimise residual
    const auto TmTN = golden_section_minimize(safe_residual, bracket[0], bracket[1]);
    // std::cout << "TmTN=" << TmTN << ", bracket=[" << bracket[0] << ", " << bracket[1] << "]\n";
    if (TmTN < 0.0) throw std::runtime_error("find_TmTN_veff failed (TmTN < 0)!");

    return TmTN;
}

double FluidProfile::T2TN_residual(const deriv_func& dydv, double TmTN, const size_t n) const {
    const auto vm = (mode_ == 0) ? vw_ : std::sqrt(veff_params_->csq_b(TmTN)); // vm=vw (deflagration), vm=cm (hybrid)
    const auto wmwN = veff_params_->wb_val(TmTN) / veff_params_->wN();

    // calculate fluid profile from wall to shock
    const auto [v_sol, y_sol] = deflagration_profile(dydv, vm, wmwN, TmTN, false, n);

    // fluid behind shock
    const auto xi_sh = y_sol.back()[0];
    const auto v1UF = v_sol.back();
    const auto v1 = mu(xi_sh, abs(v1UF));
    // const auto w1wN = y_sol.back()[1];
    const auto T1TN = y_sol.back()[2];   

    // matching at shock to get T2/TN
    auto shock_matching_helper = [this, v1, T1TN, xi_sh] (double T2TN_guess) {
        return matching_eqs_shock2(v1, T1TN, xi_sh, T2TN_guess)[1]; // v2 = xi_sh
    };

    const auto T2TN = newton_solve_1d(shock_matching_helper, 1.0);

    return abs(T2TN - 1.0);
}

std::pair<std::vector<double>, std::vector<state_type>> FluidProfile::deflagration_profile(const deriv_func& dydv, double vm, double wmwN, double TmTN, const bool test_shock, const size_t n) const {
    // matching at wall: vm,Tm -> vp,Tp
    auto wall_matching_helper = [this, vm, TmTN] (std::array<double, 2>& yp_guess) {
        return matching_eqs_wall(yp_guess[0], yp_guess[1], vm, TmTN); // yp = {vp, TpTN}
    };

    const std::array<double, 2> yp_guess = {0.9 * vw_, 1.1 * TmTN}; // {vp, Tp}
    const auto sol = newton_solve_2d(wall_matching_helper, yp_guess);

    // fluid in front of wall
    const auto vp = sol[0];
    const auto vpUF = mu(vw_, abs(vp));

    const auto wpwN = w_from_matching(wmwN, vm, vp);
    const auto TpTN = sol[1];

    const state_type yp = {vw_, wpwN, TpTN}; // {xi_w, wpwN, TpTN}

    // solve EoM - integrate from v(xi)=vpUF -> v(xi)=0
    const auto [v_sol_tmp, y_sol_tmp] = rk4_solver(dydv, vpUF, 1e-10, yp, n);

    // find shock & truncate solution
    const auto sh_idx = std::min(find_shock_idx(v_sol_tmp, y_sol_tmp, test_shock), v_sol_tmp.size() - 1);
    const std::vector<double> v_sol(v_sol_tmp.begin(), v_sol_tmp.begin() + sh_idx + 1);
    const std::vector<state_type> y_sol(y_sol_tmp.begin(), y_sol_tmp.begin() + sh_idx + 1);

    return {v_sol, y_sol};
}

// TO DO:
// check resi is sufficiently small
size_t FluidProfile::find_shock_idx(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol, const bool test_shock, const double tol) const {
    const auto pN = veff_params_->pN();
    const auto eN = veff_params_->eN();

    // testing purposes
    // if (test_shock) {
    //     test_shock_veff(v_sol, y_sol);
    // }

    std::vector<double> resi_vals;
    for (int i = 0; i < v_sol.size(); i++) {
        const auto xi_sh = y_sol[i][0];
        if (xi_sh <= vw_) continue;
        
        const auto v1UF = v_sol[i];
        const auto v1 = mu(xi_sh, abs(v1UF)); // v1=mu(xi_sh, v1UF)
        const auto T1TN = y_sol[i][2];

        const auto v2 = xi_sh; // v2=xi_sh

        // shock condition mu(xi_sh, v(xi_sh)) xi_sh = (p1-pN)/(e1-eN)
        resi_vals.push_back(abs(matching_eqs_shock(pN, eN, v2, v1, T1TN)[0]));
    }

    if (resi_vals.size() == 0)
        throw std::runtime_error("find_shock_idx failed (no shock found in fluid profile)!");

    // index where residual is minimised
    const auto it = std::min_element(resi_vals.begin(), resi_vals.end());
    const auto idx = std::distance(resi_vals.begin(), it);

    if (test_shock && resi_vals[idx] > tol) {
        std::cout << "Warning: Shock residual above tolerance (" << tol << "). Residual=" << resi_vals[idx] << "\n";
    }

    return idx;
}

// TO DO: update error handling
std::pair<double, state_type> FluidProfile::get_IC_detonation_veff() const {
    // uses matching conditions to get vm, Tm/TN from vp, Tp/TN
    //      vp * vm = (pm - pp) / (em - ep)
    //      vp / vm = (em + pp) / (ep + pm)

    const auto xi0 = vw_;
    const auto vp = -vw_;
    const auto TpTN = 1.0;
    const auto wpwN = 1.0;

    std::function<std::array<double, 2>(std::array<double, 2>)> matching_helper = [this, vp, TpTN] (std::array<double, 2> ym_guess) {
        return matching_eqs_wall(vp, TpTN, ym_guess[0], ym_guess[1]); // ym_guess = {vm_guess, TmTN_guess}
    };
    
    // solve matching eqs
    const std::array<double, 2> vm_TmTN_guess = {vp, TpTN}; // {vm, TmTN}
    const auto sol = newton_solve_2d(matching_helper, vm_TmTN_guess);

    // fluid behind wall
    const auto vm = sol[0];
    if (abs(vm) >= abs(vp)) throw std::invalid_argument("|vm| must be < |vp|=vw for detonation");

    const auto vmUF = mu(vw_, abs(vm));
    if (vmUF < 0.0) throw std::invalid_argument("vmUF must be > 0 for detonation");

    const auto TmTN = sol[1];
    if (TmTN <= TpTN) throw std::invalid_argument("Tm/TN must be > Tp/TN for detonation");

    const auto wmwN = w_from_matching(wpwN, vp, vm);
    if (wmwN <= wpwN) throw std::invalid_argument("wm/wN must be > wp/wN for detonation");

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
            throw std::invalid_argument("Hydrodynamic mode must be: 0 (deflagration), 1 (hybrid) or 2 (detonation)");
    }

    std::cout << "Solving fluid profile for hydrodynamic mode=";
    if (mode_ == 0) {
        std::cout << "deflagration";
    } else if (mode_ == 1) {
        std::cout << "hybrid";
    } else {
        std::cout << "detonation";
    }
    std::cout << "\n";

    // wrapper for hydrodynamic EoM
    auto dydxi = [this] (double xi, const state_type& y) -> state_type {
        return dydxi_vec(xi, y, vw_, cpsq_, cmsq_);
    };

    // uncomment this if det/hyb fluid profile looks wrong/discontinuous
    auto dydv = [this] (double v, const state_type& y) -> state_type {
        return dydv_vec(v, y, vw_, cpsq_, cmsq_);
    };

    double xi0, xif;
    const auto dlt = 0.001; // wall and shocks are discontinuities so start integration just before them
    std::vector<state_type> y_sol_tmp;
    state_type y0; 
    prof_type xi_sol_tmp, v_sol_tmp, w_sol_tmp, T_sol_tmp, la_sol_tmp;

    const auto w_start_val = 1.0; // w+/wN (det), w2/wN (deflag/hybrid)
    const auto T_start_val = 1.0; // T+/TN (det), T2/TN (deflag/hybrid);
    const auto la_start_val = 0.0;
    
    double w_end_val, T_end_val, la_end_val;

    if (mode_ < 2) { // deflagration & hybrid
        // hybrid and deflagration ICs the same for xi_w < xi < xi_sh
        // v(xi_sh) = v1UF, w(xi_sh) = w1wN, T(xi_sh) = T1TN

        const auto ics = get_IC_deflagration(dydxi);
        const auto xi_sh = ics.first;
        y0 = ics.second;

        const auto v1UF = y0[0];
        const auto w1wN = y0[1];
        const auto T1TN = y0[2];

        // xi0 = xi_sh - dlt;
        // xif = vw_ + dlt;
        xi0 = xi_sh;
        xif = vw_;

        // solver
        const auto sol = rk4_solver(dydxi, xi0, xif, y0, n);
        xi_sol_tmp = sol.first;
        y_sol_tmp = sol.second;

        // fill v, w vectors
        for (size_t i = 0; i < xi_sol_tmp.size(); i++) {
            v_sol_tmp.push_back(y_sol_tmp[i][0]);
            w_sol_tmp.push_back(y_sol_tmp[i][1]);
            T_sol_tmp.push_back(y_sol_tmp[i][2]);

            // std::cout << T_sol_tmp[i] << "\n";

            la_sol_tmp.push_back(lambda_s(w_sol_tmp[i]));
        }

        const auto vpUF = v_sol_tmp.back();
        // if (vpUF >= vw_) throw std::invalid_argument("vpUF must be < vw");

        const auto wpwN = w_sol_tmp.back();
        const auto TpTN = T_sol_tmp.back();
        
        // check alp okay
        // Note: need alp for this so must do root-finding BEFORE determining if alp is good/bad
        const auto alp = get_alp_wall(vpUF, vw_);
        
        if (alp >= alN_) throw std::invalid_argument("alpha_+ must be < alpha_N");
        if (alp < alp_min_) throw std::invalid_argument("alpha_+ too small for shock");
        if (alp > alp_max_) throw std::invalid_argument("alpha_+ too large for shock");

        if (mode_ == 0) {
            // fix end-value for enthalpy
            const auto vm = -vw_;
            const auto vp = vp_from_matching(vm, alp);
            const auto wpwN = w_sol_tmp.back(); // w(xi_w + dlt) = w+/wN

            w_end_val = w_from_matching(wpwN, vp, vm); // wm/wN, from matching condition at wall
            T_end_val = get_TmTN(w_end_val);
            la_end_val = lambda_b(w_end_val);

            std::cout << "Deflagration profile:\n"
                      << "  vm = " << vm << ", vmUF = " << mu(vw_, abs(vm)) << "\n"
                      << "  wmwN = " << w_end_val << ", TmTN = " << T_end_val << "\n"
                      << "  vp = " << vp << ", vpUF = " << vpUF << "\n"
                      << "  wpwN = " << wpwN << ", TpTN = " << T_sol_tmp.back() << "\n"
                      << "  v1 = " << mu(xi_sh, abs(y0[0])) << ", v1UF = " << y0[0] << "\n"
                      << "  w1wN = " << y0[1] << ", T1TN = " << y0[2] << "\n"
                      << "  xi_sh = " << xi_sh << "\n";
                      

        } else { // hybrid
            // consistency check
            if (vw_ >= vJ_det(alp)) {
                std::cout << "Warning: Hybrid condition violated (vw >= vJ_det). Hydrodynamic mode may be incorrect!\n";
                // throw std::runtime_error("Fluid profile failed due to edge case: Incorrect hydrodynamic mode detected!");
            }
            
            // initial conditions for rarefaction wave
            const auto xi0_rf = vw_;
            const auto xif_rf = std::sqrt(cmsq_);

            const auto vm = -std::sqrt(cmsq_);
            const auto vmUF = mu(vw_, abs(vm));

            const auto vp = mu(vw_, abs(vpUF));
            const auto wmwN = w_from_matching(wpwN, vp, vm);

            const auto TmTN = get_TmTN(wmwN);

            const state_type y0_rf = {xi0_rf, wmwN, TmTN};
            const auto [v_sol_rf_tmp, y_sol_rf_tmp] = rk4_solver(dydv, vmUF, 1e-10, y0_rf, n);

            // combine rarefaction wave with shockwave part of solution
            for (size_t i = 0; i < v_sol_rf_tmp.size(); i++) {
                v_sol_tmp.push_back(v_sol_rf_tmp[i]);
                xi_sol_tmp.push_back(y_sol_rf_tmp[i][0]);
                w_sol_tmp.push_back(y_sol_rf_tmp[i][1]);
                T_sol_tmp.push_back(y_sol_rf_tmp[i][2]);

                la_sol_tmp.push_back(lambda_b(y_sol_rf_tmp[i][1]));
            }

            xif = xif_rf; // update xif value to behind rarefaction wave
            w_end_val = w_sol_tmp.back();
            T_end_val = T_sol_tmp.back();
            la_end_val = la_sol_tmp.back();

            std::cout << "Hybrid profile:\n"
                      << "  vm = " << vm << ", vmUF=" << mu(vw_, abs(vm)) << "\n"
                      << "  wmwN = " << w_end_val << ", TmTN = " << TmTN << "\n"
                      << "  vp = " << vp << ", vpUF = " << vpUF << "\n"
                      << "  wpwN = " << wpwN << ", TpTN = " << TpTN << "\n"
                      << "  v1 = " << mu(xi_sh, abs(v1UF)) << ", v1UF = " << v1UF << "\n"
                      << "  w1wN = " << w1wN << ", T1TN = " << T1TN << "\n"
                      << "  xi_sh = " << xi_sh << "\n"
                      << "  w_end = " << w_end_val << ", T_end = " << T_end_val << "\n";
        }

    } else { // detonation
        // cm < xi < xi_w
        xi0 = vw_;
        xif = std::sqrt(cmsq_);

        // solver
        const auto ics = get_IC_detonation();
        const auto vmUF = ics.first;
        const auto y0 = ics.second; // {xi0, w0, T0}

        const auto sol = rk4_solver(dydv, vmUF, 1e-10, y0, n);
        v_sol_tmp = sol.first;
        y_sol_tmp = sol.second;

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

        std::cout << "Detonation profile:\n"
                      << "  vm = " << mu(vw_, abs(vmUF)) << ", vmUF=" << vmUF << "\n"
                      << "  wmwN = " << y0[1] << ", TmTN = " << y0[2] << "\n"
                      << "  w_end = " << w_end_val << ", T_end = " << T_end_val << "\n";
    }
    
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
        std::cout << "deflagration";
    } else if (mode_ == 1) {
        std::cout << "hybrid";
    } else if (mode_ == 2){
        std::cout << "detonation";
    } else {
        throw std::invalid_argument("Hydrodynamic mode must be: 0 (deflagration), 1 (hybrid) or 2 (detonation)");
    }
    std::cout << "\n";

    // wrapper for hydrodynamic EoM
    auto dydv = [this] (double v, const state_type& y) -> state_type {
        return dydv_vec(v, y, vw_, cpsq_, cmsq_);
    };

    const auto eN = veff_params_->eN();
    const auto wN_inv = 1.0 / veff_params_->wN();

    double xi0, xif;
    std::vector<state_type> y_sol_tmp;
    state_type y0; 
    prof_type xi_sol_tmp, v_sol_tmp, w_sol_tmp, T_sol_tmp, la_sol_tmp;

    const auto w_start_val = 1.0; // w+/wN (det), w2/wN (deflag/hybrid)
    const auto T_start_val = 1.0; // T+/TN (det), T2/TN (deflag/hybrid);
    const auto la_start_val = 0.0;
    
    double w_end_val, T_end_val, la_end_val;

    if (mode_ < 2) { // deflagration & hybrid
        // hybrid and deflagration ICs the same for xi_w < xi < xi_sh
        // v(xi_sh) = v1UF, w(xi_sh) = w1wN, T(xi_sh) = T1TN

        const auto TmTN = find_TmTN_veff(dydv);
        const auto vm = (mode_ == 0) ? vw_ : std::sqrt(veff_params_->csq_b(TmTN)); // vm=vw (deflagration), vm=cm (hybrid)
        const auto wmwN = veff_params_->wb_val(TmTN) / veff_params_->wN();

        const auto prof_tmp = deflagration_profile(dydv, vm, wmwN, TmTN, true, n);
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

        if (mode_ == 0) { // deflagration
            xif = vw_;

            // fix end values 
            w_end_val = wmwN;
            T_end_val = TmTN;
            la_end_val = lambda_b_veff(T_end_val, eN, wN_inv); // lambda just behind wall (broken phase)

            std::cout << "Deflagration profile:\n"
                      << "  vm = " << vm << ", vmUF=" << mu(vw_, abs(vm)) << "\n"
                      << "  wmwN = " << wmwN << ", TmTN = " << TmTN << "\n"
                      << "  vp = " << mu(vw_, abs(vpUF)) << ", vpUF = " << vpUF << "\n"
                      << "  wpwN = " << wpwN << ", TpTN = " << T_sol_tmp.back() << "\n"
                      << "  v1 = " << mu(xi0, abs(v_sol_tmp.front())) << ", v1UF = " << v_sol_tmp.front() << "\n"
                      << "  w1wN = " << w_sol_tmp.front() << ", T1TN = " << T_sol_tmp.front() << "\n"
                      << "  xi_sh = " << xi0 << "\n";

        } else { // hybrid
            // initial conditions for rarefaction wave
            const auto xi0_rf = vw_;

            const auto vmUF = mu(vw_, abs(vm));
            const state_type y0 = {xi0_rf, wmwN, TmTN}; // {xi0, wmwN, TmTN}
            const auto [v_sol_rf_tmp, y_sol_rf_tmp] = rk4_solver(dydv, vmUF, 1e-10, y0, n);

            // combine rarefaction wave with shockwave part of solution
            for (size_t i = 0; i < v_sol_rf_tmp.size(); i++) {
                v_sol_tmp.push_back(v_sol_rf_tmp[i]);
                xi_sol_tmp.push_back(y_sol_rf_tmp[i][0]);
                w_sol_tmp.push_back(y_sol_rf_tmp[i][1]);
            
                const auto T_sol_rf = y_sol_rf_tmp[i][2];
                T_sol_tmp.push_back(T_sol_rf);
                la_sol_tmp.push_back(lambda_b_veff(T_sol_rf, eN, wN_inv));
            }

            xif = std::sqrt(veff_params_->csq_b(TmTN)); // cmsq = csq_b(TmTN)
            w_end_val = w_sol_tmp.back();
            T_end_val = T_sol_tmp.back();
            la_end_val = la_sol_tmp.back();

            std::cout << "Hybrid profile:\n"
                      << "  vm = " << vm << ", vmUF=" << vmUF << "\n"
                      << "  wmwN = " << wmwN << ", TmTN = " << TmTN << "\n"
                      << "  vp = " << mu(vw_, abs(vpUF)) << ", vpUF = " << vpUF << "\n"
                      << "  wpwN = " << wpwN << ", TpTN = " << T_sol_tmp.back() << "\n"
                      << "  v1 = " << mu(xi0, abs(v_sol_tmp.front())) << ", v1UF = " << v_sol_tmp.front() << "\n"
                      << "  w1wN = " << w_sol_tmp.front() << ", T1TN = " << T_sol_tmp.front() << "\n"
                      << "  xi_sh = " << xi0 << "\n"
                      << "  w_end = " << w_end_val << ", T_end = " << T_end_val << "\n";
        }

    } else { // detonation
        // cm < xi < xi_w
        // v0 = v(xi_w) = vm(UF)
        // w0 = w(xi_w) = wm/wN
        // T0 = T(xi_w) = Tm/TN

        const auto ics = get_IC_detonation_veff();
        const auto vmUF = ics.first;
        const auto y0 = ics.second; // {xi_w, wmwN, TmTN}

        // solver
        const auto sol = rk4_solver(dydv, vmUF, 1e-10, y0, n);
        v_sol_tmp = sol.first;
        y_sol_tmp = sol.second;

        // fill v, w vectors
        for (size_t i = 0; i < v_sol_tmp.size(); i++) {
            xi_sol_tmp.push_back(y_sol_tmp[i][0]);
            w_sol_tmp.push_back(y_sol_tmp[i][1]);

            const auto T_sol = y_sol_tmp[i][2];
            T_sol_tmp.push_back(T_sol);
            la_sol_tmp.push_back(lambda_b_veff(T_sol, eN, wN_inv));
        }

        xi0 = vw_;
        xif = std::sqrt(veff_params_->csq_b(y0[2])); // cmsq = csq_b(Tm/TN)  

        w_end_val = w_sol_tmp.back();
        T_end_val = T_sol_tmp.back();
        la_end_val = la_sol_tmp.back();

        std::cout << "Detonation profile:\n"
                  << "  vm = " << mu(vw_, abs(vmUF)) << ", vmUF=" << vmUF << "\n"
                  << "  wmwN = " << y0[1] << ", TmTN = " << y0[2] << "\n"
                  << "  w_end = " << w_end_val << ", T_end = " << T_end_val << "\n";
    }

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
    // plt::xlim(0.56, 0.61);
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
    // plt::xlim(0.56, 0.61);
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