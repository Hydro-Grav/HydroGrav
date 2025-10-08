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
- change everything to {v,T,w} format (currently is {v,w,T})
- go through and check if any varaibles can be removed from class (might need to add eN, wN_inv for lambda calc)
- return std::numeric_limits<double>::quiet_NaN(); rather than throwing exception for veff shock finding!! (innappropriate to use exceptions here)
- move dydxi as a global function in profile since it is called so often
- allocate memory to y_sol_tmp etc (inside solve) -> get rid of .push_back() and change to indexing
- add dflt ctor for uninitialised FluidProfile (for testing)
*/

/* NOTES
- detonation part of hybrid is unstable because of IC v(xi=vw)=cs -> mu^2/cs^2 - 1 = 0 -> dvdxi -> inf
- reparameterising in terms of v can (kind of) fix this, although have to integrate from vmUF -> 1e-5 (zero causes dxidv=0 so xi->inf)
- simpler to offset IC slightly (e.g. use vmUF*0.95 rather than vmUF)
- WARNING: This may lead to innacuracies if |vw-cm| >~ |vw-xi_sh| (since det part of hybrid is large)
- if this is an issue then use reparameterisation of eom in terms of v (commented out for now)
- deflagration shock calculation needs to be redone entirely if using reparameterisation for deflagrations too (unnecessary)
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

state_type dydxi_vec(double xi, const state_type& y, double vw, double cmsq, double cpsq) {
    const auto v = y[0];
    const auto w = y[1];
    const auto T = y[2];
    const auto csq = (xi < vw) ? cmsq : cpsq; // change to csq(T) when implemented

    return { dvdxi(xi, v, csq), dwdxi(xi, v, w, csq), dTdxi(xi, v, T, csq) };
}

// double dxidv(double xi, double v, const double csq) {
//     const auto mu_val = mu(xi, abs(v));
//     return gammaSq(v) * (1.0 - v * xi) * (mu_val * mu_val / csq - 1.0) * xi / (2.0 * v);
// }

// double dwdv(double xi, double v, double w, const double csq) {
//     return w * gammaSq(v) * mu(xi, v) * (1.0 + 1.0 / csq);
// }

// double dTdv(double xi, double v, double T, const double csq) {
//     return T * gammaSq(v) * mu(xi, v);
// }

// state_type dydv_vec(double xi, const state_type& y, double vw, double cmsq, double cpsq) {
//     const auto v = y[0];
//     const auto w = y[1];
//     const auto T = y[2];
//     const auto csq = (xi < vw) ? cmsq : cpsq; // change to csq(T) when implemented

//     return { dxidv(xi, v, csq), dwdv(xi, v, w, csq), dTdv(xi, v, T, csq) };
// }
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
    : params_(params),
      cpsq_(params.cpsq()), cmsq_(params.cmsq()),
      cp_(std::sqrt(cpsq_)), cm_(std::sqrt(cmsq_)),
      vw_(params.vw()), alN_(params.alN()),
      alp_min_(std::numeric_limits<double>::quiet_NaN()), 
      alp_max_(std::numeric_limits<double>::quiet_NaN()),
      mode_(),
      xi_vals_(), v_vals_(), w_vals_(), T_vals_(), la_vals_()
    {
        std::vector<prof_type> profiles;

        if (params.eos_model() == "veff") { // veff eos
            std::cout << "Calculating fluid profile using generic equation of state from Veff\n";
            std::cout << "Warning: alN stored in PTParams is not used for Veff EoS!\n";

            // replace with get_mode_veff when implemented!!
            mode_ = get_mode_bag(vw_, cmsq_, alN_); // placeholder for now!!

            // delete when get_mode_veff is implemented
            if (mode_ == 0 || mode_ == 1) {
                const auto alp_minmax = get_alp_minmax(vw_);
                alp_min_ = alp_minmax[0];
                alp_max_ = alp_minmax[1];

                // alN > alp > alp_min (can't properly constrain from above since we need alp)
                if (alN_ <= alp_min_) throw std::invalid_argument("alN too small for shock!");
            }

            // calculate fluid profiles v(xi), w(xi), la(xi) (remove from conditional when solve_profile finished)
            profiles = solve_profile_veff(n);
        } else {
            std::cout << "Calculating fluid profile using Bag equation of state\n";

            // define hydrodynamic mode
            mode_ = get_mode_bag(vw_, cmsq_, alN_);

            // check alN large enough for shock (deflag/hybrid only)
            if (mode_ == 0 || mode_ == 1) {
                const auto alp_minmax = get_alp_minmax(vw_);
                alp_min_ = alp_minmax[0];
                alp_max_ = alp_minmax[1];

                std::cout << "alp_min=" << alp_min_ << ", alp_max=" << alp_max_ << "\n";

                // alN > alp > alp_min (can't properly constrain from above since we need alp)
                if (alN_ <= alp_min_) throw std::invalid_argument("alN too small for shock!");
            }

            // calculate fluid profiles v(xi), w(xi), la(xi) (remove from conditional when solve_profile finished)
            profiles = solve_profile(n);
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
    const auto mu = 1.0 + 1.0 / cpsq_;
    return std::pow(w1wN, 1.0 / mu);
}

// TO DO: update ap/am ratio calculation
double FluidProfile::get_TmTN(double wmwN) const {
    const auto mu = 1.0 + 1.0 / cpsq_;
    const auto nu = 1.0 + 1.0 / cmsq_;
    const auto r = 1.0; // ap/am ratio

    const auto fac = (mu / nu) * r * wmwN;
    return std::pow(fac, 1.0 / nu) * std::pow(params_.TN(), mu / nu - 1.0);
}

// update for mu nu - done
double FluidProfile::v1UF_from_shock(double xi_sh) const {
    if (xi_sh < std::sqrt(cpsq_) || xi_sh > 1.0) {
        // shock condition (cp < xi_sh < 1) relaxed here since for some vw & alN, xi_shock VERY close to bounds
        // so root-finder takes xi_sh_min = cp, xi_sh_max = 1
        throw std::invalid_argument("shock must be supersonic and less than speed of light (cp < xi_sh < 1)");
    }

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
    // const double xi_sh_max = 1.0 - 1e-5; // need to make this closer to 1 for extreme case of hybrids with xi_sh very close to 1
    const double xi_sh_max = 1.0 - 1e-15;

    return root_finder(residual, xi_sh_min, xi_sh_max, 1e-7, 100);
}

std::pair<double, state_type> FluidProfile::get_IC_deflagration(const deriv_func& dydxi) const {
    const auto xi_sh = find_shock(dydxi);

    const auto v1UF = v1UF_from_shock(xi_sh);
    if (abs(v1UF) >= 1.0) throw std::invalid_argument("Deflagration IC failed: v1UF must be <1!");

    const auto w1wN = w1wN_from_matching(xi_sh);
    if (w1wN <= 0.0) throw std::invalid_argument("Deflagration IC failed: w1wN must be >0!");

    // const auto T1TN = ToTN(w1wN, cpsq_, cpsq_); // c1sq = c2sq = cpsq
    const auto T1TN = get_T1TN(w1wN);
    if (T1TN <= 0.0) throw std::invalid_argument("Deflagration IC failed: T1TN must be >0!");

    const std::array<double, 3> y1 = {v1UF, w1wN, T1TN};
    return {xi_sh, y1};
}

state_type FluidProfile::get_IC_detonation() const {
    // initial conditions:
    //     v0 = v(xi_w) = vm(UF)
    //     w0 = w(xi_w) = wm/wN
    //     T0 = T(xi_w) = Tm/TN
    
    const auto vp = -vw_;
    const auto wpwN = 1.0; // w+ = wN
    const auto TpTN = 1.0;

    const auto alp = alN_; // alpha_+ = alpha_N

    const auto vm = vm_from_matching(vp, alp);
    if (abs(vm) >= vw_) throw std::invalid_argument("Detonation IC failed: vm<vw required!");

    const auto vmUF = mu(vw_, abs(vm));
    if (abs(vmUF) >= 1.0) throw std::invalid_argument("Detonation IC failed: vmUF must be <1!");
    
    const auto wmwN = wm_from_matching(wpwN, vp, vm);
    if (wmwN <= wpwN) throw std::invalid_argument("Detonation IC failed: wm>wp required!");

    // const auto TmTN = ToTN(wmwN, cpsq_, cmsq_);
    const auto TmTN = get_TmTN(wmwN);
    if (TmTN <= TpTN) throw std::invalid_argument("Detonation IC failed: Tm>Tp required!");
    
    return {vmUF, wmwN, TmTN};
}

/*****************************************************************************************/
/**************************************** Veff EoS ***************************************/

// change these to lambda functions in solve?
double FluidProfile::lambda_s_veff(double ToTN, const double eN, const double wN_inv) const {
    // la(xi) = (es(T(xi)) - eN) / wN
    const auto es_T = params_.es_val(ToTN); // es(T/TN)
    return (es_T - eN) * wN_inv;
}

double FluidProfile::lambda_b_veff(double ToTN, const double eN, const double wN_inv) const {
    // la(xi) = (eb(T(xi)) - eN) / wN
    const auto eb_T = params_.eb_val(ToTN); // eb(T/TN)
    return (eb_T - eN) * wN_inv;
}

state_type FluidProfile::test_shock_matching(const deriv_func& dydxi, double xi_sh) const {
    // fluid in front of shock
    const auto v2 = xi_sh;
    const auto w2wN = 1.0;
    const auto T2TN = 1.0;

    const auto pN = params_.pN();
    const auto eN = params_.eN();

    // matching across shock
    std::function<std::array<double, 2>(std::array<double, 2>)> shock_matching_helper = [this, v2, pN, eN] (std::array<double, 2> y0) {
        return matching_eqs_shock(pN, eN, v2, y0[0], y0[1]); // y0 = {v1, T1TN}
    };
    
    // guess needs to be slightly larger than {v2, T2TN} since residual goes to zero for v1=v2, T1=T2
    // const auto v1_guess = v2 * (1.0 + 1e-3);
    // const auto T1TN_guess = T2TN * (1.0 + 1e-3);
    const std::array<double, 2> v1_T1TN_guess = {v2 * (1.0 + 1e-3), T2TN * (1.0 + 1e-3)}; // {v1, T1TN}
    const auto shock_sol = newton_solve_2d(shock_matching_helper, v1_T1TN_guess);

    // fluid behind shock
    const auto v1 = shock_sol[0];
    const auto v1UF = mu(xi_sh, abs(v1));
    if (abs(v1UF) >= 1.0 || abs(v1) >= 1.0) throw std::invalid_argument("Shock matching: v1 must be <1!");
    if (v1 >= v2) throw std::invalid_argument("Shock matching: v1 must be less than v2!");

    const auto w1wN = wm_from_matching(w2wN, v2, v1);
    if (w1wN <= 0.0) throw std::invalid_argument("Shock matching: w1wN must be >0!");
    if (w1wN <= w2wN) throw std::invalid_argument("Shock matching: w1wN must be greater than w2wN!");

    const auto T1TN = shock_sol[1];
    if (T1TN <= 0.0) throw std::invalid_argument("Shock matching: T1TN must be >0!");
    if (T1TN <= T2TN) throw std::invalid_argument("Shock matching: T1TN must be greater than T2TN!"); 

    return {v1UF, w1wN, T1TN};
}

std::pair<state_type, state_type> FluidProfile::test_wall_matching(const deriv_func& dydxi, double xi_sh, state_type& y0, const int n) const {
    // solve EoM from shock to wall
    // NOTE: technically integration should start just behind shock and end just in front of wall, but residual for matching eqs is unstable (only 1 point where they both go to zero)
    const auto xi0 = xi_sh;
    const auto xif = vw_;
    const auto [xi_sol, y_sol] = rk4_solver(dydxi, xi0, xif, y0, n);

    // fluid in front of wall
    const auto vpUF = y_sol.back()[0]; // vpUF = v(xi_w) (endpoint of integration)
    const auto vp = mu(vw_, abs(vpUF));
    if (abs(vpUF) >= 1.0 || abs(vp) >= 1.0) throw std::invalid_argument("Wall matching: vp must be <1!");

    const auto wpwN = y_sol.back()[1];
    if (wpwN <= 0.0) throw std::invalid_argument("Wall matching: wpwN must be >0!");

    const auto TpTN = y_sol.back()[2];
    if (TpTN <= 0.0) throw std::invalid_argument("Wall matching: TpTN must be >0!");

    // matching across wall
    std::function<std::array<double, 2>(std::array<double, 2>)> wall_matching_helper = [this, vp, TpTN] (std::array<double, 2> ym) {
        return matching_eqs_wall(vp, TpTN, ym[0], ym[1]); // ym = {vm, TmTN}
    };

    const std::array<double, 2> vm_TmTN_guess = {cm_, 1.0}; // {vm, TmTN}
    const auto wall_sol = newton_solve_2d(wall_matching_helper, vm_TmTN_guess);

    // fluid behind wall
    const auto vm = wall_sol[0];
    if (vm <= vp) throw std::invalid_argument("Wall matching: vm must be greater than vp!");

    const auto vmUF = mu(vw_, abs(vm));
    if (abs(vmUF) >= 1.0) throw std::runtime_error("Wall matching: vmUF > 1!");

    const auto wmwN = wm_from_matching(wpwN, vp, vm);
    if (wmwN <= 0.0) throw std::invalid_argument("Wall matching: wmwN must be >0!");
    if (wmwN >= wpwN) throw std::invalid_argument("Wall matching: wmwN must be less than wpwN!");

    const auto TmTN = wall_sol[1];
    if (TmTN <= 0.0) throw std::invalid_argument("Wall matching: TmTN must be >0!");
    if (TmTN >= TpTN) throw std::invalid_argument("Wall matching: TmTN must be less than TpTN!");

    const state_type yp = {vpUF, wpwN, TpTN};
    const state_type ym = {vmUF, wmwN, TmTN};

    return {yp, ym};
}

std::array<double, 2> FluidProfile::matching_eqs_shock(double pN, double eN, double v2, double v1, double T1TN) const {    
    if (T1TN < params_.TTN_min() || T1TN > params_.TTN_max()) {
        throw std::out_of_range("T1/TN is called out of bounds for spline!");
    }

    // const auto p2 = params_.ps_val(T2TN); // p_2, e_2
    // const auto e2 = params_.es_val(T2TN);
    const auto p2 = pN;
    const auto e2 = eN;

    const auto p1 = params_.ps_val(T1TN); // p_1, e_1
    const auto e1 = params_.es_val(T1TN);

    // const auto eq1 = v2 * v1 * (e1 - e2) - (p1 - p2);
    // const auto eq2 = v1 * (e1 + p2) - v2 * (e2 + p1);
    const auto eq1 = v1 * v2 - (p1 - p2) / (e1 - e2);
    const auto eq2 = v1 / v2 - (e2 + p1) / (e1 + p2);

    return {eq1, eq2};
}

std::array<double, 2> FluidProfile::matching_eqs_wall(double vp, double TpTN, double vm, double TmTN) const {        
    if (TmTN < params_.TTN_min() || TmTN > params_.TTN_max()) {
        throw std::out_of_range("Tm/TN is called out of bounds for spline!");
    }

    const auto pp = params_.ps_val(TpTN); // p_+, e_+
    const auto ep = params_.es_val(TpTN);

    const auto pm = params_.pb_val(TmTN); // p_-, e_-
    const auto em = params_.eb_val(TmTN);

    // using regular form seems to make detonation IC calculation go out of bounds for Tm
    const auto eq1 = vp * vm * (em - ep) - (pm - pp);
    const auto eq2 = vm * (em + pp) - vp * (ep + pm);
    // const auto eq1 = vp * vm - (pm - pp) / (em - ep);
    // const auto eq2 = vm / vp - (ep + pm) / (em + pp);

    return {eq1, eq2};
}

// keep as void or output just xi_sh and y0={v1UF, w1wN, T1TN}?
// better error handling - change catch (...) to specific throws
void FluidProfile::get_IC_deflagration_veff(const deriv_func& dydxi, double& xi_sh, state_type& y1, state_type& yp, state_type& ym) const {
    /************************ CLOCK / PROFILER *************************/
    const auto ti = std::chrono::high_resolution_clock::now();
    /******************************************************************/

    const double xi_sh_min = std::max(std::sqrt(cpsq_), vw_); // xi_sh > cp > xi_w (deflag), xi_sh > xi_w > cp (hybrid)
    const double xi_sh_max = 1.0;

    // residual functions to find xi_sh satisfying both matching conditions
    std::function<double(double)> vm_resi;
    if (mode_ == 0) { // deflagration: residual = vm - vw
        vm_resi = [this] (double vmUF) { return abs(mu(vw_, abs(vmUF)) - vw_); };
    } else { // hybrid: residual = vm - cm
        vm_resi = [this] (double vmUF) { return abs(mu(vw_, abs(vmUF)) - cm_); };
    }

    auto shock_resi = [this, &vm_resi, &dydxi] (double xi_sh) { 
        auto y1 = test_shock_matching(dydxi, xi_sh);
        auto ypym = test_wall_matching(dydxi, xi_sh, y1);
        return vm_resi(ypym.second[0]);
    };

    // finds interval containing minima of residual vector
    auto check_minimum = [&](const std::vector<double>& xs,
                         const std::vector<double>& res) -> std::pair<double,double>
    {
        if (xs.size() < 3)
            return { std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN() };

        size_t n = xs.size();
        if (res[n-3] > res[n-2] && res[n-1] > res[n-2]) {
            return { xs[n-3], xs[n-1] };  // bracket around minimum
        }

        return { std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN() };
    };

    double xi_sh_guess = xi_sh_min;
    const auto dx = 5e-5; // make this smaller if  can't find bracket for shock

    state_type y1_guess;
    std::pair<state_type, state_type> ypym_guess;
    std::vector<double> pass_case_xi, pass_case_resi;
    std::pair<double, double> bracket;

    // find interval containing xi_sh
    // NOTE: Only works if solution for xi_sh is unique!! (pretty sure this is true)  
    while (xi_sh_guess < xi_sh_max) {
        try {
            y1_guess = test_shock_matching(dydxi, xi_sh_guess);
            ypym_guess = test_wall_matching(dydxi, xi_sh_guess, y1_guess);

            pass_case_xi.push_back(xi_sh_guess);
            pass_case_resi.push_back(vm_resi(ypym_guess.second[0]));

            // break when interval containing residual minima is found
            bracket = check_minimum(pass_case_xi, pass_case_resi);
            if (!std::isnan(bracket.first)) break;

        } catch (...) {} // do nothing if matching fails

        xi_sh_guess += dx;
    }

    if (pass_case_xi.size() == 0 || std::isnan(bracket.first)) {
        throw std::runtime_error("No bracketed interval found for shock!");
    }

    // search for xi_sh that minimises residual
    xi_sh = golden_section_minimize(shock_resi, bracket.first, bracket.second);
    // xi_sh = brent_minimize(shock_resi, bracket.first, bracket.second);

    y1 = test_shock_matching(dydxi, xi_sh); // {v1UF, w1wN, T1TN}
    const auto ypym = test_wall_matching(dydxi, xi_sh, y1);
    yp = ypym.first; // {vpUF, wpwN, TpTN}
    ym = ypym.second; // {vmUF, wmwN, TmTN}

    std::cout << "Found shock front!\n";

    return;
}

// TO DO: update error handling
state_type FluidProfile::get_IC_detonation_veff() const {
    // uses matching conditions to get vm, Tm/TN from vp, Tp/TN
    //      vp * vm = (pm - pp) / (em - ep)
    //      vp / vm = (em + pp) / (ep + pm)

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
    const auto vmUF = mu(vw_, abs(vm));
    const auto TmTN = sol[1];
    const auto wmwN = wm_from_matching(wpwN, vp, vm);

    // do better error handling later
    if (vm >= vw_) {
        throw std::invalid_argument("vm must be < vw for detonation");
    }
    if (TmTN <= TpTN) {
        throw std::invalid_argument("Tm/TN must be > Tp/TN for detonation");
    }
    if (wmwN <= wpwN) {
        throw std::invalid_argument("wm/wN must be > wp/wN for detonation");
    }

    return {vmUF, wmwN, TmTN};
}

/*****************************************************************************************/
/************************************* Both Bag/Veff *************************************/
// from matching condition wp * vp * gamma_p^2 = wm * vm * gamma_m^2
double FluidProfile::wm_from_matching(double wp, double vp, double vm) const {
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
    // auto dydv = [this] (double xi, const state_type& y) -> state_type {
    //     return dydv_vec(xi, y, vw_, cpsq_, cmsq_);
    // };

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

        xi0 = xi_sh - dlt;
        xif = vw_ + dlt;

        // solver
        const auto sol = rk4_solver(dydxi, xi0, xif, y0, n);
        xi_sol_tmp = sol.first;
        y_sol_tmp = sol.second;

        // fill v, w vectors
        for (size_t i = 0; i < xi_sol_tmp.size(); i++) {
            v_sol_tmp.push_back(y_sol_tmp[i][0]);
            w_sol_tmp.push_back(y_sol_tmp[i][1]);
            T_sol_tmp.push_back(y_sol_tmp[i][2]);

            la_sol_tmp.push_back(lambda_s(w_sol_tmp[i]));
        }

        const auto vpUF = v_sol_tmp.back();
        // if (vpUF >= vw_) throw std::invalid_argument("vpUF must be < vw");

        const auto wpwN = w_sol_tmp.back();
        
        // check alp okay
        // Note: need alp for this so must do root-finding BEFORE determining if alp is good/bad
        // BAG MODEL ONLY
        const auto alp = get_alp_wall(vpUF, vw_);

        // std::cout << "alp=" << alp << ", alp_min=" << alp_min_ << ", alp_max=" << alp_max_ << ", alN=" << alN_ << "\n";
        // std::cout << std::setprecision(10) << "xi_sh=" << xi_sh << "\n";
        
        if (alp >= alN_) throw std::invalid_argument("alpha_+ must be < alpha_N");
        if (alp < alp_min_) throw std::invalid_argument("alpha_+ too small for shock");
        if (alp > alp_max_) throw std::invalid_argument("alpha_+ too large for shock");

        if (mode_ == 0) {
            // fix end-value for enthalpy
            const auto vm = -vw_;
            const auto vp = vp_from_matching(vm, alp);
            const auto wpwN = w_sol_tmp.back(); // w(xi_w + dlt) = w+/wN

            w_end_val = wm_from_matching(wpwN, vp, vm); // wm/wN, from matching condition at wall
            // T_end_val = ToTN(w_end_val, cpsq_, cmsq_); // Tm/TN
            T_end_val = get_TmTN(w_end_val);
            la_end_val = lambda_b(w_end_val);

        } else { // hybrid
            // initial conditions for rarefaction wave
            const auto xi0_rf = vw_ - dlt;
            const auto xif_rf = std::sqrt(cmsq_) + dlt;

            state_type y0_rf;

            const auto vm = -std::sqrt(cmsq_);
            // const auto vmUF = mu(vw_, abs(vm));
            const auto vmUF = mu(vw_, abs(vm));
            y0_rf[0] = vmUF * 0.95; // need this to offset singularity in dv/dxi

            const auto vp = mu(vw_, abs(vpUF));
            
            const auto wmwN = wm_from_matching(wpwN, vp, vm);
            y0_rf[1] = wmwN;

            // const auto TmTN = ToTN(wmwN, cpsq_, cmsq_);
            const auto TmTN = get_TmTN(wmwN);
            y0_rf[2] = TmTN;

            const auto [xi_sol_rf_tmp, y_sol_rf_tmp] = rk4_solver(dydxi, xi0_rf, xif_rf, y0_rf, n);

            // combine rarefaction wave with shockwave part of solution
            for (size_t i = 0; i < xi_sol_rf_tmp.size(); i++) {
                xi_sol_tmp.push_back(xi_sol_rf_tmp[i]);
                v_sol_tmp.push_back(y_sol_rf_tmp[i][0]);
                w_sol_tmp.push_back(y_sol_rf_tmp[i][1]);
                T_sol_tmp.push_back(y_sol_rf_tmp[i][2]);

                la_sol_tmp.push_back(lambda_b(y_sol_rf_tmp[i][1]));
            }

            // uncomment this if det part of fluid profile looks wrong/discontinuous
            // const auto [v_sol_rf_tmp, y_sol_rf_tmp] = rk4_solver(dydv, vmUF, 1e-8, {xi0_rf, wmwN, TmTN}, n);

            // // combine rarefaction wave with shockwave part of solution
            // for (size_t i = 0; i < v_sol_rf_tmp.size(); i++) {
            //     v_sol_tmp.push_back(v_sol_rf_tmp[i]);
            //     xi_sol_tmp.push_back(y_sol_rf_tmp[i][0]);
            //     w_sol_tmp.push_back(y_sol_rf_tmp[i][1]);
            //     T_sol_tmp.push_back(y_sol_rf_tmp[i][2]);

            //     la_sol_tmp.push_back(lambda_b(y_sol_rf_tmp[i][1]));
            // }

            xif = xif_rf; // update xif value to behind rarefaction wave
            w_end_val = w_sol_tmp.back();
            T_end_val = T_sol_tmp.back();
            la_end_val = la_sol_tmp.back();
        }

    } else { // detonation
        // cm < xi < xi_w
        xi0 = vw_ - dlt;
        xif = std::sqrt(cmsq_) + dlt;

        y0 = get_IC_detonation();

        // solver
        const auto sol = rk4_solver(dydxi, xi0, xif, y0, n);
        xi_sol_tmp = sol.first;
        y_sol_tmp = sol.second;

        // fill v, w vectors
        // Optimisation Note: pre-allocate memory for these vectors to speed up?
        for (size_t i = 0; i < xi_sol_tmp.size(); i++) {
            v_sol_tmp.push_back(y_sol_tmp[i][0]);
            w_sol_tmp.push_back(y_sol_tmp[i][1]);
            T_sol_tmp.push_back(y_sol_tmp[i][2]);

            la_sol_tmp.push_back(lambda_b(w_sol_tmp[i]));
        }

        // uncomment this if fluid profile looks wrong/discontinuous
        // const auto sol = rk4_solver(dydv, vmUF, 0.0, {vw_, wmwN, TmTN}, n);
        // v_sol_tmp = sol.first;
        // y_sol_tmp = sol.second;

        // // fill v, w vectors
        // // Optimisation Note: pre-allocate memory for these vectors to speed up?
        // for (size_t i = 0; i < v_sol_tmp.size(); i++) {
        //     xi_sol_tmp.push_back(y_sol_tmp[i][0]);
        //     w_sol_tmp.push_back(y_sol_tmp[i][1]);
        //     T_sol_tmp.push_back(y_sol_tmp[i][2]);

        //     la_sol_tmp.push_back(lambda_b(w_sol_tmp[i]));
        // }

        w_end_val = w_sol_tmp.back(); // not sure why
        T_end_val = T_sol_tmp.back();
        la_end_val = la_sol_tmp.back();
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
    // auto dydv = [this] (double xi, const state_type& y) -> state_type {
    //     return dydv_vec(xi, y, vw_, cpsq_, cmsq_);
    // };

    const auto eN = params_.eN();
    const auto wN_inv = 1.0 / params_.wN();

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

        // yp unused currently
        double xi_sh;
        state_type y1, yp, ym; // y = {vUF, wwN, TTN}

        // find shock front and ICs
        get_IC_deflagration_veff(dydxi, xi_sh, y1, yp, ym);
        
        xi0 = xi_sh - dlt;
        xif = vw_ + dlt;

        // a little redundant since get_IC_deflagration_veff() has to calc this anyway (need to calc these again to store in *_sol_tmp vectors though)
        const auto sol = rk4_solver(dydxi, xi0, xif, y1, n);
        xi_sol_tmp = sol.first;
        y_sol_tmp = sol.second;

        // fill v, w vectors
        for (size_t i = 0; i < xi_sol_tmp.size(); i++) {
            v_sol_tmp.push_back(y_sol_tmp[i][0]);
            w_sol_tmp.push_back(y_sol_tmp[i][1]);
            // T_sol_tmp.push_back(y_sol_tmp[i][2]);
            
            const auto T_sol = y_sol_tmp[i][2];
            T_sol_tmp.push_back(T_sol);
            la_sol_tmp.push_back(lambda_s_veff(T_sol, eN, wN_inv)); // lambda in symmetric phase
        }

        const auto vpUF = v_sol_tmp.back();
        const auto wpwN = w_sol_tmp.back();
        const auto TpTN = T_sol_tmp.back();

        if (mode_ == 0) {
            // fix end values 
            w_end_val = ym[1]; // wmwN
            T_end_val = ym[2]; // TmTN
            la_end_val = lambda_b_veff(T_end_val, eN, wN_inv); // lambda just behind wall (broken phase)

        } else { // hybrid
            // initial conditions for rarefaction wave
            const auto xi0_rf = vw_ - dlt;
            const auto xif_rf = std::sqrt(cmsq_) + dlt;

            ym[0] = ym[0] * 0.95; // need this to offset singularity in dv/dxi
            const auto [xi_sol_rf_tmp, y_sol_rf_tmp] = rk4_solver(dydxi, xi0_rf, xif_rf, ym, n);
            // combine rarefaction wave with shockwave part of solution
            for (size_t i = 0; i < xi_sol_rf_tmp.size(); i++) {
                xi_sol_tmp.push_back(xi_sol_rf_tmp[i]);
                v_sol_tmp.push_back(y_sol_rf_tmp[i][0]);
                w_sol_tmp.push_back(y_sol_rf_tmp[i][1]);

                const auto T_sol_rf = y_sol_rf_tmp[i][2];
                T_sol_tmp.push_back(T_sol_rf);
                la_sol_tmp.push_back(lambda_b_veff(T_sol_rf, eN, wN_inv));
            }

            // uncomment this if det part of fluid profile looks wrong/discontinuous
            // const auto [v_sol_rf_tmp, y_sol_rf_tmp] = rk4_solver(dydv, ym[0], 1e-8, {xi0_rf, ym[1], ym[2]}, n);

            // // combine rarefaction wave with shockwave part of solution
            // for (size_t i = 0; i < v_sol_rf_tmp.size(); i++) {
            //     v_sol_tmp.push_back(v_sol_rf_tmp[i]);
            //     xi_sol_tmp.push_back(y_sol_rf_tmp[i][0]);
            //     w_sol_tmp.push_back(y_sol_rf_tmp[i][1]);
            
            //     const auto T_sol_rf = y_sol_rf_tmp[i][2];
            //     T_sol_tmp.push_back(T_sol_rf);
            //     la_sol_tmp.push_back(lambda_b_veff(T_sol_rf, eN, wN_inv));
            // }

            xif = xif_rf; // update xif value to behind rarefaction wave
            w_end_val = w_sol_tmp.back();
            T_end_val = T_sol_tmp.back();
            la_end_val = la_sol_tmp.back();
        }

    } else { // detonation
        // cm < xi < xi_w
        // v0 = v(xi_w) = vm(UF)
        // w0 = w(xi_w) = wm/wN
        // T0 = T(xi_w) = Tm/TN

        xi0 = vw_ - dlt;
        xif = std::sqrt(cmsq_) + dlt;        

        y0 = get_IC_detonation_veff(); // {vmUF, wmwN, TmTN}

        // solver
        const auto sol = rk4_solver(dydxi, xi0, xif, y0, n);
        xi_sol_tmp = sol.first;
        y_sol_tmp = sol.second;

        // fill v, w vectors
        // Optimisation Note: pre-allocate memory for these vectors to speed up?
        for (size_t i = 0; i < xi_sol_tmp.size(); i++) {
            v_sol_tmp.push_back(y_sol_tmp[i][0]);
            w_sol_tmp.push_back(y_sol_tmp[i][1]);

            const auto T_sol = y_sol_tmp[i][2];
            T_sol_tmp.push_back(T_sol);
            la_sol_tmp.push_back(lambda_b_veff(T_sol, eN, wN_inv));
        }

        // uncomment this if fluid profile looks wrong/discontinuous
        // const auto sol = rk4_solver(dydv, y0[0], 1e-8, {vw_, y0[1], y0[2]}, n);
        // v_sol_tmp = sol.first;
        // y_sol_tmp = sol.second;

        // // fill v, w vectors
        // for (size_t i = 0; i < v_sol_tmp.size(); i++) {
        //     xi_sol_tmp.push_back(y_sol_tmp[i][0]);
        //     w_sol_tmp.push_back(y_sol_tmp[i][1]);

        //     const auto T_sol = y_sol_tmp[i][2];
        //     T_sol_tmp.push_back(T_sol);
        //     la_sol_tmp.push_back(lambda_b_veff(T_sol, eN, wN_inv));
        // }

        w_end_val = w_sol_tmp.back(); // not sure why
        T_end_val = T_sol_tmp.back();
        la_end_val = la_sol_tmp.back();
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

// void FluidProfile::format_profiles(prof_type xi_sol, prof_type v_sol, prof_type w_sol, prof_type T_sol, prof_type la_sol) {
//     // xi0, xif, n, w_end_val,, T_end_val,, la_end_val need to be in scope

//     // define start & end points where profile=const (outside integration)
//     const auto xi_start = linspace(0.99, xi0, n); // backwards integration
//     const auto xi_end = linspace(xif, 0.01, n);

//     const prof_type v_start(n, 0.0);
//     const prof_type v_end = v_start;

//     const prof_type w_start(n, 1.0); // wpwN=1 (det), w2wN=1 (def/hyb)
//     const prof_type w_end(n, w_end_val);

//     const prof_type T_start(n, 1.0); // TpTN=1 (det), T2TN=1 (def/hyb)
//     const prof_type T_end(n, T_end_val);

//     const prof_type la_start(n, 0.0); // la(xi>vw)=(eN-eN)/wN=0 (det), la(xi>xi_sh)=(eN-eN)/wN=0 def/hyb)
//     const prof_type la_end(n, la_end_val);

//     // concatenate xi vals
//     xi_sol.insert(xi_sol.end(), xi_start.begin(), xi_start.end());
//     xi_sol.insert(xi_sol.end(), xi_sol_tmp.begin(), xi_sol_tmp.end());
//     xi_sol.insert(xi_sol.end(), xi_end.begin(), xi_end.end());

//     // concatenate v(xi) vals
//     v_sol.insert(v_sol.end(), v_start.begin(), v_start.end());
//     v_sol.insert(v_sol.end(), v_sol_tmp.begin(), v_sol_tmp.end());
//     v_sol.insert(v_sol.end(), v_end.begin(), v_end.end());

//     // concatenate w(xi) vals
//     w_sol.insert(w_sol.end(), w_start.begin(), w_start.end());
//     w_sol.insert(w_sol.end(), w_sol_tmp.begin(), w_sol_tmp.end());
//     w_sol.insert(w_sol.end(), w_end.begin(), w_end.end());

//     // concatenate T(xi) vals
//     T_sol.insert(T_sol.end(), T_start.begin(), T_start.end());
//     T_sol.insert(T_sol.end(), T_sol_tmp.begin(), T_sol_tmp.end());
//     T_sol.insert(T_sol.end(), T_end.begin(), T_end.end());

//     // concatenate la(xi) vals
//     la_sol.insert(la_sol.end(), la_start.begin(), la_start.end());
//     la_sol.insert(la_sol.end(), la_sol_tmp.begin(), la_sol_tmp.end());
//     la_sol.insert(la_sol.end(), la_end.begin(), la_end.end());    

//     // reformat from backwards integration
//     std::reverse(xi_sol.begin(), xi_sol.end());
//     std::reverse(v_sol.begin(), v_sol.end());
//     std::reverse(w_sol.begin(), w_sol.end());
//     std::reverse(T_sol.begin(), T_sol.end());
//     std::reverse(la_sol.begin(), la_sol.end());

//     return;
// }

// template <typename EOS>
// std::vector<prof_type> FluidProfile::solve_profile_new(const EOS& eos, int n) {
//     // check valid hydrodynamic mode
//     if (!(mode_ == 0 || mode_ == 1 || mode_ == 2)) {
//             throw std::invalid_argument("Hydrodynamic mode must be: 0 (deflagration), 1 (hybrid) or 2 (detonation)");
//     }

//     std::cout << "Solving fluid profile for hydrodynamic mode=";
//     if (mode_ == 0) {
//         std::cout << "deflagration";
//     } else if (mode_ == 1) {
//         std::cout << "hybrid";
//     } else {
//         std::cout << "detonation";
//     }
//     std::cout << "\n";

//     // wrapper for hydrodynamic EoM
//     auto dydxi = [this] (double xi, const state_type& y) -> state_type {
//         return dydxi_vec(xi, y, vw_, cpsq_, cmsq_);
//     };

//     // uncomment this if det/hyb fluid profile looks wrong/discontinuous
//     // auto dydv = [this] (double xi, const state_type& y) -> state_type {
//     //     return dydv_vec(xi, y, vw_, cpsq_, cmsq_);
//     // };

//     const auto eN = params_.eN();
//     const auto wN_inv = 1.0 / params_.wN();

//     double xi0, xif;
//     const auto dlt = 0.001; // wall and shocks are discontinuities so start integration just before them
//     std::vector<state_type> y_sol_tmp;
//     state_type y0; 
//     prof_type xi_sol_tmp, v_sol_tmp, w_sol_tmp, T_sol_tmp, la_sol_tmp;
    
//     double w_end_val, T_end_val, la_end_val;

//     if (mode_ < 2) { // deflagration/hybrid
//         if (mode_ = 0) { // deflagration
//             const auto ics = eos.get_IC_deflagration(dydxi);


//         } else { // hybrid

//         }
//     } else { // detonation
//          const auto ics = eos.get_IC_detonation();
//     }

//     prof_type xi_sol, v_sol, w_sol, T_sol, la_sol;
//     format_profiles(xi_sol, v_sol, w_sol, T_sol, la_sol);
    
//     return {xi_sol, v_sol, w_sol, T_sol, la_sol};
// }


} // namespace Hydrodynamics