// profile.cpp
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <cassert>

#include "ap.h"
#include "interpolation.h"

#include "profile.hpp"
#include "phasetransition.hpp"
#include "hydrodynamics.hpp"
#include "physics.hpp"
#include "maths_ops.hpp"

#ifdef ENABLE_MATPLOTLIB
#include "matplotlibcpp.h"
#endif

/*
TO DO:
- update generate_streamplot_data() to include points of interest (fixed pts, detonation/deflag/hybrid regions)
- fix calc of w profile in generate_streamplot_data()
- modify plot() so you can individually plot v,w,la too (i.e. separate plot functions)
- write summary of detonation/deflag/hybrid at top of ctor
- change from using y0 to using FluidState class!
- update xi_start and xi_end so they have same spacing as xi_vals -> makes it easier for integrating v(xi) i think!
*/

/*
Veff stuff
- update w1wN
- update xi_shock
- update v1UF_from_shock
- update la(xi)
- check shooting method is okay
- change speed of sound to use dp/dT and de/dT (use interpolator?)

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
    const auto mu_val = mu(xi, v);
    const auto denom = gammaSq(v) * (1.0 - v * xi) * (mu_val * mu_val / csq - 1.0);
    return (2.0 * v / xi) / denom;
}

double dwdxi(double xi, double v, double w, const double csq) {
    return w * gammaSq(v) * mu(xi, v) * (1.0 + 1.0 / csq) * dvdxi(xi, v, csq);
}

double dTdxi(double xi, double v, double T, const double csq) {
    return T * gammaSq(v) * mu(xi, v) * dvdxi(xi, v, csq);
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

// ctor for passing Veff
// make sure it passes in T/TN, P(T)/wN, e(T)/wN!!
FluidProfile::FluidProfile(const PhaseTransition::PTParams& params, const state_type& veff_T_vals, const state_type& veff_p_vals, const state_type& veff_e_vals, const size_t n)
    : eos_(), // does this need to be a member var?
      params_(params),
      cpsq_(params.cpsq()), cmsq_(params.cmsq()),
      vw_(params.vw()), alN_(params.alN()),
      alp_min_(std::numeric_limits<double>::quiet_NaN()), 
      alp_max_(std::numeric_limits<double>::quiet_NaN()),
      mode_(),
      xi0_(), xif_(),
      y0_(),
      veff_T_vals_(veff_T_vals), veff_p_vals_(veff_p_vals), veff_e_vals_(veff_e_vals), veff_w_vals_(),
      ToTN_interp_(),
      xi_vals_(), v_vals_(), w_vals_(), T_vals_(), la_vals_()
    {
        if (alN_ <= 0.0) throw std::invalid_argument("alN must be > 0");

        // define hydrodynamic mode
        mode_ = get_mode(vw_, cmsq_, alN_);

        // check alN large enough for shock (deflag/hybrid only)
        if (mode_ == 0 || mode_ == 1) {
            const auto alp_minmax = get_alp_minmax(vw_);
            alp_min_ = alp_minmax[0];
            alp_max_ = alp_minmax[1];

            // alN > alp > alp_min (can't properly constrain from above since we need alp)
            if (alN_ <= alp_min_) throw std::invalid_argument("alN too small for shock!");
        }

        if (veff_T_vals_.empty() || veff_p_vals_.empty() || veff_e_vals_.empty()) { // bag eos
            std::cout << "Calculating fluid profile using Bag equation of state\n";
            eos_ = "bag";
        } else { // veff eos
            std::cout << "Calculating fluid profile using generic equation of state from Veff\n";
            eos_ = "veff";

            const auto nT = veff_T_vals_.size();
            if (veff_p_vals_.size() != nT && veff_e_vals_.size() != nT) {
                throw std::runtime_error("Temperature, pressure and energy density vectors must have the same size.");
            }

            // store w(T) vals
            state_type veff_w_vals_(nT);
            for (size_t i = 0; i < nT; ++i) {
                veff_w_vals_[i] = veff_e_vals_[i] + veff_p_vals_[i]; // w=e+p
            }

            // interpolate T/TN(w/wN) for generic EoS
            alglib::real_1d_array ToTN_vals, wowN_vals;
            ToTN_vals.setcontent(nT, veff_T_vals_.data());
            wowN_vals.setcontent(nT, veff_w_vals_.data());

            alglib::spline1dbuildcubic(ToTN_vals, wowN_vals, ToTN_interp_);
        }

        // calculate fluid profiles v(xi), w(xi), la(xi)
        const auto prof = solve_profile(n);

        xi_vals_ = prof[0];
        v_vals_ = prof[1];
        w_vals_ = prof[2];
        T_vals_ = prof[3];
        la_vals_ = prof[4];

        std::cout << "Fluid profile constructed!\n";
    }

// ctor for bag model
FluidProfile::FluidProfile(const PhaseTransition::PTParams& params, const size_t n)
    : FluidProfile(params, {}, {}, {}, n) {}

// Public functions
void FluidProfile::write(const std::string& filename) const {
    std::cout << "Writing fluid profile to disk... ";

    std::ofstream file(filename);
    file << "xi,v,w,la\n";

    for (size_t i = 0; i < xi_vals_.size(); ++i) {
        file << xi_vals_[i] << "," << v_vals_[i] << "," << w_vals_[i] << "," << la_vals_[i] << "," << T_vals_[i] << "\n";
    }
    file.close();

    std::cout << "Fluid profile saved to " << filename << "!\n";

    return;
}

#ifdef ENABLE_MATPLOTLIB
void FluidProfile::plot(const std::string& filename) const {
    namespace plt = matplotlibcpp;

    plt::figure_size(2400, 600);

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
std::vector<state_type> FluidProfile::read(const std::string& filename) const {
    std::cout << "Warning: Read fluid profile does not check PT parameters of input file. Manual entry of PT parameters required!\n";

    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Could not open file " + filename);
    }

    std::string line;
    std::getline(file, line); // Skip header

    state_type xi_vals, v_vals, w_vals, la_vals;

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::array<double, 4> values;
        std::string token;

        for (auto& val : values) {
            if (!std::getline(ss, token, ',')) {
                throw std::runtime_error("Malformed line in " + filename + ": " + line);
            }
            val = std::stod(token);
        }

        xi_vals.push_back(values[0]);
        v_vals.push_back(values[1]);
        w_vals.push_back(values[2]);
        la_vals.push_back(values[3]);
    }

    return {xi_vals, v_vals, w_vals, la_vals};
}

int FluidProfile::get_mode(double vw, double cmsq, double alN) const {
    const auto vwsq = vw * vw;

    if (vwsq < cmsq) return 0; // deflagration
    if (vw < vJ_det(alN)) return 1; // hybrid
    return 2; // detonation
}

// generic for Veff
double FluidProfile::vJ_det(double alp) const {
    return calc_vp(std::sqrt(cmsq_), alp); // vJ(alp) = vp(|vm|=cm, alp)
}

// generic for Veff (not sure when to use which sign though)
double FluidProfile::calc_vp(double vm, double alp) const { // vp(vm,alp) from matching eqs
    const auto sgn = 1.0;
    const auto fac1 = 1.0 / (2.0 * (1.0 + alp));
    const auto fac2 = 1.0 / (3.0 * vm);

    return fac1 * (fac2 + vm + sgn * std::sqrt((fac2 - vm) * (fac2 - vm) + 4.0 * alp * alp + 8.0 * alp / 3.0));
}

// generic for Veff (not sure when to use which sign though)
double FluidProfile::calc_vm(double vp, double alp) const { // inverse of vp(vm,alp)
    // sgn=1 for detonation, not sure what condition specifically fixes this though
    // sgn=-1 for deflag/hybrid?
    const auto vp_abs = abs(vp);
    const auto sgn = 1.0;
    const auto fac = (1.0 + alp) * vp_abs + (1.0 - 3.0 * alp) / (3.0 * vp_abs);

    return 0.5 * (fac + sgn * std::sqrt(fac * fac - 4.0 / 3.0));
}

// generic for Veff
// from matching condition w+*v+*gamma+^2=w-*v-*gamma-^2
double FluidProfile::calc_wm(double wp, double vp, double vm) const {
    return wp * abs(vp) * (1.0 - vm * vm) / (abs(vm) * (1.0 - vp * vp));
}

// bag model only
double FluidProfile::calc_w1wN(double xi_sh) const { // w1/wN
    // alpha_1 w1 = alpha_N wN
    const auto xi_sh_sq = xi_sh * xi_sh;
    return (9.0 * xi_sh_sq - 1.0) / (3.0 * (1.0 - xi_sh_sq));
}

// generic for Veff
// passes in cpsq, cmsq since it is reused to get T1/TN (at shock front) with c1sq=c2sq=cpsq
// not sure if (1+cpsq)/(1+cmsq) is correct
double FluidProfile::calc_ToTN(double wmwN, double cpsq, double cmsq) const {
    if (eos_ == "bag") {
        // w=(4/3)*a*T^4 -> w/wN = (a/aN) * (T/TN)^4
        const auto fac = wmwN * (1.0 + cpsq) / (1.0 + cmsq);
        return std::pow(fac, 0.25);
    }
    return alglib::spline1dcalc(ToTN_interp_, wmwN); // call interpolator for generic EoS
}

// bag model only
double FluidProfile::xi_shock(double v1UF) const {
    
    // const auto fac0 = 0.5 * (1.0 - cpsq_) * v1UF;
    // return fac0 + std::sqrt(fac0 * fac0 + cpsq_);

    const auto fac = v1UF / 3.0;
    const auto xi_sh = fac + std::sqrt(fac * fac + 1.0 / 3.0);
    if (xi_sh <= std::sqrt(cpsq_)) throw std::runtime_error("shock failed (xi_sh <= c+)");

    return xi_sh;

    // if (model_ == "bag") { // add one for mu nu model too?
    // } else { // generic EoS
    //     throw std::runtime_error("Veff support not fully implemented yet");

    //     // pressure/energy density at shock
    //     // p1 = p(T(xi_sh))
    //     const auto p1 = 0.0; // behind shock
    //     const auto e1 = 0.0;

    //     const auto p2 = 0.0; // in front of shock
    //     const auto e2 = 0.0;

    //     const auto v2v1 = (p2 - p1) / (e2 - e1); // v2*v1 from matching cond.
    //     const auto v2v1_rat = (e1 + p2) / (e2 + p1); // v2/v1
    // }
}

// bag model only
double FluidProfile::v1UF_from_shock(double xi_sh) const {
    if (xi_sh < std::sqrt(cpsq_) || xi_sh > 1.0) {
        // condition relaxed in if statement since for some vw & alN, xi_shock VERY close to bounds
        // so root-finder takes xi_sh_min = cp, xi_sh_max = 1
        throw std::invalid_argument("shock must be supersonic and less than speed of light (cp < xi_sh < 1)");
    }

    if (xi_sh == std::sqrt(cpsq_)) return 1e-10; // avoid numerical precision errors
    return (3.0 * xi_sh * xi_sh - 1.0) / (2.0 * xi_sh);
}

// generic for Veff
// might need to fix al_min = 0 if numerical precision causes it to be slightly negative
std::vector<double> FluidProfile::get_alp_minmax(double vw) const {
    const auto cp = std::sqrt(cpsq_);

    const auto vm = std::min(cp, vw); // vw for deflag (vw < cm), cm for hybrid (cp < vw)
    const auto vp_min = 0.0;
    const auto vp_max = vm; // |v+| < |v-|
    
    // same as get_alp_wall but using vp, vm
    auto get_alp = [] (double vp, double vm) {
        return gammaSq(vp) * (vp * vp - vp * vm - vp / (3.0 * vm) + 1.0 / 3.0);
    };
    
    const auto al_max = get_alp(vp_min, vm);
    const auto al_min = get_alp(vp_max, vm);

    return {al_min, al_max};
}

// generic for Veff
// alpha_+ from wall condition
double FluidProfile::get_alp_wall(double vpUF, double vw) const {
    return gammaSq(vpUF) * vpUF * (2.0 * vw * vpUF + 1.0 - 3.0 * vw * vw) / (3.0 * vw);
}

// generic for Veff
double FluidProfile::alN_residual_func(double xi_sh, const deriv_func& dydxi) const {
    // initial conditions
    const auto xi0 = xi_sh - 0.001;
    const auto v1UF = v1UF_from_shock(xi_sh);
    const auto w1wN = calc_w1wN(xi_sh);

    const std::vector<double> y0 = {v1UF, w1wN}; // v0 = v(xi_sh) = v1UF

    // solve fluid EoM to get vpUF
    // WARNING: choosing num steps too small gives bad result!
    const auto [xi_sol, y_sol] = rk4_solver(dydxi, xi0, xif_, y0, 10000);
    const auto vpUF = y_sol.back()[0]; // vpUF = v(xi_w) (endpoint of integration)
    const auto wpwN = y_sol.back()[1];
    
    // calc alpha_N from wall constraint
    const auto alN_wall = wpwN * get_alp_wall(vpUF, vw_);

    // std::cout << "log(alN_rat)=" << std::log(std::abs(alN_wall / alN_)) << ", v1UF=" << v1UF << ", xi_sh=" << xi_sh << ", vpUF = " << vpUF << "\n";
    
    // return alN_wall - alN_;
    return std::log(std::abs(alN_wall / alN_)); // doesn't always work for some vw, alN
}

// bag model only
double FluidProfile::get_la_behind_wall(double w) const {
    // la(xi)=(3/4)*(w(xi)/wN - 1 - alN) behind bubble wall (detonations)
    return 0.75 * (w - 1.0 - alN_);
}

// bag model only
double FluidProfile::get_la_front_wall(double w) const {
    // la(xi)=(3/4)*(w(xi)/wN - 1) in front of bubble wall (deflagrations)
    return 0.75 * (w - 1.0);
}

// generic for Veff
double FluidProfile::find_shock(const deriv_func& dydxi) const {
    // Root-finding algorithm for initial condition v0 = v(xi_sh) = v1UF

    // residual function f(xi_sh) = alN_calc - alN_actual
    std::function<double(double)> residual = [this, &dydxi] (double xi_sh) {
        return alN_residual_func(xi_sh, dydxi);
    };

    // cp < xi_sh < 1 (shock must be supersonic and less than speed of light)
    const double xi_sh_min = std::max(std::sqrt(cpsq_), vw_); // xi_sh > cp > xi_w (deflag), xi_sh > xi_w > cp (hybrid)
    const double xi_sh_max = 1.0 - 1e-5; // need to make this closer to 1 for extreme case of hybrids with xi_sh very close to 1

    // const auto xi_sh_vals = linspace(xi_sh_min, xi_sh_max, 500);
    // const auto xi_sh_vals = linspace(xi_sh_min, xi_sh_min + 0.03, 500);
    // std::vector<double> res_vals;
    // for (int i = 0; i < xi_sh_vals.size(); i++) {
    //     const auto val = residual(xi_sh_vals[i]);
    //     res_vals.push_back(val);
    //     // if (val < 0.0) {
    //     //     std::cout << "res<0 for xi_sh[" << i << "]=" << xi_sh_vals[i] << "\n";
    //     // }
    // }

    // plt::figure_size(800, 600);
    // plt::plot(xi_sh_vals, res_vals, "k-");
    // // plt::xlim(xi_sh_min - 0.01, xi_sh_max + 0.01);
    // plt::xlim(xi_sh_min - 0.01, xi_sh_min + 0.04);
    // // plt::ylim(-1.0, 1.0);
    // plt::grid(true);
    // plt::save("../residual.png");

    return root_finder(residual, xi_sh_min, xi_sh_max, 1e-7, 100);
    // return find_smallest_root(residual, xi_sh_min, xi_sh_max); // very slow
}

std::vector<state_type> FluidProfile::solve_profile(int n) {
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

    // wrapper for hydrodynamic EoM - move into struct?
    auto dydxi = [this](double xi, const state_type& y) -> state_type {
        const auto v = y[0];
        const auto w = y[1];
        const auto csq = (xi < vw_) ? cmsq_ : cpsq_;

        return { dvdxi(xi, v, csq), dwdxi(xi, v, w, csq), dTdxi(xi, v, w, csq) };
    };

    const auto dlt = 0.001; // wall and shocks are discontinuities so start integration just before them
    std::vector<state_type> y_sol_tmp;
    state_type xi_sol_tmp, v_sol_tmp, w_sol_tmp, T_sol_tmp, la_sol_tmp;

    const auto w_start_val = 1.0; // w+/wN (det), w2/wN (deflag/hybrid)
    const auto T_start_val = 1.0; // T+/TN (det), T2/TN (deflag/hybrid);
    const auto la_start_val = 0.0;
    
    double w_end_val, T_end_val, la_end_val;

    if (mode_ < 2) { // deflagration & hybrid
        // hybrid and deflagration ICs the same for xi_w < xi < xi_sh
        xif_ = vw_ + dlt;

        const auto xi_sh = find_shock(dydxi);
        const auto v1UF = v1UF_from_shock(xi_sh);
        y0_.push_back(v1UF);

        xi0_ = xi_shock(v1UF) - dlt;

        // initial condition w(xi_sh) = w1/wN
        const auto w1wN = calc_w1wN(xi0_);
        y0_.push_back(w1wN);

        const auto T1TN = calc_ToTN(w1wN, cpsq_, cpsq_); // c1sq = c2sq = cpsq
        y0_.push_back(T1TN);

        // solver
        const auto sol = rk4_solver(dydxi, xi0_, xif_, y0_, n);
        xi_sol_tmp = sol.first;
        y_sol_tmp = sol.second;

        // fill v, w vectors
        for (size_t i = 0; i < xi_sol_tmp.size(); i++) {
            v_sol_tmp.push_back(y_sol_tmp[i][0]);
            w_sol_tmp.push_back(y_sol_tmp[i][1]);
            T_sol_tmp.push_back(y_sol_tmp[i][2]);

            la_sol_tmp.push_back(get_la_front_wall(w_sol_tmp[i]));
        }

        const auto vpUF = v_sol_tmp.back();
        // if (vpUF >= vw_) throw std::invalid_argument("vpUF must be < vw");

        const auto wpwN = w_sol_tmp.back();
        
        // check alp okay
        // Note: need alp for this so must do root-finding BEFORE determining if alp is good/bad
        const auto alp = get_alp_wall(vpUF, vw_);

        // std::cout << "alp=" << alp << ", alp_min=" << alp_min_ << ", alp_max=" << alp_max_ << ", alN=" << alN_ << "\n";
        // std::cout << std::setprecision(10) << "xi_sh=" << xi_sh << "\n";
        
        if (alp >= alN_) throw std::invalid_argument("alpha_+ must be < alpha_N");
        if (alp < alp_min_) throw std::invalid_argument("alpha_+ too small for shock");
        if (alp > alp_max_) throw std::invalid_argument("alpha_+ too large for shock");

        if (mode_ == 0) {
            // fix end-value for enthalpy
            const auto vm = -vw_;
            const auto vp = calc_vp(vm, alp);
            const auto wpwN = w_sol_tmp.back(); // w(xi_w + dlt) = w+/wN

            w_end_val = calc_wm(wpwN, vp, vm); // wm/wN, from matching condition at wall
            T_end_val = calc_ToTN(w_end_val, cpsq_, cmsq_); // Tm/TN
            la_end_val = get_la_behind_wall(w_end_val);

        } else { // hybrid
            // initial conditions for rarefaction wave
            // const auto xi0_rf = vw_ - dlt;
            const auto xi0_rf = vw_ + dlt;
            const auto xif_rf = std::sqrt(cmsq_) + dlt;

            std::vector<double> y0_rf;

            const auto vm = -std::sqrt(cmsq_);
            const auto vmUF = mu(vw_, abs(vm));
            y0_rf.push_back(vmUF);

            const auto vp = mu(vw_, abs(vpUF));
            
            const auto wmwN = calc_wm(wpwN, vp, vm);
            y0_rf.push_back(wmwN);

            const auto TmTN = calc_ToTN(wmwN, cpsq_, cmsq_);
            y0_rf.push_back(TmTN);

            const auto [xi_sol_rf_tmp, y_sol_rf_tmp] = rk4_solver(dydxi, xi0_rf, xif_rf, y0_rf, n);

            // combine rarefaction wave with shockwave part of solution
            for (size_t i = 0; i < xi_sol_rf_tmp.size(); i++) {
                xi_sol_tmp.push_back(xi_sol_rf_tmp[i]);
                v_sol_tmp.push_back(y_sol_rf_tmp[i][0]);
                w_sol_tmp.push_back(y_sol_rf_tmp[i][1]);
                T_sol_tmp.push_back(y_sol_rf_tmp[i][2]);

                la_sol_tmp.push_back(get_la_behind_wall(y_sol_rf_tmp[i][1]));
            }

            xif_ = xif_rf; // update xif value to behind rarefaction wave
            w_end_val = w_sol_tmp.back();
            T_end_val = T_sol_tmp.back();
            la_end_val = la_sol_tmp.back();
        }

        // std::cout << "v_start=" << v_sol_tmp.front() << ", v_end=" << v_sol_tmp.back() << "\n";
        // std::cout << "w_start=" << w_sol_tmp.front() << ", w_end=" << w_sol_tmp.back() << "\n";
        // std::cout << "T_start=" << T_sol_tmp.front() << ", T_end=" << T_sol_tmp.back() << "\n";
        // std::cout << "la_start=" << la_sol_tmp.front() << ", la_end=" << la_sol_tmp.back() << "\n";

    } else { // detonation
        // cm < xi < xi_w
        xi0_ = vw_ - dlt;
        xif_ = std::sqrt(cmsq_) + dlt;

        // alpha_+ = alpha_N
        const auto alp = alN_;

        // initial condition v0 = v(xi_w) = vm(UF)
        const auto vp = -vw_;
        const auto vm = calc_vm(vp, alp);
        const auto vmUF = mu(vw_, vm);
        y0_.push_back(vmUF);

        // initial condition w0 = w(xi_w) = wm/wN
        const auto wpwN = 1.0; // w+ = wN
        const auto wmwN = calc_wm(wpwN, vp, vm);
        y0_.push_back(wmwN);

        // initial condition T0 = T(xi_w) = Tm/TN
        // const auto TmTN = calc_ToTN(wmwN, cpsq_, cmsq_);
        const auto TmTN = 1.0;
        y0_.push_back(TmTN);

        // solver
        // [xi_sol_tmp, y_sol_tmp] = rk4_solver(dydxi, xi0_, xif_, y0_, n);
        const auto sol = rk4_solver(dydxi, xi0_, xif_, y0_, n);
        xi_sol_tmp = sol.first;
        y_sol_tmp = sol.second;

        // fill v, w vectors
        // Optimisation Note: pre-allocate memory for these vectors to speed up?
        for (size_t i = 0; i < xi_sol_tmp.size(); i++) {
            v_sol_tmp.push_back(y_sol_tmp[i][0]);
            w_sol_tmp.push_back(y_sol_tmp[i][1]);
            T_sol_tmp.push_back(y_sol_tmp[i][2]);

            la_sol_tmp.push_back(get_la_behind_wall(w_sol_tmp[i]));
        }

        w_end_val = w_sol_tmp.back(); // not sure why
        T_end_val = T_sol_tmp.back();
        la_end_val = la_sol_tmp.back();

        // std::cout << "v_start=" << v_sol_tmp.front() << ", v_end=" << v_sol_tmp.back() << "\n";
        // std::cout << "w_start=" << w_sol_tmp.front() << ", w_end=" << w_sol_tmp.back() << "\n";
        // std::cout << "T_start=" << T_sol_tmp.front() << ", T_end=" << T_sol_tmp.back() << "\n";
        // std::cout << "la_start=" << la_sol_tmp.front() << ", la_end=" << la_sol_tmp.back() << "\n";
    }
    

    // define start & end points where profile=const (outside integration)
    const auto xi_start = linspace(0.99, xi0_, n); // backwards integration
    const auto xi_end = linspace(xif_, 0.01, n);

    const state_type v_start(n, 0.0);
    const state_type v_end = v_start;

    const state_type w_start(n, w_start_val);
    const state_type w_end(n, w_end_val);

    const state_type T_start(n, T_start_val);
    const state_type T_end(n, T_end_val);

    const state_type la_start(n, la_start_val);
    const state_type la_end(n, la_end_val);

    state_type xi_sol, v_sol, w_sol, T_sol, la_sol;

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

} // namespace Hydrodynamics