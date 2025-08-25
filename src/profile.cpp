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
namespace plt = matplotlibcpp;
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
FluidProfile::FluidProfile(const PhaseTransition::PTParams& params, const state_type& veff_T_vals, const state_type& veff_ps_vals, const state_type& veff_pb_vals, const state_type& veff_es_vals, const state_type& veff_eb_vals, const size_t n)
    : eos_(),
      params_(params),
      cpsq_(params.cpsq()), cmsq_(params.cmsq()),
      vw_(params.vw()), alN_(params.alN()),
      alp_min_(std::numeric_limits<double>::quiet_NaN()), 
      alp_max_(std::numeric_limits<double>::quiet_NaN()),
      mode_(),
      xi0_(), xif_(),
      y0_(),
      veff_TTN_vals_(),
      veff_ps_vals_(veff_ps_vals), veff_pb_vals_(veff_pb_vals),
      veff_es_vals_(veff_es_vals), veff_eb_vals_(veff_eb_vals),
      xi_vals_(), v_vals_(), w_vals_(), T_vals_(), la_vals_()
    {
        if (veff_T_vals.empty()) { // bag eos
            std::cout << "Calculating fluid profile using Bag equation of state\n";
            eos_ = "bag";

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

        
        } else { // veff eos
            std::cout << "Calculating fluid profile using generic equation of state from Veff\n";
            std::cout << "Warning: alN stored in PTParams is not used for Veff EoS!\n";
            eos_ = "veff";

            const auto nT = veff_T_vals.size();
            if (veff_ps_vals_.size() != nT && veff_pb_vals_.size() != nT && 
                veff_es_vals_.size() != nT && veff_eb_vals_.size() != nT) {
                throw std::runtime_error("Temperature, pressure and energy density vectors must have the same size.");
            }

            // replace with get_mode_veff when implemented!!
            mode_ = get_mode(vw_, cmsq_, alN_); // placeholder for now!!

            // delete when get_mode_veff is implemented
            if (mode_ == 0 || mode_ == 1) {
                const auto alp_minmax = get_alp_minmax(vw_);
                alp_min_ = alp_minmax[0];
                alp_max_ = alp_minmax[1];

                // alN > alp > alp_min (can't properly constrain from above since we need alp)
                if (alN_ <= alp_min_) throw std::invalid_argument("alN too small for shock!");
            }

            // store T/TN vals
            const auto TN_inv = 1.0 / params_.TN();
            for (size_t i = 0; i < nT; i++) {
                veff_TTN_vals_.push_back(veff_T_vals[i] * TN_inv);
            }

            // store w(T/TN) vals
            state_type veff_ws_vals_(nT), veff_wb_vals_(nT);
            for (size_t i = 0; i < nT; ++i) {
                veff_ws_vals_[i] = veff_es_vals_[i] + veff_ps_vals_[i]; // w=e+p
                veff_wb_vals_[i] = veff_eb_vals_[i] + veff_pb_vals_[i];
            }

            // construct interpolating functions p(T/TN), e(T/TN) in s/b phases
            alglib::real_1d_array veff_TTN_array, veff_ps_array, veff_es_array, veff_ws_array, veff_pb_array, veff_eb_array, veff_wb_array;
            veff_TTN_array.setcontent(nT, veff_TTN_vals_.data());

            veff_ps_array.setcontent(nT, veff_ps_vals_.data()); // symmetric phase
            veff_es_array.setcontent(nT, veff_es_vals_.data());
            veff_ws_array.setcontent(nT, veff_ws_vals_.data());
            alglib::spline1dbuildcubic(veff_TTN_array, veff_ps_array, veff_ps_interp_);
            alglib::spline1dbuildcubic(veff_TTN_array, veff_es_array, veff_es_interp_);
            alglib::spline1dbuildcubic(veff_TTN_array, veff_ws_array, veff_ws_interp_);

            veff_pb_array.setcontent(nT, veff_pb_vals_.data()); // broken phase
            veff_eb_array.setcontent(nT, veff_eb_vals_.data());
            veff_wb_array.setcontent(nT, veff_wb_vals_.data());
            alglib::spline1dbuildcubic(veff_TTN_array, veff_pb_array, veff_pb_interp_);
            alglib::spline1dbuildcubic(veff_TTN_array, veff_eb_array, veff_eb_interp_);
            alglib::spline1dbuildcubic(veff_TTN_array, veff_wb_array, veff_wb_interp_);
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
    : FluidProfile(params, {}, {}, {}, {}, {}, n) {}

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

void FluidProfile::plot_thermo(const std::string& filename) const {
    if (eos_ != "veff")
        throw std::runtime_error("plot_thermo can only be called when using Veff");

    plt::figure_size(2400, 1600);

    plt::subplot2grid(3, 2, 0, 0);
    plt::plot(veff_TTN_vals_, veff_es_vals_);
    plt::xlabel("T/TN");
    plt::ylabel("es(T/TN)");
    plt::grid(true);

    plt::subplot2grid(3, 2, 0, 1);
    plt::plot(veff_TTN_vals_, veff_eb_vals_);
    plt::xlabel("T/TN");
    plt::ylabel("eb(T/TN)");
    plt::grid(true);

    plt::subplot2grid(3, 2, 1, 0);
    plt::plot(veff_TTN_vals_, veff_ps_vals_);
    plt::xlabel("T/TN");
    plt::ylabel("ps(T/TN)");
    plt::grid(true);

    plt::subplot2grid(3, 2, 1, 1);
    plt::plot(veff_TTN_vals_, veff_pb_vals_);
    plt::xlabel("T/TN");
    plt::ylabel("pb(T/TN)");
    plt::grid(true);

    plt::subplot2grid(3, 2, 2, 0);
    plt::plot(veff_TTN_vals_, veff_ws_vals_);
    plt::xlabel("T/TN");
    plt::ylabel("ws(T/TN)");
    plt::grid(true);

    plt::subplot2grid(3, 2, 2, 1);
    plt::plot(veff_TTN_vals_, veff_wb_vals_);
    plt::xlabel("T/TN");
    plt::ylabel("wb(T/TN)");
    plt::grid(true);

    plt::save(filename);

    return;
}
#endif


// calculate csq(T)
// state_type FluidProfile::get_csq() const {
//     const size_t n = veff_T_vals_.size();
//     if (n == 0) {
//         std::cerr << "Cannot call get_csq() for Bag EoS!\n";
//         return {};
//     }

//     state_type csq_vals(n);
//     double p, dpdT, dp2dT2;
//     double e, dedT, de2dT2;

//     for (size_t i = 0; i < n; ++i) {
//         const auto T = veff_T_vals_[i];
//         spline1ddiff(veff_p_interp_, T, p, dpdT, dp2dT2);
//         spline1ddiff(veff_e_interp_, T, e, dedT, de2dT2);

//         csq_vals[i] = dpdT / dedT;
//     }

//     return csq_vals;
// }

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

// not finished - hybrids not done yet
// should only be used in region where EoM is solved (otherwise need to change def of e_xi to check if xi corresponds to s or b phase)
state_type FluidProfile::get_lambda(state_type T_vals) const {
    // la(xi) = (e(xi) - eN) / wN
    const auto eN = alglib::spline1dcalc(veff_es_interp_, 1.0); // eN = e(T/TN=1)
    const auto wN = alglib::spline1dcalc(veff_ws_interp_, 1.0); // wN = w(T/TN=1)
    const auto wN_inv = 1.0 / wN;

    // define helper for computing e(xi)=e(T(xi))
    std::function<double(size_t)> e_xi;

    // add error handling if i > T_vals.size()
    if (mode_ == 0) { // deflagration - all xi in symmetric phase (xi_w < xi < xi_sh)
        e_xi = [&] (size_t i) { return alglib::spline1dcalc(veff_es_interp_, T_vals[i]); };
    } else if (mode_ == 1) { // hybrid
        e_xi = [&] (size_t i) { return 0.0; };
    } else { // detonantion - all xi in broken phase (cm < xi < xi_w)
        e_xi = [&] (size_t i) { return alglib::spline1dcalc(veff_eb_interp_, T_vals[i]); };
    }

    const size_t n = T_vals.size();
    state_type la_vals(n);

    for (size_t i = 0; i < n; i++) {
        la_vals[i] = (e_xi(i) - eN) * wN_inv;
    }

    return la_vals;
}

double FluidProfile::veff_residual_func(double xi_sh, const deriv_func& dydxi) const {
    // fluid in front of shock
    const auto v2 = -xi_sh;
    const auto w2wN = 1.0;
    const auto T2TN = 1.0;

    const auto TTN_min = veff_TTN_vals_.front();
    const auto TTN_max = veff_TTN_vals_.back();

    std::cout << "test1 - veff_residual_func called!\n";

    // matching across shock
    // std::function<std::vector<double>(std::vector<double>)> shock_matching_helper = [this, v2, T2TN] (std::vector<double> y0) {
    //     return matching_eqs_shock(v2, T2TN, y0[0], y0[1]); // y0 = {v1, T1TN}
    // };

    std::function<double(double)> shock_matching_helper2 = [this, v2, T2TN] (double T1TN) {
        return matching_eqs_shock2(v2, T2TN, T1TN);
    };

    // testing stuff
    const auto p2 = alglib::spline1dcalc(veff_ps_interp_, T2TN);
    const auto e2 = alglib::spline1dcalc(veff_es_interp_, T2TN);

    const auto T1TN_vals = linspace(TTN_min, TTN_max, 500);
    std::vector<double> resi_vals;
    for (const auto& T1TN : T1TN_vals) {
        resi_vals.push_back(shock_matching_helper2(T1TN));
    }

    plt::figure_size(800, 600);
    plt::plot(T1TN_vals, resi_vals);
    plt::grid(true);
    plt::save("residual_sh_matching.png");
    // testing stuff end
    
    // const auto shock_sol = newton_solve(shock_matching_helper, {v2, T2TN}); // initial guess v2, T2TN
    // const auto T1TN = root_finder(shock_matching_helper2, TTN_min, TTN_max);
    double T1TN;
    try {
        T1TN = root_finder(shock_matching_helper2, TTN_min, TTN_max);
    } catch (...) {
        T1TN = T2TN + 0.01;
    }
    
    const auto p1 = alglib::spline1dcalc(veff_ps_interp_, T1TN);
    const auto e1 = alglib::spline1dcalc(veff_es_interp_, T1TN);
    // const auto p2 = alglib::spline1dcalc(veff_ps_interp_, T2TN);
    // const auto e2 = alglib::spline1dcalc(veff_es_interp_, T2TN);

    const auto v1 = (p1 - p2) * (e2 + p1) / ((e1 - e2) * (e1 + p2));
    const auto v1UF = mu(xi_sh, v1);
    const auto w1wN = calc_wm(w2wN, v2, v1);

    std::cout << "test2 - shock_sol found for given xi_sh!\n";
    
    // fluid behind shock
    // const auto v1 = shock_sol[0];
    // const auto v1UF = mu(xi_sh, v1);
    // const auto w1wN = calc_wm(w2wN, v2, v1);
    // const auto T1TN = shock_sol[1];

    // solve EoM from shock to wall
    const auto xi0 = xi_sh - 0.001;
    const auto y0 = {v1UF, w1wN, T1TN};
    const auto [xi_sol, y_sol] = rk4_solver(dydxi, xi0, xif_, y0, 10000);

    std::cout << "test3 - EoM solved!\n";

    const auto vpUF = y_sol.back()[0]; // vpUF = v(xi_w) (endpoint of integration)
    const auto vp = mu(vw_, abs(vpUF));
    const auto wpwN = y_sol.back()[1];
    const auto TpTN = y_sol.back()[2];

    std::cout << "vpUF=" << vpUF << ", vp=" << vp << ", wpwN=" << wpwN << ", TpTN=" << TpTN << "\n";

    // matching across wall
    // std::function<std::vector<double>(std::vector<double>)> wall_matching_helper = [this, vp, TpTN] (std::vector<double> y0) {
    //     return matching_eqs_wall(vp, TpTN, y0[0], y0[1]); // y0 = {vm, TmTN}
    // };

    std::function<double(double)> wall_matching_helper2 = [this, vp, TpTN] (double TmTN) {
        return matching_eqs_wall2(vp, TpTN, TmTN);
    };

    // testing stuff
    const auto pp = alglib::spline1dcalc(veff_ps_interp_, TpTN);
    const auto ep = alglib::spline1dcalc(veff_es_interp_, TpTN);

    const auto TmTN_vals = linspace(TTN_min, TTN_max, 500);
    std::vector<double> resi_vals2;
    for (const auto& TmTN : TmTN_vals) {
        resi_vals2.push_back(wall_matching_helper2(TmTN));
    }

    plt::figure_size(800, 600);
    plt::plot(TmTN_vals, resi_vals2);
    plt::grid(true);
    plt::save("residual_wall_matching.png");
    // testing stuff end

    // const auto wall_sol = newton_solve(wall_matching_helper, {vp, TpTN}); // initial guess vp, TpTN
    // const auto TmTN = root_finder(wall_matching_helper2, TTN_min, TTN_max);
    double TmTN;
    try {
        TmTN = root_finder(wall_matching_helper2, TTN_min, TTN_max);
    } catch (...) {
        TmTN = TpTN + 0.01;
    }

    std::cout << "test4 - wall_sol found!\n";

    const auto pm = alglib::spline1dcalc(veff_pb_interp_, TmTN);
    const auto em = alglib::spline1dcalc(veff_eb_interp_, TmTN);
    // const auto pp = alglib::spline1dcalc(veff_ps_interp_, TpTN);
    // const auto ep = alglib::spline1dcalc(veff_es_interp_, TpTN);

    const auto vm = (pm - pp) * (ep + pm) / ((em - ep) * (em + pp));

    // fluid behind wall
    // const auto vm = wall_sol[0];
    // const auto wmwN = calc_wm(wpwN, vp, vm);
    // const auto TmTN = wall_sol[1];

    std::cout << "res=" << vm-vw_ << "\n";

    // residual
    return std::abs(vm) - vw_;
}

double FluidProfile::find_shock_veff(const deriv_func& dydxi) const {
    auto residual = [this, &dydxi] (double xi_sh) {
        return veff_residual_func(xi_sh, dydxi);
    };

    const double xi_sh_min = std::max(std::sqrt(cpsq_), vw_);
    const double xi_sh_max = 1.0;

    // return newton_solve_1d(residual, xi_sh_min);
    return root_finder(residual, xi_sh_min, xi_sh_max, 1e-7, 100);
}

double FluidProfile::find_shock(const deriv_func& dydxi) const {
    // Root-finding algorithm for initial condition v0 = v(xi_sh) = v1UF
    
    // f(xi_sh) = alN_calc - alN_actual
    auto residual = [this, &dydxi] (double xi_sh) {
        return alN_residual_func(xi_sh, dydxi);
    };

    // cp < xi_sh < 1 (shock must be supersonic and less than speed of light)
    const double xi_sh_min = std::max(std::sqrt(cpsq_), vw_); // xi_sh > cp > xi_w (deflag), xi_sh > xi_w > cp (hybrid)
    const double xi_sh_max = 1.0 - 1e-5; // need to make this closer to 1 for extreme case of hybrids with xi_sh very close to 1

    return root_finder(residual, xi_sh_min, xi_sh_max, 1e-7, 100);
}

std::vector<double> FluidProfile::matching_eqs_shock(double v2, double T2TN, double v1, double T1TN) const {    
    if (T1TN < veff_TTN_vals_.front() || T1TN > veff_TTN_vals_.back()) {
        throw std::out_of_range("T1/TN is called out of bounds for spline!");
    }

    const auto p2 = alglib::spline1dcalc(veff_ps_interp_, T2TN); // p_2, e_2
    const auto e2 = alglib::spline1dcalc(veff_es_interp_, T2TN);

    const auto p1 = alglib::spline1dcalc(veff_ps_interp_, T1TN); // p_1, e_1
    const auto e1 = alglib::spline1dcalc(veff_es_interp_, T1TN);

    const auto eq1 = v2 * v1 * (e1 - e2) - (p1 - p2);
    const auto eq2 = v1 * (e1 + p2) - v2 * (e2 + p1);
    // const auto eq1 = v2 * v1 - (p1 - p2) / (e1 - e2);
    // const auto eq2 = v1 / v2 - (e2 + p1) / (e1 + p2);

    return {eq1, eq2};
}

double FluidProfile::matching_eqs_shock2(double v2, double T2TN, double T1TN) const {    
    if (T1TN < veff_TTN_vals_.front() || T1TN > veff_TTN_vals_.back()) {
        throw std::out_of_range("T1/TN is called out of bounds for spline!");
    }

    const auto p2 = alglib::spline1dcalc(veff_ps_interp_, T2TN); // p_2, e_2
    const auto e2 = alglib::spline1dcalc(veff_es_interp_, T2TN);

    const auto p1 = alglib::spline1dcalc(veff_ps_interp_, T1TN); // p_1, e_1
    const auto e1 = alglib::spline1dcalc(veff_es_interp_, T1TN);

    return v2 * v2 - (p1 - p2) * (e1 + p2) / ((e1 - e2) * (e2 + p1));
    // return v2 * v2 * (e1 - e2) * (e2 + p1) - (p1 - p2) * (e1 + p2); 
}

std::vector<double> FluidProfile::matching_eqs_wall(double vp, double TpTN, double vm, double TmTN) const {    
    if (TmTN < veff_TTN_vals_.front() || TmTN > veff_TTN_vals_.back()) {
        throw std::out_of_range("Tm/TN is called out of bounds for spline!");
    }

    const auto pp = alglib::spline1dcalc(veff_ps_interp_, TpTN); // p_+, e_+
    const auto ep = alglib::spline1dcalc(veff_es_interp_, TpTN);

    const auto pm = alglib::spline1dcalc(veff_pb_interp_, TmTN); // p_-, e_-
    const auto em = alglib::spline1dcalc(veff_eb_interp_, TmTN);

    const auto eq1 = vp * vm * (em - ep) - (pm - pp);
    const auto eq2 = vm * (em + pp) - vp * (ep + pm);

    return {eq1, eq2};
}

double FluidProfile::matching_eqs_wall2(double vp, double TpTN, double TmTN) const {    
    if (TmTN < veff_TTN_vals_.front() || TmTN > veff_TTN_vals_.back()) {
        throw std::out_of_range("Tm/TN is called out of bounds for spline!");
    }

    const auto pp = alglib::spline1dcalc(veff_ps_interp_, TpTN); // p_+, e_+
    const auto ep = alglib::spline1dcalc(veff_es_interp_, TpTN);

    const auto pm = alglib::spline1dcalc(veff_pb_interp_, TmTN); // p_-, e_-
    const auto em = alglib::spline1dcalc(veff_eb_interp_, TmTN);

    return vp * vp - (pm - pp) * (em + pp) / ((em - ep) * (ep + pm));
}

std::vector<double> FluidProfile::get_IC_detonation_veff(double vp, double TpTN) const {
    // uses matching conditions to get vm, Tm/TN from vp, Tp/TN
    //      vp * vm = (pm - pp) / (em - ep)
    //      vp / vm = (em + pp) / (ep + pm)

    if (TpTN < veff_TTN_vals_.front() || TpTN > veff_TTN_vals_.back()) {
        throw std::out_of_range("Tp/TN is called out of bounds for spline");
    }
    // const auto pp = alglib::spline1dcalc(veff_ps_interp_, TpTN); // p_+, e_+
    // const auto ep = alglib::spline1dcalc(veff_es_interp_, TpTN);

    std::function<std::vector<double>(std::vector<double>)> matching_helper = [this, vp, TpTN] (std::vector<double> y0) {
        return matching_eqs_wall(vp, TpTN, y0[0], y0[1]); // y0 = {vm, TmTN}
    };
    
    // solve matching eqs
    const auto sol = newton_solve(matching_helper, {vp, TpTN}); // initial guess vp, TpTN

    return {sol[0], sol[1]}; // {vm, Tm/TN}
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
    auto dydxi = [this] (double xi, const state_type& y) -> state_type {
        const auto v = y[0];
        const auto w = y[1];
        const auto T = y[2];
        const auto csq = (xi < vw_) ? cmsq_ : cpsq_; // change to csq(T) when implemented

        return { dvdxi(xi, v, csq), dwdxi(xi, v, w, csq), dTdxi(xi, v, T, csq) };
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

        const auto xi_sh = (eos_ == "veff") ? find_shock_veff(dydxi) : find_shock(dydxi);
        const auto v1UF = v1UF_from_shock(xi_sh); // only valid for bag
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

        // initial conditions:
        //     v0 = v(xi_w) = vm(UF)
        //     w0 = w(xi_w) = wm/wN
        //     T0 = T(xi_w) = Tm/TN
        const auto vp = -vw_;
        const auto wpwN = 1.0; // w+ = wN
        const auto TpTN = 1.0; // T+ = TN
        double vm, wmwN, TmTN;

        if (eos_ == "bag") {
            const auto alp = alN_; // alpha_+ = alpha_N

            vm = calc_vm(vp, alp);
            wmwN = calc_wm(wpwN, vp, vm);
            TmTN = calc_ToTN(wmwN, cpsq_, cmsq_);
        } else {
            const auto ics = get_IC_detonation_veff(vp, TpTN);
            vm = ics[0];
            TmTN = ics[1];
            // wmwN = 1.45995;
            wmwN = calc_wm(wpwN, vp, vm);
        }

        // do better error handling later
        if (vm >= vw_) {
            throw std::invalid_argument("vm must be < vw for detonation");
        }
        if (wmwN <= wpwN) {
            throw std::invalid_argument("wm/wN must be > wp/wN for detonation");
        }
        if (TmTN <= TpTN) {
            throw std::invalid_argument("Tm/TN must be > Tp/TN for detonation");
        }

        const auto vmUF = mu(vw_, abs(vm));
        y0_.push_back(vmUF);
        y0_.push_back(wmwN);
        y0_.push_back(TmTN);

        std::cout << "v0=" << y0_[0] << ", w0=" << y0_[1] << ", T0=" << y0_[2] << "\n";

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

            // la_sol_tmp.push_back(get_la_behind_wall(w_sol_tmp[i]));
        }

        // tmp for testing - use better implementation later
        if (eos_ == "bag") {
            for (size_t i = 0; i < xi_sol_tmp.size(); i++) {
                la_sol_tmp.push_back(get_la_behind_wall(w_sol_tmp[i]));
            }
        } else {
            la_sol_tmp = get_lambda(T_sol_tmp);
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