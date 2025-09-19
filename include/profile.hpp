// profile.hpp
#ifndef INCLUDE_PROFILE_HPP_H
#define INCLUDE_PROFILE_HPP_H

#include <vector>
#include <array>
#include <string>

#include "ap.h"
#include "interpolation.h"

#include "maths_ops.hpp"
#include "phasetransition.hpp"

/*
TO DO:
- move generate_stream_plot() elsewhere? independent of initial conditions so the same for all instances of FluidProfile
- add some safety thing for private vals and prof in FluidProfile to stop them from changing (helper static function to precompute in initialiser list)
- add dflt ctor for FluidProfile that just uses dflt ctor for PTParams
- add error handling if P(T) or e(T) called using FluidProfile ctor for Bag model
- calculate cpsq, cmsq internally in FluidProfile (rather than as input) when passing in p(T) and e(T) to avoid wrong input
*/

namespace Hydrodynamics {

/**
 * @brief Computes the Lorentz factor between the wall frame and the universe frame.
 *
 * @param xi Fluid velocity in the wall frame.
 * @param v Fluid velocity in the universe frame.
 * 
 * @return Lorentz factor (mu).
 */
double mu(double xi, double v);

double dvdxi(double xi, double v, const double csq);
double dwdxi(double xi, double v, double w, const double csq);
double dTdxi(double xi, double v, double T, const double csq);

double dxidv(double xi, double v, const double csq);
double dxidw(double xi, double v, double w, const double csq);
double dxidT(double xi, double v, double T, const double csq);

void generate_streamplot_data(const PhaseTransition::PTParams& params, int xi_pts, int y_pts, const std::string& filename);
void generate_streamplot_data(const PhaseTransition::PTParams& params);

// struct FluidState {
//     double v;   // velocity
//     double w;   // enthalpy
//     // double T;   // temperature

//     FluidState operator+(const FluidState& other) const {
//         // return {v + other.v, w + other.w, T + other.T};
//         return {v + other.v, w + other.w};
//     }
//     FluidState operator*(double scalar) const {
//         // return {v * scalar, w * scalar, T * scalar};
//         return {v * scalar, w * scalar};
//     }
// };

// struct FluidODE {
//   double dvdxi(double xi, double v, const double csq);
//   double dwdxi(double xi, double v, double w, const double csq);

//   double dxi_dtau(double xi, double v, const double csq);
//   double dv_dtau(double xi, double v, const double csq);
//   double dw_dtau(double xi, double v, double w, const double csq);
// };

/**
 * @class FluidProfile
 * @brief Represents the hydrodynamic profile of a bubble wall in a first-order phase transition.
 */
class FluidProfile {
  public:
    FluidProfile(const PhaseTransition::PTParams& params, const size_t n=5000);
    // write dflt ctor?

    // are these needed?
    PhaseTransition::PTParams params() const { return params_; }; // PT parameters
    std::string eos() const { return eos_; }; // Equation of state (bag or veff)

    // unused
    double xi0() const { return xi0_; };
    double xif() const { return xif_; };
    
    prof_type xi_vals() const { return xi_vals_; }; // Vector of xi=r/t
    prof_type v_vals() const { return v_vals_; }; // v(xi)
    prof_type w_vals() const { return w_vals_; }; // w(xi)
    prof_type la_vals() const { return la_vals_; }; // la(xi)
    prof_type T_vals() const { return T_vals_; }

    void write(const std::string& filename = "bubble_prof.csv") const; // write bubble profile to disk
    
    #ifdef ENABLE_MATPLOTLIB
    void plot(const std::string& filename = "bubble_prof.png") const; // Plots bubble profiles
    #endif

  private:
    std::string eos_; // bag or veff

    const PhaseTransition::PTParams params_; // local copy of PT parameters
    const double cpsq_, cp_, cmsq_, cm_, vw_, alN_;
    double alp_min_, alp_max_;
    
    int mode_; // hydrodynamic mode (deflagration=0, hybrid=1, detonation=2)
    double xi0_, xif_; // initial/final xi
    
    prof_type xi_vals_, v_vals_, w_vals_, T_vals_, la_vals_; // xi, v(xi), w(xi), la(x)

    int get_mode(double vw, double cmsq, double alN) const;
    double vJ_det(double alp) const;

    double calc_vm(double vp, double alpha_p) const;
    double calc_vp(double vm, double alpha_p) const;
    double calc_wm(double wp, double vp, double vm) const;
    double calc_w1wN(double xi_sh) const;
    double calc_ToTN(double wowN, double cpsq, double cmsq) const;

    double xi_shock(double v1UF) const; // position of shock front
    double v1UF_from_shock(double xi_sh) const;
    std::array<double, 2> get_alp_minmax(double vw) const;
    double get_alp_wall(double vpUF, double vw) const;

    double alN_residual_func(double xi_sh, const deriv_func& dydxi, const int n=1000) const;
    double veff_residual_func(double xi_sh, const deriv_func& dydxi) const;

    double get_la_behind_wall(double w) const;
    double get_la_front_wall(double w) const;

    double lambda_s(double ToTN, const double eN, const double wN_inv) const;
    double lambda_b(double ToTN, const double eN, const double wN_inv) const;

    double find_shock(const deriv_func& dydxi) const;
    double find_shock_veff(const deriv_func& dydxi) const;

    state_type test_shock_matching(const deriv_func& dydxi, double xi_sh) const;
    std::pair<state_type, state_type> test_wall_matching(const deriv_func& dydxi, double xi_sh, state_type& y0, const int n=1000) const;

    double matching_residual_veff(double vp, double pp, double ep, double TmTN) const;
    std::array<double, 2> matching_eqs_wall(double vp, double TpTN, double vm, double TmTN) const;
    std::array<double, 2> matching_eqs_shock(double pN, double eN, double v2, double v1, double T1TN) const;
    double matching_eqs_shock2(double v2, double T2TN, double T1TN) const;
    double matching_eqs_wall2(double vp, double TpTN, double TmTN) const;

    void get_IC_deflagration_veff(const deriv_func& dydxi, double& xi_sh, state_type& y1, state_type& yp, state_type& ym) const;
    state_type get_IC_detonation_veff() const;

    // put number of integration points in input file? seems bad to hardcode
    std::vector<prof_type> solve_profile(int n=100);
    std::vector<prof_type> solve_profile_veff(int n=100);

    // testing purposes ONLY
    std::vector<prof_type> read(const std::string& filename) const; // read bubble profile from disk
};

} // namespace Hydrodynamics

#endif // INCLUDE_PROFILE_HPP_H
