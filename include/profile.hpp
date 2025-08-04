// profile.hpp
#ifndef INCLUDE_PROFILE_HPP_H
#define INCLUDE_PROFILE_HPP_H

#include <vector>
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

using state_type = std::vector<double>;

/**
 * @class FluidProfile
 * @brief Represents the hydrodynamic profile of a bubble wall in a first-order phase transition.
 */
class FluidProfile {
  public:
    FluidProfile(const PhaseTransition::PTParams& params, const state_type& T_vals, const state_type& pres_vals, const state_type& edens_vals, const size_t n=5000);
    FluidProfile(const PhaseTransition::PTParams& params, const size_t n=5000); // ctor
    // write dflt ctor?

    // are these needed?
    PhaseTransition::PTParams params() const { return params_; }; // PT parameters
    std::string eos() const { return eos_; }; // Equation of state (bag or veff)

    double xi0() const { return xi0_; };
    double xif() const { return xif_; };
    std::vector<double> init_state() const { return y0_; }; // Initial state vector {xi0, v0, w0}
    
    state_type xi_vals() const { return xi_vals_; }; // Vector of xi=r/t
    state_type v_vals() const { return v_vals_; }; // v(xi)
    state_type w_vals() const { return w_vals_; }; // w(xi)
    state_type la_vals() const { return la_vals_; }; // la(xi)

    void write(const std::string& filename = "bubble_prof.csv") const; // write bubble profile to disk
    
    #ifdef ENABLE_MATPLOTLIB
    void plot(const std::string& filename = "bubble_prof.png") const; // Plots bubble profiles
    #endif

  private:
    std::string eos_; // bag or veff

    const PhaseTransition::PTParams params_; // local copy of PT parameters
    const double cpsq_, cmsq_, vw_, alN_;
    double alp_min_, alp_max_;
    
    int mode_; // hydrodynamic mode (deflagration=0, hybrid=1, detonation=2)
    double xi0_, xif_; // initial/final xi
    std::vector<double> y0_; // initial conditions {v0, w0}

    // Veff stuff
    const state_type veff_T_vals_, veff_p_vals_, veff_e_vals_, veff_w_vals_; // T, P(T), e(T) for fluid profile using generic EoS
    alglib::spline1dinterpolant ToTN_interp_;
    
    state_type xi_vals_, v_vals_, w_vals_, T_vals_, la_vals_; // xi, v(xi), w(xi), la(x)

    int get_mode(double vw, double cmsq, double alN) const;
    double vJ_det(double alp) const;

    double calc_vm(double vp, double alpha_p) const;
    double calc_vp(double vm, double alpha_p) const;
    double calc_wm(double wp, double vp, double vm) const;
    double calc_w1wN(double xi_sh) const;
    double calc_ToTN(double wmwN, double cpsq, double cmsq) const;

    double xi_shock(double v1UF) const; // position of shock front
    double v1UF_from_shock(double xi_sh) const;
    std::vector<double> get_alp_minmax(double vw) const;
    double get_alp_wall(double vpUF, double vw) const;

    double alN_residual_func(double xi_sh, const deriv_func& dydxi) const;

    double get_la_behind_wall(double w) const;
    double get_la_front_wall(double w) const;

    double find_shock(const deriv_func& dydxi) const;
    // put number of integration points in input file? seems bad to hardcode
    std::vector<state_type> solve_profile(int n=100);
    // state_type calc_lambda_vals(state_type w_vals) const;

    // testing purposes ONLY
    std::vector<state_type> read(const std::string& filename) const; // read bubble profile from disk

};

} // namespace Hydrodynamics

#endif // INCLUDE_PROFILE_HPP_H
