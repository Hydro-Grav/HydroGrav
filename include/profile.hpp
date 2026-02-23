// profile.hpp
#ifndef INCLUDE_PROFILE_HPP_H
#define INCLUDE_PROFILE_HPP_H

#include <vector>
#include <array>
#include <string>
#include <variant>

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

using prof_type = std::vector<double>;
using state_type = std::array<double, 3>;
using deriv_func = std::function<std::array<double, 3>(double, const std::array<double, 3>&)>;

/**
 * @brief Computes the Lorentz factor between the wall frame and the universe frame.
 *
 * @param xi Fluid velocity in the wall frame.
 * @param v Fluid velocity in the universe frame.
 * 
 * @return Lorentz factor (mu).
 */
double mu(double xi, double v);

double dxidv(double xi, double v, const double csq);
double dxidw(double xi, double v, double w, const double csq);
double dxidT(double xi, double v, double T, const double csq);
state_type dydv_vec(double v, const state_type& y, double vw, double cmsq, double cpsq);

void generate_streamplot_data(const PhaseTransition::PTParams& params, int xi_pts, int y_pts, const std::string& filename);
void generate_streamplot_data(const PhaseTransition::PTParams& params);

/**
 * @class FluidProfile
 * @brief Represents the hydrodynamic profile of a bubble wall in a first-order phase transition.
 */
class FluidProfile {
  public:
    FluidProfile(const PhaseTransition::PTParams& params, const size_t n=5000);

    // return pointer to base class
    const PhaseTransition::PTParams* params() const { return params_; }; // PT parameters

    // unused?
    std::string eos() const { return eos_; }; // Equation of state (bag or veff)
    
    prof_type xi_vals() const { return xi_vals_; }; // Vector of xi=r/t
    prof_type v_vals() const { return v_vals_; }; // v(xi)
    prof_type w_vals() const { return w_vals_; }; // w(xi)
    prof_type la_vals() const { return la_vals_; }; // la(xi)
    prof_type T_vals() const { return T_vals_; }

    double xi_max() const { return xi_max_integrate_; } // xi_sh
    double xi_min() const { return xi_min_integrate_; } // vw (def), cm (det/hyb)

    int mode() const { return mode_; }; // Hydrodynamic mode (0=deflagration, 1=hybrid, 2=detonation)
    std::string mode_str() const;

    void write(const std::string& filename = "fp.csv") const; // write bubble profile to disk
    
    #ifdef ENABLE_MATPLOTLIB
    void plot(const std::string& filename = "fp.png", double xi_min = 0.0, double xi_max = 1.0) const; // Plots bubble profiles
    #endif

  private:
    const PhaseTransition::PTParams* params_; // local copy of PT parameters  
    const PhaseTransition::PTParams_Bag* bag_params_ = nullptr;
    const PhaseTransition::PTParams_Veff* veff_params_ = nullptr;

    std::string eos_; // bag or veff
    
    const double cpsq_, cmsq_, vw_, alN_;
    double alp_min_, alp_max_;
    int mode_; // hydrodynamic mode (deflagration=0, hybrid=1, detonation=2)
    prof_type xi_vals_, v_vals_, w_vals_, T_vals_, la_vals_; // xi, v(xi), w(xi), la(x)
    double xi_min_integrate_, xi_max_integrate_; // start/endpoints of profile for integration

    /************************** Bag EoS **************************/
    int get_mode_bag(double vw, double cmsq, double alN) const;
    double vJ_det(double alp) const;
    std::array<double, 2> get_alp_minmax(double vw) const;

    // matching at wall
    // double get_alp_wall(double vpUF, double vw) const;
    double get_alp_wall(double vp, double vm) const;
    double vp_from_matching(double vm, double alpha_p) const;
    double vm_from_matching(double vp, double alpha_p) const;
    double get_TmTN(double wmwN) const;

    // matching at shock
    // double v1UF_from_shock(double xi_sh) const;
    std::pair<double, double> wT_from_shock(double xi_sh) const;

    // deflagrations
    double find_vpUF(const deriv_func& dydv, const size_t n=1000) const;
    std::pair<std::vector<double>, std::vector<state_type>> deflagration_profile(const deriv_func& dydv, double vpUF, const bool test_resi=false, const size_t n=1000) const;
    double alN_residual(const deriv_func& dydv, double vpUF, double vm, const size_t n=1000) const;
    std::pair<std::vector<double>, std::vector<state_type>> deflagration_profile_internal(const deriv_func& dydv, double vpUF, double wpwN, double TpTN, const bool test_resi=false, const size_t n=1000) const;
    size_t find_shock_idx(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol, const bool test_resi=false, const double tol=1e-5) const;

    // dev
    void test_alN_residual(const deriv_func& dydv, double vm, const size_t n) const;
    void test_shock_bag(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol) const;

    // detonations
    std::pair<double, state_type> get_IC_detonation() const;

    // lambda profiles
    double lambda_b(double wowN) const;
    double lambda_s(double wowN) const;

    /************************** Veff EoS *************************/  
    int get_mode_veff(double vw, double cmsq) const;

    // matching eqs
    std::array<double, 2> matching_eqs_wall(double vp, double TpTN, double vm, double TmTN) const;
    std::array<double, 2> matching_eqs_shock(double v1, double T1TN, double v2, double T2TN) const;
    
    // deflagrations
    double find_TmTN_veff(const deriv_func& dydv) const;
    double T2TN_residual_veff(const deriv_func& dydv, double TmTN, const size_t n=1000) const;
    std::pair<std::vector<double>, std::vector<state_type>> deflagration_profile_veff(const deriv_func& dydv, double vm, double wmwN, double TmTN, const bool test_resi=false, const size_t n=1000) const;
    size_t find_shock_idx_veff(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol, const bool test_resi=false) const;
    
    // dev
    void test_residual_veff(const deriv_func& dydv, const size_t n) const;
    void test_shock_veff(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol) const;
    
    // detonations
    std::pair<double, state_type> get_IC_detonation_veff() const;

    // lambda profiles
    double lambda_s_veff(double ToTN, const double eN, const double wN_inv) const;
    double lambda_b_veff(double ToTN, const double eN, const double wN_inv) const;

    /********************** Both Bag/Veff EoS ********************/  
    double w_from_matching(double wp, double vp, double vm) const;

    // solve profiles
    std::vector<prof_type> solve_profile(int n=5000);
    std::vector<prof_type> solve_profile_veff(int n=5000);
};

#ifdef ENABLE_MATPLOTLIB
void plot_profiles(const FluidProfile& fp_bag, const FluidProfile& fp_munu, const FluidProfile& fp_veff, const std::string& filename="fp_combined.png", const double xi_min=0.0, const double xi_max=1.0);
#endif

} // namespace Hydrodynamics

#endif // INCLUDE_PROFILE_HPP_H
