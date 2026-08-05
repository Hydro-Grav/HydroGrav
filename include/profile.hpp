/**
 * @file profile.hpp
 * @brief Constructs fluid profiles using the bag/mu-nu equations of state
 *        and from generic Veff(T).
 *
 * This header provides the 'FluidProfile' class to construct fluid
 * profiles by solving the fluid equations of motion and enforcing 
 * matching conditions at the wall/shock.
 */
#ifndef INCLUDE_PROFILE_HPP_H
#define INCLUDE_PROFILE_HPP_H

#include <vector>
#include <array>
#include <string>
#include <variant>

#include "ap.h"
#include "interpolation.h"

#include "maths.hpp"
#include "phasetransition.hpp"

/**
 * @namespace Hydrodynamics
 * @brief Hydrodynamic calculations for FOPT.
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

/**
 * @brief Derivative \f$d\xi/dv\f$ along a fluid profile characteristic.
 *
 * @param xi  Self-similar coordinate \f$\xi = r/t\f$.
 * @param v   Fluid velocity in the universe frame.
 * @param csq Local speed of sound squared.
 * @return \f$d\xi/dv\f$.
 */
double dxidv(double xi, double v, const double csq);

/**
 * @brief Derivative \f$d\xi/dw\f$ along a fluid profile characteristic.
 *
 * @param xi  Self-similar coordinate.
 * @param v   Fluid velocity.
 * @param w   Enthalpy density.
 * @param csq Local speed of sound squared.
 * @return \f$d\xi/dw\f$.
 */
double dxidw(double xi, double v, double w, const double csq);

/**
 * @brief Derivative \f$d\xi/dT\f$ along a fluid profile characteristic.
 *
 * @param xi  Self-similar coordinate.
 * @param v   Fluid velocity.
 * @param T   Temperature.
 * @param csq Local speed of sound squared.
 * @return \f$d\xi/dT\f$.
 */
double dxidT(double xi, double v, double T, const double csq);

/**
 * @brief Equations of motion for the fluid profile ODE system.
 *
 * Returns the derivative vector \f$(d\xi/dv,\, dw/dv,\, dT/dv)\f$ at
 * the current state, suitable for passing to an RK4 solver.
 *
 * @param v     Current velocity (independent variable).
 * @param y     State vector \f$(\xi, w, T)\f$.
 * @param vw    Wall velocity.
 * @param cmsq  Speed of sound squared in the broken phase.
 * @param cpsq  Speed of sound squared in the symmetric phase.
 * @return Derivative state vector.
 */
state_type dydv_vec(double v, const state_type& y, double vw, double cmsq, double cpsq);

/**
 * @class FluidProfile
 * @brief Represents the hydrodynamic profile of a bubble wall in a first-order phase transition.
 *
 * Given a set of phase transition parameters, this class solves the relativistic
 * fluid equations to obtain the self-similar velocity \f$v(\xi)\f$, enthalpy
 * \f$w(\xi)\f$, and temperature \f$T(\xi)\f$ profiles across the bubble wall.
 * Both bag and effective-potential (Veff) equations of state are supported.
 */
class FluidProfile {
  public:
    /**
     * @brief Construct a fluid profile by solving the hydrodynamic ODEs.
     *
     * @param params  Phase transition parameters (bag or Veff EoS).
     * @param n       Number of integration steps (default: 5000).
     * @param dev_log Enable verbose developer logging (default: false).
     */
    FluidProfile(const PhaseTransition::PTParams& params, const size_t n=5000, const bool dev_log=false);

    /// Pointer to the phase transition parameters used to construct this profile.
    const PhaseTransition::PTParams* params() const { return params_; }; // PT parameters

    // unused?
    std::string eos() const { return eos_; }; // Equation of state (bag or veff)
    
    /// Vector of self-similar coordinate values \f$\xi = r/t\f$.
    prof_type xi_vals() const { return xi_vals_; };
    /// Fluid velocity profile \f$v(\xi)\f$.
    prof_type v_vals() const { return v_vals_; };
    /// Enthalpy density profile \f$w(\xi)\f$ normalised to the nucleation value.
    prof_type w_vals() const { return w_vals_; };
    /// Lambda profile \f$\lambda(\xi)\f$.
    prof_type la_vals() const { return la_vals_; };
    /// Temperature profile \f$T(\xi)/T_N\f$.
    prof_type T_vals() const { return T_vals_; }

    /// Outer boundary of the profile integration domain (shock position \f$\xi_{\rm sh}\f$).
    double xi_max() const { return xi_max_integrate_; } // xi_sh
    /// Inner boundary of the integration domain (\f$v_w\f$ for deflagrations, \f$c_-\f$ for hybrids/detonations).
    double xi_min() const { return xi_min_integrate_; } // vw (def), cm (det/hyb)

    double cpsq() const { return cpsq_; }
    double cmsq() const { return cmsq_; }

    /**
     * @brief Integer code identifying the hydrodynamic mode.
     * @return 0 = deflagration, 1 = hybrid, 2 = detonation.
     */
    int mode() const { return mode_; }; // Hydrodynamic mode (0=deflagration, 1=hybrid, 2=detonation)
    /// String for the hydrodynamic mode.
    std::string mode_str() const;

    /**
     * @brief Write the fluid profile to a CSV file.
     * @param filename Output filename (default: "fp.csv").
     */
    void write(const std::string& filename = "fp.csv") const; // write bubble profile to disk

    /**
     * @brief Flag indicating whether the shock front converged during integration.
     *
     * If @c false the profile uses mu-nu shock matching condition as a fallback.
     */
    bool shock_flag() const { return shock_flag_; }
    
    #ifdef ENABLE_MATPLOTLIB
    void plot(const std::string& filename = "fp.png", double xi_min = 0.0, double xi_max = 1.0) const; // Plots bubble profiles
    #endif

  private:
    const size_t n_;                                              ///< Number of points used to construct profile
    const PhaseTransition::PTParams* params_;                     ///< Pointer to the phase transition parameters.
    const PhaseTransition::PTParams_Bag* bag_params_ = nullptr;   ///< Downcast pointer for bag EoS (null if Veff).
    const PhaseTransition::PTParams_Veff* veff_params_ = nullptr; ///< Downcast pointer for Veff EoS (null if bag).

    std::string eos_; ///< Equation of state identifier: "bag" or "veff".
    
    double cpsq_;       ///< Speed of sound squared in the symmetric phase.
    double cmsq_;       ///< Speed of sound squared in the broken phase.
    const double vw_;   ///< Wall velocity.
    const double alN_;  ///< Phase transition strength parameter at the nucleation temperature.
    double alp_min_;    ///< Minimum physical value of \f$\alpha_+\f$ (lower bound for root search).
    double alp_max_;    ///< Maximum physical value of \f$\alpha_+\f$ (upper bound for root search).
    int mode_;          ///< Hydrodynamic mode: 0 = deflagration, 1 = hybrid, 2 = detonation.
    prof_type xi_vals_, v_vals_, w_vals_, T_vals_, la_vals_; ///< Profile arrays: \f$\xi\f$, \f$v(\xi)\f$, \f$w(\xi)\f$, \f$T(\xi)\f$, \f$\lambda(\xi)\f$.
    double xi_min_integrate_; ///< Lower limit of the profile integration domain.
    double xi_max_integrate_; ///< Upper limit of the profile integration domain (shock position).

    bool shock_flag_;    ///< Set to @c false if shock convergence failed during Veff integration.
    const bool dev_log_; ///< If @c true, emit verbose diagnostic output during profile construction.

    /************************** Bag EoS **************************/

    /**
     * @brief Determine the hydrodynamic mode for the bag/mu-nu equation of state.
     *
     * Compares the wall velocity against the Jouguet velocity to classify
     * the solution as deflagration, hybrid, or detonation.
     *
     * @param vw   Wall velocity.
     * @param cmsq Speed of sound squared in the broken phase.
     * @param alN  Phase transition strength parameter.
     * @return Integer mode code: 0 = deflagration, 1 = hybrid, 2 = detonation.
     */
    int get_mode_bag(double vw, double cmsq, double cpsq, double alN) const;

    /**
     * @brief Compute the Jouguet detonation velocity for a given \f$\alpha_+\f$.
     *
     * @param alp Strength parameter \f$\alpha_+\f$ immediately ahead of the wall.
     * @return Jouguet velocity \f$v_J(\alpha_+)\f$.
     */
    double vJ_det(double alp) const;

    double vJ_inv_det(double alp) const;

    /**
     * @brief Compute the physically allowed range of \f$\alpha_+\f$ for a given wall velocity.
     *
     * @param vw Wall velocity.
     * @return Array {alpha_min, alpha_max}.
     */
    std::array<double, 2> get_alp_minmax(double vw) const;

    /**
     * @brief Compute \f$\alpha_+\f$ at the wall from the fluid velocities on each side.
     *
     * Uses the bag-model matching conditions to infer the local strength
     * parameter from the velocities immediately ahead (\f$v_+\f$) and
     * behind (\f$v_-\f$) the wall.
     *
     * @param vp Fluid velocity in the symmetric phase just ahead of the wall (wall frame).
     * @param vm Fluid velocity in the broken phase just behind the wall (wall frame).
     * @return \f$\alpha_+\f$.
     */
    double alp_from_matching(double vp, double vm) const;

    /**
     * @brief Compute \f$v_+\f$ from \f$v_-\f$ and \f$\alpha_+\f$ via the bag matching conditions.
     *
     * @param vm      Fluid velocity behind the wall (wall frame).
     * @param alp Strength parameter \f$\alpha_+\f$.
     * @return \f$v_+\f$.
     */
    double vp_from_matching(double vm, double alp, int sgn) const;

    /**
     * @brief Compute \f$v_-\f$ from \f$v_+\f$ and \f$\alpha_+\f$ via the bag matching conditions.
     *
     * @param vp      Fluid velocity ahead of the wall (wall frame).
     * @param alp Strength parameter \f$\alpha_+\f$.
     * @return \f$v_-\f$.
     */
    double vm_from_matching(double vp, double alp, int sgn) const;

    /**
     * @brief Compute the temperature ratio \f$T_-/T_N\f$ from the enthalpy ratio \f$w_-/w_N\f$.
     *
     * @param wmwN Ratio of enthalpy density behind the wall to the nucleation value.
     * @return \f$T_-/T_N\f$.
     */
    double TmTN_from_matching(double wmwN) const;
    double T1TN_from_matching(double w1wN, double w2wN, double T2TN) const;

    /**
     * @brief Compute the enthalpy ratio and temperature immediately behind shock.
     *
     * Applies the Rankine–Hugoniot shock-matching conditions at shock
     * position \f$\xi_{\rm sh}\f$.
     *
     * @param xi_sh Self-similar shock position.
     * @return Pair {w_shock/w_N, T_shock/T_N}.
     */
    std::pair<double, double> wT_from_shock(double xi_sh) const;

    /**
     * @brief Find the fluid velocity \f$v_+^{\rm UF}\f$ (universe frame) ahead of the wall
     *        for a deflagration, using a root search over the \f$\alpha_N\f$ residual.
     *
     * @param dydv Fluid ODE right-hand side.
     * @param n    Number of integration steps (default: 1000).
     * @return \f$v_+^{\rm UF}\f$.
     */
    double find_vpUF(const deriv_func& dydv, const size_t n=1000) const;

    /**
     * @brief Integrate the deflagration profile and locate the shock front.
     *
     * Integrates the fluid ODEs from the wall outward and identifies
     * the shock position by monitoring the Rankine–Hugoniot residual.
     *
     * @param dydv      Fluid ODE right-hand side.
     * @param vpUF      Universe-frame fluid velocity ahead of the wall.
     * @param test_resi If @c true, print diagnostic residual information.
     * @param n         Number of ODE integration steps.
     * @param tol       Shock-finding convergence tolerance.
     * @return Tuple {v_solution, y_solution, shock_found_flag}.
     */
    std::tuple<std::vector<double>, std::vector<state_type>, bool> deflagration_profile(const deriv_func& dydv, double vpUF, const bool test_resi=false, const size_t n=1000) const;

    /**
     * @brief Compute the \f$\alpha_N\f$ residual for a deflagration given trial initial conditions.
     *
     * Used internally as the objective function for the root search in
     * @c find_vpUF().
     *
     * @param dydv  Fluid ODE right-hand side.
     * @param vpUF  Trial universe-frame velocity ahead of the wall.
     * @param vm    Fluid velocity behind the wall (wall frame).
     * @param n     Number of integration steps.
     * @return Residual \f$\alpha_N^{\rm computed} - \alpha_N^{\rm target}\f$.
     */
    double alN_residual(const deriv_func& dydv, double vpUF, double vm, const size_t n=1000) const;

    /**
     * @brief Integrate the deflagration profile from explicit internal initial conditions.
     *
     * @param dydv      Fluid ODE right-hand side.
     * @param vpUF      Universe-frame velocity ahead of the wall.
     * @param wpwN      Enthalpy ratio \f$w_+/w_N\f$.
     * @param TpTN      Temperature ratio \f$T_+/T_N\f$.
     * @param test_resi If @c true, print diagnostic residual information.
     * @param n         Number of integration steps.
     * @return Tuple {v_solution, y_solution, shock_found_flag}.
     */
    std::tuple<std::vector<double>, std::vector<state_type>, bool> deflagration_profile_internal(const deriv_func& dydv, double vpUF, double wpwN, double TpTN, const bool test_resi=false, const size_t n=1000) const;

    /**
     * @brief Locate the shock front index in a solved deflagration profile.
     *
     * Scans the integrated solution and returns the index at which the
     * Rankine–Hugoniot shock conditions are first satisfied.
     *
     * @param v_sol     Vector of velocity values along the ODE solution.
     * @param y_sol     Corresponding state vectors.
     * @param test_resi If @c true, print diagnostic information.
     * @return Index of the shock in @p v_sol / @p y_sol.
     */
    size_t find_shock_idx(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol, const bool test_resi=false) const;
    size_t find_shock_idx_inv(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol) const;

    /** @brief Developer diagnostic: print the \f$\alpha_N\f$ residual over a range of \f$v_+^{\rm UF}\f$ values. */
    void test_alN_residual(const deriv_func& dydv, double vm, const size_t n) const;
    /** @brief Developer diagnostic: verify that the bag-model shock conditions are satisfied. */
    void test_shock_bag(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol) const;
    void test_shock_bag_inv(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol) const;

    /**
     * @brief Compute the initial conditions for a bag/mu-nu model detonation profile.
     *
     * Returns the self-similar starting coordinate and state vector
     * \f$(\xi_0, v_0, w_0, T_0)\f$ immediately behind the detonation front.
     *
     * @return Pair {xi_0, state_0}.
     */
    std::pair<double, state_type> get_IC_detonation() const;

    /**
     * @brief Lambda profile \f$\lambda_b\f$ in the broken phase (bag/mu-nu EoS).
     * @param wowN Enthalpy ratio \f$w/w_N\f$.
     * @return \f$\lambda_b(w/w_N)\f$.
     */
    double lambda_b(double wowN) const;

    /**
     * @brief Lambda profile \f$\lambda_s\f$ in the symmetric phase (bag/mu-nu EoS).
     * @param wowN Enthalpy ratio \f$w/w_N\f$.
     * @return \f$\lambda_s(w/w_N)\f$.
     */
    double lambda_s(double wowN) const;

    /************************** Veff EoS *************************/

    /**
     * @brief Determine the hydrodynamic mode for the Veff equation of state.
     *
     * @param vw   Wall velocity.
     * @param cmsq Speed of sound squared in the broken phase evaluated at \f$T_N\f$.
     * @return Integer mode code: 0 = deflagration, 1 = hybrid, 2 = detonation.
     */
    int get_mode_veff(double vw, double cmsq) const;

    /**
     * @brief Evaluate the Veff wall-matching equations.
     *
     * Returns the two-component residual of the junction conditions
     * \f$(p_+ = p_-,\; T_+ v_+ = T_- v_-)\f$ at the bubble wall.
     *
     * @param vp   Fluid velocity in the symmetric phase (wall frame).
     * @param TpTN Temperature ratio \f$T_+/T_N\f$.
     * @param vm   Fluid velocity in the broken phase (wall frame).
     * @param TmTN Temperature ratio \f$T_-/T_N\f$.
     * @return Two-component residual vector.
     */
    std::array<double, 2> matching_eqs_wall(double vp, double TpTN, double vm, double TmTN) const;

    /**
     * @brief Evaluate the Veff shock-matching equations.
     *
     * Returns the two-component residual of the Rankine–Hugoniot conditions
     * across the shock front.
     *
     * @param v1   Fluid velocity ahead of the shock (universe frame).
     * @param T1TN Temperature ratio ahead of the shock.
     * @param v2   Fluid velocity behind the shock (universe frame).
     * @param T2TN Temperature ratio behind the shock.
     * @return Two-component residual vector.
     */
    std::array<double, 2> matching_eqs_shock(double v1, double T1TN, double v2, double T2TN) const;

    /**
     * @brief Find \f$T_-/T_N\f$ for a Veff deflagration by root-finding on the shock residual.
     *
     * @param dydv    Fluid ODE right-hand side.
     * @param fallback If @c true, use a fallback bisection algorithm.
     * @return \f$T_-/T_N\f$ at the wall.
     */
    double find_TmTN_veff(const deriv_func& dydv, const bool fallback=false) const;

    /**
     * @brief Compute the universe-frame velocity residual at the shock for a given \f$T_-/T_N\f$.
     *
     * Used as the objective in the root search performed by @c find_TmTN_veff().
     *
     * @param dydv    Fluid ODE right-hand side.
     * @param TmTN    Trial temperature ratio \f$T_-/T_N\f$.
     * @param fallback If @c true, use the fallback shock-finding algorithm.
     * @param n       Number of ODE integration steps.
     * @return Residual of the universe-frame velocity at the shock.
     */
    double v1_residual_veff(const deriv_func& dydv, double TmTN, const bool fallback=false, const size_t n=1000) const;

    /**
     * @brief Solve the Veff wall-matching conditions given the broken-phase state.
     *
     * Given \f$v_-\f$, \f$T_-/T_N\f$, and \f$w_-/w_N\f$, determines the
     * symmetric-phase state \f$(v_+, T_+/T_N, w_+/w_N)\f$ by satisfying the
     * junction conditions.
     *
     * @param vm   Fluid velocity behind the wall (wall frame).
     * @param TmTN Temperature ratio \f$T_-/T_N\f$.
     * @param wmwN Enthalpy ratio \f$w_-/w_N\f$.
     * @return Tuple {v_+, T_+/T_N, w_+/w_N}.
     */
    std::tuple<double, double, double> wall_matching_veff(const double vm, const double TmTN, const double wmwN) const;

    /**
     * @brief Locate the shock front index in a Veff deflagration solution.
     *
     * @param v_sol   Integrated velocity values.
     * @param y_sol   Integrated state vectors.
     * @param fallback If @c true, use the fallback shock-detection criterion.
     * @return Index of the shock in the solution arrays.
     */
    size_t find_shock_idx_veff(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol, const bool fallback=false) const;

    /**
     * @brief Integrate the full Veff deflagration profile from the wall to the shock.
     *
     * @param dydv      Fluid ODE right-hand side.
     * @param vm        Fluid velocity behind the wall (wall frame).
     * @param wmwN      Enthalpy ratio \f$w_-/w_N\f$.
     * @param TmTN      Temperature ratio \f$T_-/T_N\f$.
     * @param fallback  If @c true, use fallback algorithms for shock finding.
     * @param test_resi If @c true, print diagnostic residual information.
     * @param n         Number of ODE integration steps.
     * @return Pair {v_solution, y_solution}.
     */
    std::pair<std::vector<double>, std::vector<state_type>> deflagration_profile_veff(const deriv_func& dydv, double vm, double wmwN, double TmTN, const bool fallback=false, const bool test_resi=false, const size_t n=1000) const;

    /**
     * @brief Check whether the Veff shock conditions have converged to the required tolerance.
     *
     * Updates @c shock_flag_ as a side-effect.
     *
     * @param v_sol Integrated velocity values.
     * @param y_sol Integrated state vectors.
     * @return @c true if the shock is converged.
     */
    bool check_shock_convergence(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol);

    /**
     * @brief Fallback convergence check with a relaxed tolerance.
     *
     * @param v_sol Integrated velocity values.
     * @param y_sol Integrated state vectors.
     * @return @c true if the shock satisfies the fallback criterion.
     */
    bool check_shock_convergence_fallback(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol) const;

    /** @brief Developer diagnostic: print the Veff velocity residual over a scan of \f$T_-/T_N\f$. */
    void test_residual_veff(const deriv_func& dydv, const size_t n) const;
    /** @brief Developer diagnostic: verify that the Veff shock conditions are satisfied at the detected shock. */
    void test_shock_veff(const std::vector<double>& v_sol, const std::vector<state_type>& y_sol) const;

    /**
     * @brief Compute the initial conditions for a Veff detonation profile.
     *
     * Searches for \f$v_-\f$ in [vm_min, vm_max] satisfying the wall
     * matching conditions, then returns the corresponding initial state.
     *
     * @param vm_min Lower bound of the \f$v_-\f$ search interval.
     * @param vm_max Upper bound of the \f$v_-\f$ search interval.
     * @return Pair {xi_0, state_0}.
     */
    std::pair<double, state_type> get_IC_detonation_veff(const double vm_min=0.0, const double vm_max=1.0) const;

    /**
     * @brief Lambda profile \f$\lambda_s\f$ in the symmetric phase (Veff EoS).
     * @param ToTN    Temperature ratio \f$T/T_N\f$.
     * @param eN      Energy density at the nucleation temperature.
     * @param wN_inv  Inverse of the enthalpy density at nucleation.
     * @return \f$\lambda_s(T/T_N)\f$.
     */
    double lambda_s_veff(double ToTN, const double eN, const double wN_inv) const;

    /**
     * @brief Lambda profile \f$\lambda_b\f$ in the broken phase (Veff EoS).
     * @param ToTN    Temperature ratio \f$T/T_N\f$.
     * @param eN      Energy density at the nucleation temperature.
     * @param wN_inv  Inverse of the enthalpy density at nucleation.
     * @return \f$\lambda_b(T/T_N)\f$.
     */
    double lambda_b_veff(double ToTN, const double eN, const double wN_inv) const;

    /********************** Both Bag/Veff EoS ********************/

    /**
     * @brief Compute the enthalpy ratio \f$w_-/w_N\f$ behind the wall from the matching conditions.
     *
     * Uses the junction condition \f$w_- = w_+ v_+ / v_-\f$ (in the wall frame),
     * valid for both bag and Veff equations of state.
     *
     * @param wp Enthalpy density ahead of the wall (wall frame).
     * @param vp Fluid velocity ahead of the wall (wall frame).
     * @param vm Fluid velocity behind the wall (wall frame).
     * @return Enthalpy density behind the wall.
     */
    double w_from_matching(double wp, double vp, double vm) const;

    double find_vp_inv_det(const deriv_func& dydv) const;
    double inv_det_residual(const deriv_func& dydv, const double vpUF) const;
    void test_inv_det_residual(const deriv_func& dydv, const size_t n) const;
    double wpwN_from_alp(const double alp) const;

    std::vector<prof_type> clean_profiles(const prof_type& v_sol_tmp, const std::vector<state_type>& y_sol_tmp,
                                          const double w_start_val, const double w_end_val,
                                          const double T_start_val, const double T_end_val) const;

    /**
     * @brief Solve for the full bag/mu-nu model fluid profile and return all profile arrays.
     *
     * @param n Number of integration steps (default: 5000).
     * @return Vector of profile arrays {xi, v, w, T, lambda}.
     */
    std::vector<prof_type> solve_profile(int n=5000);

    /**
     * @brief Solve for the full Veff fluid profile and return all profile arrays.
     *
     * @param n Number of integration steps (default: 5000).
     * @return Vector of profile arrays {xi, v, w, T, lambda}.
     */
    std::vector<prof_type> solve_profile_veff(int n=5000);

    std::vector<prof_type> solve_inverse_profile(int n=5000);
};

#ifdef ENABLE_MATPLOTLIB
void plot_profiles(const FluidProfile& fp_bag, const FluidProfile& fp_munu, const FluidProfile& fp_veff, const std::string& filename="fp_combined.png", const double xi_min=0.0, const double xi_max=1.0);
#endif

} // namespace Hydrodynamics

#endif // INCLUDE_PROFILE_HPP_H
