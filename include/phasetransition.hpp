// PhaseTransition.hpp
#ifndef INCLUDE_PHASETRANSITION_HPP_H
#define INCLUDE_PHASETRANSITION_HPP_H

#include <string>
#include <iostream>
#include <vector>

#include "ap.h"
#include "interpolation.h"

#include "physics.hpp"

/*
TO DO:
- add option to input a Veff -> derive fluid dynamics from this?
- change PTParams to single ctor with default arguments (see slide 23, wk 4)
- initialise all the new variables i've added in ctor
- add deflag, hybrid, detonation type in PTParams
- add vp, vm, wp, wm as parameters in PTParams? currently uses function to calculate since they depend on Veff
  - these are calculated from a specific matching condition using bag model - generalise to other models in future
*/

namespace PhaseTransition {

// DO NOT CHANGE DEFAULT VALS
struct dflt_universe { // in units hbar = c = kB = 1
  // values today
  static constexpr double T0 = 2.34914e-13; // 2.725 K / (1.16e+13 K/GeV) = 2.349e-13 GeV
  static constexpr double H0 = 1.44328e-42; // 67.8 km/s/Mpc = 2.09502e21 s^-1 = 1.44328e-42 GeV
  static constexpr double g0 = 3.91;

  // values at start of PT
  static constexpr double Ts = 50.0; // GeV
  static constexpr double gs = 106.75;
};

/**
 * @class Universe
 * @brief Class representing the universe parameters used in the phase transition calculations.
 * 
 * This class holds the current temperature, Hubble constant, and degrees of freedom of the universe.
 */

class Universe {
  public:
    // ctors
    Universe();
    Universe(double Ts, double gs);
    Universe(double T0, double Ts, double g0, double gs, double H0);

    // params today (0) and at start of PT (s)
    double T0() const { return T0_; } // temperature of universe
    double Ts() const { return Ts_; }

    double H0() const { return H0_; } // Hubble constant
    double Hs() const { return Hs_; }

    double g0() const { return g0_; } // number of dof
    double gs() const { return gs_; }

    // print params
    void print() const;
    friend std::ostream& operator<<(std::ostream& os, const Universe& p);

  private:
    const double T0_, Ts_, g0_, gs_, H0_, Hs_;
};

// unused
const Universe& default_universe();

// DO NOT CHANGE DEFAULT VALS
struct dflt_PTParams {
  static constexpr double vw = 0.8;              // Wall velocity
  static constexpr double alN = 0.1;           // PT strength 
  static constexpr double beta = 1e-12;            // Transition rate param
  static constexpr double Rs = std::pow(8 * M_PI, 1. / 3.) * vw / beta;   
  static constexpr double dtau = 10.0 * Rs;         // PT duration
  static constexpr double TN = dflt_universe::Ts;     // Nucleation temperature
  static constexpr double cpsq = 1.0 / 3.0;      // speed of sound squared (symmetric phase)
  static constexpr double cmsq = cpsq;           // speed of sound squared (broken phase
  static constexpr const char* nuc_type = "exp"; // bubble nucleation type
};


/* 
units:
[vw] = dimensionless (0 < vw < 1)
[alN] = dimensionless (alN > 0)
[beta] = GeV
[dtau] = 1/GeV
*/

  /**
 * @class PTParams
 * @brief Class representing the parameters of the phase transition (PT).
 * 
 * This class encapsulates the parameters needed to describe the phase transition dynamics,
 * including speeds of sound, wall velocity, strength of the transition, and bubble nucleation type.
 */
 class PTParams { // will probably need to update this later
  public:
    // ctors
    PTParams();
    PTParams(double vw, double alN);
    PTParams(double vw, double alN, double beta, double dtau, double TN, double cpsq, double cmsq, const char* nuc_type, const Universe& un);
    PTParams(double vw, double alN, const std::string& eos);
    PTParams(double vw, double alN, double beta, double dtau, double TN, double cpsq, double cmsq, const char* nuc_type, const Universe& un, const std::string& eos);

    Universe un() const { return universe_; } // universe parameters

    double cpsq() const { return cpsq_; } // speed of sound squared (symmetric phase)
    double cmsq() const { return cmsq_; } // speed of sound squared (broken phase)
    double vw() const { return vw_; } // wall velocity
    double alN() const { return alN_; } // strength parameter at nuc temp (alN_N)
    double beta() const { return beta_; } // inverse PT duration
    double Rs() const { return Rs_; } // characteristic length scale R_*
    double tau_s() const { return tau_s_; } // start time of PT
    double tau_fin() const { return tau_fin_; } // end time of PT
    double dtau() const { return dtau_; } // PT duration
    double TN() const { return TN_; } // nucleation temperature
    double wNeN_rat() const { return wNeN_rat_; } // ratio of nucleation enthalpy and energy density

    const char* nuc_type() const { return nuc_type_; } // bubble nucleation type
    const std::string eos_model() const { return eos_model_; } // equation of state model (bag or Veff)

    // for Veff eos
    std::vector<double> veff_TTN_vals() const { return veff_TTN_vals_; }
    std::vector<double> veff_ps_vals() const { return veff_ps_vals_; }
    std::vector<double> veff_pb_vals() const { return veff_pb_vals_; }
    std::vector<double> veff_es_vals() const { return veff_es_vals_; }
    std::vector<double> veff_eb_vals() const { return veff_eb_vals_; }

    double ps_val(double TTN) const { return alglib::spline1dcalc(veff_ps_interp_, TTN); }
    double pb_val(double TTN) const { return alglib::spline1dcalc(veff_pb_interp_, TTN); }
    double es_val(double TTN) const { return alglib::spline1dcalc(veff_es_interp_, TTN); }
    double eb_val(double TTN) const { return alglib::spline1dcalc(veff_eb_interp_, TTN); }
    double ws_val(double TTN) const { return alglib::spline1dcalc(veff_ws_interp_, TTN); }
    double wb_val(double TTN) const { return alglib::spline1dcalc(veff_wb_interp_, TTN); }

    double TTN_min() const { return veff_TTN_vals_.front(); }
    double TTN_max() const { return veff_TTN_vals_.back(); }

    double pN() const { return pN_; }
    double eN() const { return eN_; }
    double wN() const { return wN_; }

    friend std::ostream& operator<<(std::ostream& os, const PTParams& p);
    void print() const;

    #ifdef ENABLE_MATPLOTLIB
    void plot_thermo(const std::string& filename) const; // Plots e(T), p(T), w(T)
    #endif
  
  private:
      const Universe universe_;
      std::string eos_model_; // equation of state model (bag or Veff)
      const char *nuc_type_;
      double vw_, alN_, beta_, Rs_, tau_s_, tau_fin_, dtau_, TN_, wNeN_rat_, cpsq_, cmsq_;

      // for Veff eos
      std::vector<double> veff_TTN_vals_, veff_ps_vals_, veff_pb_vals_, veff_es_vals_, veff_eb_vals_, veff_ws_vals_, veff_wb_vals_;
      alglib::spline1dinterpolant veff_ps_interp_, veff_pb_interp_, veff_es_interp_, veff_eb_interp_, veff_ws_interp_, veff_wb_interp_;
      double pN_, eN_, wN_;


      bool is_valid_model(const char* model, const char* allowed_models[], const int n) const;
      bool is_valid_csq(double csq) const;
    };

} // namespace PhaseTransition

#endif // INCLUDE_PHASETRANSITION_HPP_H
