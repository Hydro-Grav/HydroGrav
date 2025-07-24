// PhaseTransition.hpp
#ifndef INCLUDE_PHASETRANSITION_HPP_H
#define INCLUDE_PHASETRANSITION_HPP_H

#include <string>
#include <iostream>
#include <vector>

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

// move universe class somewhere else?
struct dflt_universe {
  static constexpr double T0 = 2.41e-13; // GeV
  static constexpr double Ts = 100.0; // GeV
  static constexpr double H0 = 1.45e-42; // GeV
  static constexpr double Hs = 1.41e-14; // GeV
  static constexpr double g0 = 3.91;
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
    Universe()
      : Universe(dflt_universe::T0, dflt_universe::Ts, dflt_universe::H0, dflt_universe::Hs, dflt_universe::g0, dflt_universe::gs) {}
    Universe(double Ts, double Hs, double gs)
      : T0_(dflt_universe::T0), Ts_(Ts), H0_(dflt_universe::H0), Hs_(Hs), g0_(dflt_universe::g0), gs_(gs) {}
    Universe(double T0, double Ts, double H0, double Hs, double g0, double gs)
      : T0_(T0), Ts_(Ts), H0_(H0), Hs_(Hs), g0_(g0), gs_(gs) {}

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
    const double T0_, Ts_, H0_, Hs_, g0_, gs_;
};

const Universe& default_universe();

struct dflt_PTParams {
  static constexpr double vw = 0.8;              // Wall velocity
  static constexpr double alN = 0.1;           // PT strength 
  static constexpr double beta = 1e-12;            // Transition rate param
  static constexpr double dt = 1e12;            // PT duration
  static constexpr double wNeN_rat = 1.0 + 1./3.;       // wN/eN = 1 + pN/eN = 1 + 1/3 for bag model
  static constexpr const char* nuc_type = "exp"; // bubble nucleation type
};


/* 
units:
[vw] = dimensionless (0 < vw < 1)
[alN] = dimensionless (alN > 0)
[beta] = GeV
[dt] = 1/GeV
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
    PTParams(double vw, double alN, double beta, double dt, double wNeN_rat, const char* nuc_type, const Universe& un);

    Universe un() const { return universe_; } // universe parameters

    double cpsq() const { return cpsq_; } // speed of sound squared (symmetric phase)
    double cmsq() const { return cmsq_; } // speed of sound squared (broken phase)
    double vw() const { return vw_; } // wall velocity
    double alN() const { return alN_; } // strength parameter at nuc temp (alN_N)
    double beta() const { return beta_; } // inverse PT duration
    double Rs() const { return Rs_; } // characteristic length scale R_*
    double tau_s() const { return tau_s_; } // start time of PT
    double tau_fin() const { return tau_fin_; } // end time of PT
    double dt() const { return dt_; } // PT duration
    double wNeN_rat() const { return wNeN_rat_; } // ratio of enthalpy to energy density (wN/eN)

    const std::string eos_model() const { return eos_model_; } // equation of state model (bag or Veff)
    const char* nuc_type() const { return nuc_type_; } // bubble nucleation type

    // print params
    void print() const;
    friend std::ostream& operator<<(std::ostream& os, const PTParams& p);
  
  private:
      const Universe universe_;
      double vw_, alN_, beta_, Rs_, tau_s_, tau_fin_, dt_, wNeN_rat_, cpsq_, cmsq_;
      std::string eos_model_; // equation of state model (bag or Veff)
      const char *nuc_type_;

      void check_valid_params() const;
      bool is_valid_model(const char* model, const char* allowed_models[], const int n) const;
      bool is_valid_csq(double csq) const;
    };

} // namespace PhaseTransition

#endif // INCLUDE_PHASETRANSITION_HPP_H
