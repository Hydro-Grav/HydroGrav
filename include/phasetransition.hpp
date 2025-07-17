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
  static constexpr double T0 = 2.725;
  static constexpr double Ts = 100.0;
  static constexpr double H0 = 67.8;
  static constexpr double Hs = 1.0;
  static constexpr double g0 = 3.91;
  static constexpr double gs = 100.0;
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
    Universe(double T0, double Ts, double H0, double Hs, double g0, double gs);
    
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
  static constexpr double alpha = 0.1;           // PT strength 
  static constexpr double beta = 1.0;            // Transition rate param
  static constexpr double dtau = 10.0;            // PT duration
  static constexpr const char* model = "bag";    // equation of state
  static constexpr const char* nuc_type = "exp"; // bubble nucleation type
  static constexpr double wN = 1.71;             // enthalpy at nuc temp
};

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
    PTParams(double vw, double alpha, double beta, double dtau, double wN, const char* model, const char* nuc_type, const Universe& un);

    double cpsq() const { return cpsq_; } // speed of sound squared (symmetric phase)
    double cmsq() const { return cmsq_; } // speed of sound squared (broken phase)
    double csq() const { return cmsq_; } // keep this? - applies to bag model only
    double vw() const { return vw_; } // wall velocity
    double vcj() const { return vcj_; } // Chapman-Jouget speed
    double alphaN() const { return alpha_; } // strength parameter at nuc temp (alpha_N)
    double beta() const { return beta_; } // inverse PT duration
    double Rs() const { return Rs_; } // characteristic length scale R_*
    double tau_s() const { return tau_s_; } // start time of PT
    double tau_fin() const { return tau_fin_; } // end time of PT
    double dtau() const { return dtau_; } // PT duration
    double wN() const { return wN_; } // enthalpy at nuc temp

    const char* model() const { return model_; } // equation of state model
    const char* nuc_type() const { return nuc_type_; } // bubble nucleation type
    std::string wall_type() const { return wall_type_; } // deflagration, hybrid, or detonation

    // print params
    void print() const;
    friend std::ostream& operator<<(std::ostream& os, const PTParams& p);
  
  private:
      const Universe universe_;
      double vw_, alpha_, beta_, Rs_, tau_s_, tau_fin_, dtau_, cpsq_, cmsq_, vcj_, wN_;
      const char *model_, *nuc_type_;
      std::string wall_type_;

      bool is_valid_model(const char* model, const char* allowed_models[], const int n);
};

} // namespace PhaseTransition

#endif // INCLUDE_PHASETRANSITION_HPP_H
