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
- remove alN from PTParams base class and move to PTParams_Bag (need to change how get_mode() works in FluidProfile first)
- TN only used for mu nu and Veff (not bag) - write ctor without TN for bag?
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
  static constexpr double Hs = std::pow(4.0 * std::pow(M_PI, 3) * gs * std::pow(Ts, 4) / (45.0 * mP * mP), 0.5);
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
    Universe(double Ts, double gs, double Hs);
    Universe(double T0, double Ts, double g0, double gs, double H0, double Hs);

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

const Universe& default_universe();

constexpr double Rs_approx(double vw, double beta) { return std::pow(8 * M_PI, 1. / 3.) * vw / beta; };

// DO NOT CHANGE DEFAULT VALS
struct dflt_PTParams {
  static constexpr double vw = 0.8;              // Wall velocity
  static constexpr double alN = 0.1;           // PT strength 
  static constexpr double beta = 1e-12;            // Transition rate param
  static constexpr double Rs = Rs_approx(vw, beta);
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
*/

/**
 * @class PTParams
 * @brief Class representing the parameters of the phase transition (PT).
 * 
 * This class encapsulates the parameters needed to describe the phase transition dynamics,
 * including speeds of sound, wall velocity, strength of the transition, and bubble nucleation type.
 */
class PTParams {
  public:
    // ctor
    PTParams(double vw, double alN, double TN, double beta, double Rs, const char* nuc_type, const Universe& un);

    enum class ModelType { Bag, Veff }; // equation of state model
    virtual ModelType eos() const = 0;
    inline std::string eos_to_string() const {
        switch (eos()) {
            case ModelType::Bag:  return "Bag";
            case ModelType::Veff: return "Veff";
            default:              return "Unknown";
        }
    }

    Universe un() const { return un_; } // universe parameters

    // Fluid parameters
    double vw() const { return vw_; } // wall velocity
    double alN() const { return alN_; } // strength parameter at nuc temp (alN_N)
    double TN() const { return TN_; } // nucleation temperature
    double wNeN_rat() const { return wNeN_rat_; } // ratio of nucleation enthalpy and energy density

    // GW parameters
    double beta() const { return beta_; } // inverse PT duration
    double betaHs() const { return beta_ / un_.Hs(); } // beta/Hs
    double Rs() const { return Rs_; } // characteristic length scale R_*
    double tau_s() const { return tau_s_; } // start time of PT
    double tau_fin() const { return tau_fin_; } // end time of PT
    const char* nuc_type() const { return nuc_type_; } // bubble nucleation type

    // friend std::ostream& operator<<(std::ostream& os, const PTParams& p);
    // void print() const;
    
    virtual double cpsq() const = 0; // speed of sound squared (symmetric phase)
    virtual double cmsq() const = 0; // speed of sound squared (broken phase)

  protected:
    const Universe un_;
    double vw_, alN_, TN_, wNeN_rat_, beta_, Rs_, tau_s_, tau_fin_;
    const char *nuc_type_;

    virtual void print() const;
  
  private:
    bool is_valid_model(const char* model, const char* allowed_models[], const int n) const;
};

class PTParams_Bag : public PTParams {
  public:
    // ctors
    PTParams_Bag(double vw, double alN, double TN, double beta, double Rs, const char* nuc_type, const Universe& un);
    PTParams_Bag(double vw, double alN, double TN, double beta, double Rs, const char* nuc_type, const Universe& un, double cpsq, double cmsq);

    ModelType eos() const override { return ModelType::Bag; }

    double cpsq() const override { return cpsq_; }
    double cmsq() const override { return cmsq_; }

    void print() const override;

  private:
    const double cpsq_, cmsq_;
    bool is_valid_csq(double csq) const;
};

struct EquationOfState 
{
  std::vector<double> T_vals;
  std::vector<double> ps_vals;
  std::vector<double> pb_vals;
  std::vector<double> es_vals;
  std::vector<double> eb_vals;

  EquationOfState(
    const std::vector<double>& T,
    const std::vector<double>& ps,
    const std::vector<double>& pb,
    const std::vector<double>& es,
    const std::vector<double>& eb);

  static EquationOfState from_file(const std::string& filename);

  bool is_valid() const;
  size_t size() const { return T_vals.size(); }

private:
  void validate() const;
};

class PTParams_Veff : public PTParams {
public:
  // These use the new eos_data
  PTParams_Veff(double vw, double alN, double TN, const EquationOfState& eos_data);
  PTParams_Veff(double vw, double alN, double TN, double beta, double Rs, const char* nuc_type, const Universe& un, const EquationOfState& eos_data);

  // Backward compatibile
  PTParams_Veff(double vw, double alN, double TN, const std::string& veff_eos_filename);
  PTParams_Veff(double vw, double alN, double TN, double beta, double Rs, const char* nuc_type, const Universe& un, const std::string& veff_eos_filename);

  ModelType eos() const override { return ModelType::Veff; }

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

  double csq_s(double TTN) const { return alglib::spline1dcalc(cpsq_fit_, TTN); } // WARNING: spline can go out of bounds
  double csq_b(double TTN) const { return alglib::spline1dcalc(cmsq_fit_, TTN); }

  // estimate cpsq, cmsq (needed to determine hydrodynamic mode)
  // Note: not perfect, since cpsq = csq_s(Tp/TN), cmsq = csq_b(Tm/TN) and Tp,Tm != TN in general
  double cpsq() const override { return csq_s(1.0); }
  double cmsq() const override { return csq_b(1.0); }

  double pN() const { return pN_; }
  double eN() const { return eN_; }
  double wN() const { return wN_; }

  #ifdef ENABLE_MATPLOTLIB
  void plot_thermo(const std::string& filename="thermo.png") const; // Plots e(T), p(T), w(T)
  void plot_thermo2(const std::string& filename="thermo.png") const;
  void plot_csq(const std::string& filename="csq_veff.png") const; // Plots cs^2(T)
  #endif

  void print() const override;

private:
  std::vector<double> veff_TTN_vals_, veff_ps_vals_, veff_pb_vals_, veff_es_vals_, veff_eb_vals_, veff_ws_vals_, veff_wb_vals_;
  std::vector<double> cpsq_vals_, cmsq_vals_; // cs^2(T) values
  alglib::spline1dinterpolant veff_ps_interp_, veff_pb_interp_, veff_es_interp_, veff_eb_interp_, veff_ws_interp_, veff_wb_interp_;
  alglib::spline1dinterpolant cpsq_fit_, cmsq_fit_; // cs^2(T)
  double pN_, eN_, wN_;

  void initialize_from_eos_data(const EquationOfState& eos_data);
};

} // namespace PhaseTransition

#endif // INCLUDE_PHASETRANSITION_HPP_H
