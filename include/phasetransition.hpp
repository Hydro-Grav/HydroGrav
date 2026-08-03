/**
 * @file phasetransition.hpp
 * @brief Classes and helpers describing phase transition parameters,
 *        equations of state, and universe properties.
 *
 * This header provides `Universe` and `PTParams` classes alongside
 * default-parameter structs and EOS data structures used throughout the
 * gravitational-wave calculations.
 */
#ifndef INCLUDE_PHASETRANSITION_HPP_H
#define INCLUDE_PHASETRANSITION_HPP_H

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <optional>

#include "ap.h"
#include "interpolation.h"

#include "physics.hpp"

/**
 * @namespace PhaseTransition
 * @brief Contains classes that store phase transition parameters
 */
namespace PhaseTransition {

// DO NOT CHANGE DEFAULT VALS
/**
 * @brief Default numerical values for the universe parameters in natural units (\f$\hbar=c=k_B=1\f$).
 */
struct dflt_universe { // in units hbar = c = kB = 1
  // values today
  static constexpr double T0 = 2.34914e-13; // 2.725 K / (1.16e+13 K/GeV) = 2.349e-13 GeV
  static constexpr double H0 = 1.44328e-42; // 67.8 km/s/Mpc = 2.09502e21 s^-1 = 1.44328e-42 GeV
  static constexpr double g0 = 3.91;

  // values at start of PT
  static constexpr double Ts = 52.9772; // GeV
  static constexpr double gs = 106.75;
  static constexpr double Hs = 4.34679e-15;
  // static constexpr double Hs = std::pow(4.0 * std::pow(M_PI, 3) * gs * std::pow(Ts, 4) / (45.0 * mP * mP), 0.5);
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

/**
 * @brief Retrieve a reference to default universe object.
 */
const Universe& default_universe();

/**
 * @brief Approximate mean bubble separation \f$R_*\f$ from wall velocity and
 *        transition rate parameter.
 *
 * @param vw   Wall velocity.
 * @param beta Inverse duration parameter of the phase transition.
 * @return Approximate value of \f$R_*\f$.
 */
constexpr double Rs_approx(double vw, double beta) { return std::pow(8 * M_PI, 1. / 3.) * vw / beta; };

// DO NOT CHANGE DEFAULT VALS
/**
 * @brief Default values for the phase transition parameters.
 *
 * These values are used when the user does not supply explicit PT
 * parameters. They correspond to a benchmark point for xSM.
 */
struct dflt_PTParams {
  static constexpr double vw = 0.675122; // Wall velocity
  static constexpr double alN_bag = 0.0972393; // PT strength 
  static constexpr double alN_munu = 0.103742; // PT strength (mu nu)
  static constexpr double betaHs = 1231.05; // Transition rate param / Hs
  static constexpr double beta = betaHs * dflt_universe::Hs; // Transition rate parameter
  static constexpr double Rs = Rs_approx(vw, beta); // mean bubble separation
  static constexpr double TN = dflt_universe::Ts; // Nucleation temperature
  static constexpr double cpsq = 0.56705 * 0.56705; // speed of sound squared (symmetric phase)
  static constexpr double cmsq = 0.539046 * 0.539046; // speed of sound squared (broken phase)
  static inline const std::string nuc_type = "exp"; // bubble nucleation type
};

/**
 * @brief Lifetime distribution function of bubbles.
 *
 * Encapsulates the lifetime distribution function nu(T_tilde), where T_tilde
 * is a dimensionless time variable.
 */
struct LifetimeDistribution {

  enum class ZeroHandling {
    Throw, // throws an exception if nu <= 0
    Clamp, // clamps nu to a small positive value
  };

  std::vector<double> Ttilde_values;
  std::vector<double> nu_values;
  std::vector<double> log_nu_values;

  LifetimeDistribution(
    const std::vector<double>& Ttilde_vals,
    const std::vector<double>& nu_vals,
    ZeroHandling zero_handling = ZeroHandling::Clamp)
    : Ttilde_values(Ttilde_vals)
    , nu_values(nu_vals)
    , log_nu_values(nu_vals.size())
  {

    if (Ttilde_values.size() != nu_values.size()) {
      throw std::invalid_argument(
        "LifetimeDistribution: Ttilde and nu vectors must have the same size, got " + std::to_string(Ttilde_values.size()) + " vs "+ std::to_string(nu_values.size()));
    }
    if (Ttilde_values.size() < 2) {
      throw std::invalid_argument(
        "LifetimeDistribution: at least 2 data points are required for spline construction");
    }

    for (size_t i = 0; i < nu_values.size(); i++) {
      if (nu_values[i] <= 0.0) {
        if (zero_handling == ZeroHandling::Throw) {
          throw std::invalid_argument(
            "LifetimeDistribution: nu_values[" + std::to_string(i) + "] = "
            + std::to_string(nu_values[i]) + " is <= 0. "
            "Use ZeroHandling::Clamp if you want non-positive values replaced automatically.");
        }
        log_nu_values[i] = std::log(std::numeric_limits<double>::denorm_min());
      } else {
        log_nu_values[i] = std::log(nu_values[i]);
      }
    }

    alglib::real_1d_array Ttilde_arr, log_nu_arr;
    Ttilde_arr.setcontent(Ttilde_values.size(), Ttilde_values.data());
    log_nu_arr.setcontent(log_nu_values.size(), log_nu_values.data());

    try {
      alglib::spline1dbuildcubic(Ttilde_arr, log_nu_arr, nu_spline_);
    } catch (const alglib::ap_error& e) {
      throw std::runtime_error(
        std::string("LifetimeDistribution: error building spline: ") + e.msg);
    } catch (...) {
      throw std::runtime_error(
        "LifetimeDistribution: unknown error building spline");
    }
  }

  double operator()(double Ttilde) const {
    double T_clamped = std::clamp(Ttilde, Ttilde_values.front(), Ttilde_values.back());
    return std::exp(alglib::spline1dcalc(nu_spline_, T_clamped));
  }

  bool inDomain(double Ttilde) const {
    return Ttilde >= Ttilde_values.front() && Ttilde <= Ttilde_values.back();
  }

private:
  alglib::spline1dinterpolant nu_spline_;
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
    PTParams(double vw, double alN, double TN, double beta, double Rs, const std::string nuc_type, const Universe& un); // ctor
    virtual ~PTParams() = default; // dtor

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
    const std::string nuc_type() const { return nuc_type_; } // bubble nucleation type

    // friend std::ostream& operator<<(std::ostream& os, const PTParams& p);
    // void print() const;
    
    virtual double cpsq() const = 0; // speed of sound squared (symmetric phase)
    virtual double cmsq() const = 0; // speed of sound squared (broken phase)

    void set_lifetime_distribution(const LifetimeDistribution& lt_dist) { lt_dist_ = lt_dist; }
    std::optional<LifetimeDistribution> get_lifetime_distribution() const { return lt_dist_; }

  protected:
    const Universe un_;
    double vw_, alN_, TN_, wNeN_rat_, beta_, Rs_, tau_s_;
    std::string nuc_type_;
    std::optional<LifetimeDistribution> lt_dist_;

    virtual void print() const;
  
  private:
    bool is_valid_model(const std::string& model, const std::vector<std::string>& allowed_models) const;
};

/**
 * @class PTParams_Bag
 * @brief Phase transition parameters for the bag and mu-nu equations of state.
 */
class PTParams_Bag : public PTParams {
  public:
    // ctors
    /**
     * @brief Construct bag-model parameters with default sound speed.
     *
     * `cpsq` and `cmsq` default to 1/3 when not provided.
     */
    PTParams_Bag(double vw, double alN, double TN, double beta, double Rs, const std::string nuc_type, const Universe& un);
    /**
     * @brief Full constructor specifying speeds of sound squared.
     *
     * @param cpsq Sound speed squared in symmetric phase.
     * @param cmsq Sound speed squared in broken phase.
     */
    PTParams_Bag(double vw, double alN, double TN, double beta, double Rs, const std::string nuc_type, const Universe& un, double cpsq, double cmsq);

    ModelType eos() const override { return ModelType::Bag; }

    double cpsq() const override { return cpsq_; }
    double cmsq() const override { return cmsq_; }

    void print() const override;

  private:
    const double cpsq_, cmsq_;
    bool is_valid_csq(double csq) const;
};

/**
 * @class EquationOfState
 * @brief Container for thermodynamic data used by the Veff model.
 *
 * Holds temperature and pressure/energy density arrays along with
 * factory methods for reading from a file. The data are used to build
 * spline interpolants for various thermodynamic quantities.
 */
struct EquationOfState 
{
  std::vector<double> T_vals;
  std::vector<double> ps_vals;
  std::vector<double> pb_vals;
  std::vector<double> es_vals;
  std::vector<double> eb_vals;

  /**
   * @brief Construct an EquationOfState from tabulated data.
   *
   * @param T  Temperatures.
   * @param ps Pressures in symmetric phase.
   * @param pb Pressures in broken phase.
   * @param es Energy densities in symmetric phase.
   * @param eb Energy densities in broken phase.
   */
  EquationOfState(
    const std::vector<double>& T,
    const std::vector<double>& ps,
    const std::vector<double>& pb,
    const std::vector<double>& es,
    const std::vector<double>& eb);

  /**
   * @brief Load EOS data from a file and construct an EquationOfState.
   *
   * The file format must match that produced by the EOS preprocessing
   * utilities.
   *
   * @param filename Path to data file.
   * @return Constructed EquationOfState object.
   */
  static EquationOfState from_file(const std::string& filename);

  /**
   * @brief Check whether the EOS data vectors are non-empty and consistent.
   *
   * @return `true` if the object contains valid data.
   */
  bool is_valid() const;
  size_t size() const { return T_vals.size(); }

  void write(const std::string& filename="thermo.csv") const;

private:
  void validate() const;
};

/**
 * @class PTParams_Veff
 * @brief Phase transition parameters using a generic effective potential
 *        equation of state.
 */
class PTParams_Veff : public PTParams {
public:
  // These use the new eos_data
  /**
   * @brief Construct Veff-model parameters using EOS data.
   *
   * Beta and Rs are computed internally from the EOS.
   */
  PTParams_Veff(double vw, double alN, double TN, const EquationOfState& eos_data);
  /**
   * @brief Full constructor feeding all GW parameters explicitly.
   */
  PTParams_Veff(double vw, double alN, double TN, double beta, double Rs, const std::string nuc_type, const Universe& un, const EquationOfState& eos_data);

  // Backward compatibile
  PTParams_Veff(double vw, double alN, double TN, const std::string& veff_eos_filename);
  PTParams_Veff(double vw, double alN, double TN, double beta, double Rs, const std::string nuc_type, const Universe& un, const std::string& veff_eos_filename);

  ModelType eos() const override { return ModelType::Veff; }

  std::vector<double> veff_TTN_vals() const { return veff_TTN_vals_; }
  std::vector<double> veff_ps_vals() const { return veff_ps_vals_; }
  std::vector<double> veff_pb_vals() const { return veff_pb_vals_; }
  std::vector<double> veff_es_vals() const { return veff_es_vals_; }
  std::vector<double> veff_eb_vals() const { return veff_eb_vals_; }

  /** @brief Pressure in symmetric phase at given T/TN. */
  double ps_val(double TTN) const { return alglib::spline1dcalc(veff_ps_interp_, TTN); }
  /** @brief Pressure in broken phase at given T/TN. */
  double pb_val(double TTN) const { return alglib::spline1dcalc(veff_pb_interp_, TTN); }
  /** @brief Energy density in symmetric phase at given T/TN. */
  double es_val(double TTN) const { return alglib::spline1dcalc(veff_es_interp_, TTN); }
  /** @brief Energy density in broken phase at given T/TN. */
  double eb_val(double TTN) const { return alglib::spline1dcalc(veff_eb_interp_, TTN); }
  /** @brief Enthalpy density in symmetric phase. */
  double ws_val(double TTN) const { return alglib::spline1dcalc(veff_ws_interp_, TTN); }
  /** @brief Enthalpy density in broken phase. */
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
