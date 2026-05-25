// PhaseTransition.cpp
#include <string>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <unordered_set>
#include <cassert>
#include <fstream>

#include "ap.h"
#include "interpolation.h"

#ifdef ENABLE_MATPLOTLIB
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

#include "phasetransition.hpp"
#include "maths_ops.hpp"

/*
TO DO:
- change taus to 1/Hs_conformal (currently using 1/Hs since otherwise program breaks)
- remove alN from PTParams base class and move to PTParams_Bag (need to change how get_mode() works in FluidProfile first)
- change all spline vectors to fixed length arrays (faster)
*/

namespace PhaseTransition {

Universe::Universe()
    : Universe(dflt_universe::T0, dflt_universe::Ts, dflt_universe::g0, dflt_universe::gs, dflt_universe::H0, dflt_universe::Hs) {}

Universe::Universe(double Ts, double gs, double Hs)
    : Universe(dflt_universe::T0, Ts, dflt_universe::g0, gs, dflt_universe::H0, Hs) {}

Universe::Universe(double T0, double Ts, double g0, double gs, double H0, double Hs)
    : T0_(T0), 
      Ts_(Ts), 
      g0_(g0), 
      gs_(gs), 
      H0_(H0), 
      Hs_(Hs) {}

std::ostream& operator<<(std::ostream& os, const Universe& un) {
    os << "************** Universe parameters **************\n"
       << std::left
       << std::setw(20) << " " << std::setw(15) << "Today" << "Start of PT\n"
       << std::setw(20) << " " << std::setw(15) << "-----" << "-----------\n"
       << std::setw(20) << "Temperature:" << "T0=" << std::setw(12) << un.T0() << "Ts=" << un.Ts() << "\n"
       << std::setw(20) << "Hubble constant:" << "H0=" << std::setw(12) << un.H0() << "Hs=" << un.Hs() << "\n"
       << std::setw(20) << "Number of DoF:" << "g0=" << std::setw(12) << un.g0() << "gs=" << un.gs() << "\n"
       << "*************************************************\n";
       
    return os;
}

// some way to combine this with PTParams print()?
void Universe::print() const {
    std::cout << *this;
}

const Universe& default_universe() {
    static Universe u;
    return u;
}

/************************************ PTParams ************************************/
PTParams::PTParams(double vw, double alN, double TN, double beta, double Rs, const std::string nuc_type, const Universe& un)
    : un_(un),
      vw_(vw),
      alN_(alN),
      TN_(TN),
      wNeN_rat_(std::numeric_limits<double>::quiet_NaN()),
      beta_(beta),
      Rs_(Rs),
      tau_s_(1.0 / un.Hs()),
      nuc_type_(nuc_type) {

      /*
      Note on Hs_conformal:
      Usual definition of conformal Hubble rate is Hs_conformal = as*Hs, but here we include a redshift factor 1/a0 (i.e. so this is Hs_conformal as measured today)
      tau_s should be defined as 1.0 / Hs_conformal (don't understand why!), but this breaks calculation
      since tau_s is too large -> Si and Ci integrals are const for such large values so dlt=0 -> OmegaGW=0

      const auto asa0_rat = std::pow(un_.g0() / un_.gs(), 1./3.) * un_.T0() / un_.Ts(); // a_* / a_0
      const auto Hs_conformal = un_.Hs() * asa0_rat; // conformal Hubble rate at PT as measured today
      tau_s_ = 1.0 / Hs_conformal;
      */

      // check valid vw
      if (vw_ < 0.0 ) {
        std::cerr << "Warning: vw < 0. Taking |vw| as input instead.";
        vw_ = std::abs(vw);
      } else if (vw == 0.0 || vw >= 1.0) {
        throw std::invalid_argument("Unphysical wall velocity passed into PTParams. Must have 0 < vw < 1.");
      }

      // check valid alN
      if (alN_ <= 0.0) {
        throw std::invalid_argument("Unphysical strength parameter passed into PTParams. Must have alN > 0.");
      }

      // check valid TN
      if (TN_ <= 0.0) {
        throw std::invalid_argument("Unphysical nucleation temperature passed into PTParams. Must have TN > 0.");
      }

      // Check valid beta
      if (beta_ <= 0.0) {
        throw std::invalid_argument("Unphysical transition rate parameter passed into PTParams. Must have beta > 0.");
      }

      // check valid bubble nucleation type
      const std::vector<std::string> allowed_nuc = {"exp", "sim"};
      if (!is_valid_model(nuc_type, allowed_nuc)) {
          std::cerr << "Warning: Invalid model '" << nuc_type << "' for bubble nucleation. Using default nucleation type (" << dflt_PTParams::nuc_type << ")\n";
          nuc_type_ = dflt_PTParams::nuc_type;
      }
    }

// Protected:
void PTParams::print() const {
  std::cout << "********** Phase Transition parameters **********\n"
            << std::left
            //  << std::setw(35) << "Equation of state:" << params.model_ << "\n"
            << std::setw(35) << "Nucleation type:" << nuc_type_ << "\n"
            << std::setw(35) << "Wall velocity:" << "vw=" << vw_ << "\n"
            << std::setw(35) << "PT strength parameter:" << "alN=" << alN_ << "\n"
            << std::setw(35) << "Transition rate parameter:" << "beta=" << beta_ << "\n"
            << std::setw(35) << "Mean bubble separation:" << "Rs=" << Rs_ << "\n";
}

// Private:
bool PTParams::is_valid_model(const std::string& model, const std::vector<std::string>& allowed_models) const {
    return std::find(allowed_models.begin(), allowed_models.end(), model) != allowed_models.end();
}

/********************************** PTParams_Bag **********************************/
PTParams_Bag::PTParams_Bag(double vw, double alN, double TN, double beta, double Rs, const std::string nuc_type, const Universe& un) // Bag model
    : PTParams_Bag(vw, alN, TN, beta, Rs, nuc_type, un, 1.0 / 3.0, 1.0 / 3.0) {}

PTParams_Bag::PTParams_Bag(double vw, double alN, double TN, double beta, double Rs, const std::string nuc_type, const Universe& un, double cpsq, double cmsq) // full ctor
    : PTParams(vw, alN, TN, beta, Rs, nuc_type, un),
      cpsq_(cpsq),
      cmsq_(cmsq) {

      std::cout << "Storing phase transition parameters. Note that alN definition differs between bag and mu-nu models!\n";

      // check valid speed of sound
      if (!is_valid_csq(cpsq_)) {
        throw std::invalid_argument("Unphysical speed of sound in symmetric phase passed into PTParams_Bag. Must have 0 < cpsq < 1.");
      }
      if (!is_valid_csq(cmsq_)) {
        throw std::invalid_argument("Unphysical speed of sound in broken phase passed into PTParams_Bag. Must have 0 < cmsq < 1.");
      }

      // this is fine for mu nu since csq assumed to be constant everywhere in symmetric phase (csq=cpsq in front of shock too)
      wNeN_rat_ = 1.0 + cpsq_; // wN/eN = 1 + pN/eN = 1 + 1/3 for bag model
    }

// Public:
void PTParams_Bag::print() const {
  PTParams::print();
  std::cout << std::setw(35) << "Speed of sound (symmetric phase):" << "cpsq=" << cpsq_ << "\n"
            << std::setw(35) << "Speed of sound (broken phase):" << "cmsq=" << cmsq_ << "\n"
            << "*************************************************\n";
}

// Private:
bool PTParams_Bag::is_valid_csq(double csq) const {
  return (csq > 0.0 && csq < 1.0);
}

EquationOfState::EquationOfState(
  const std::vector<double>& T,
  const std::vector<double>& ps,
  const std::vector<double>& pb,
  const std::vector<double>& es,
  const std::vector<double>& eb
) : T_vals(T), ps_vals(ps), pb_vals(pb), es_vals(es), eb_vals(eb) 
{
  validate();
}

EquationOfState EquationOfState::from_file(const std::string& filename) 
{
  if (filename.empty()) 
  {
    throw std::invalid_argument("Equation of state filename cannot be empty");
  }

  std::cout << "Reading equation of state from file: " << filename << "\n"
            << "Note: File must be formatted as T, ps, pb, es, eb (comma separated) with header line\n";
  
  std::ifstream file(filename);
  if (!file) 
  {
    throw std::runtime_error("Could not open file " + filename);
  }

  std::string line;
  std::getline(file, line); // Skip header

  std::vector<double> T_vals, ps_vals, pb_vals, es_vals, eb_vals;

  while (std::getline(file, line)) 
  {
    std::istringstream ss(line);
    std::array<double, 5> values;
    std::string token;

    for (auto& val : values) 
    {
      if (!std::getline(ss, token, ',')) 
      {
        throw std::runtime_error("Malformed line in " + filename + ": " + line);
      }
      val = std::stod(token);
    }

    T_vals.push_back(values[0]);
    ps_vals.push_back(values[1]);
    pb_vals.push_back(values[2]);
    es_vals.push_back(values[3]);
    eb_vals.push_back(values[4]);
  }

  std::cout << "Equation of state read successfully! (" << T_vals.size() << " data points)\n";
  
  return EquationOfState(T_vals, ps_vals, pb_vals, es_vals, eb_vals);
}

void EquationOfState::validate() const 
{
  const auto n = T_vals.size();
  
  if (n < 2) 
  {
    throw std::invalid_argument("Equation of state data must contain at least 2 data points");
  }

  if (ps_vals.size() != n || pb_vals.size() != n || es_vals.size() != n || eb_vals.size() != n) 
  {
    throw std::invalid_argument("All equation of state vectors must have the same size");
  }

  // Check temperatures are monotonically increasing
  for (size_t i = 1; i < n; ++i) 
  {
    if (T_vals[i] <= T_vals[i-1]) 
    {
      throw std::invalid_argument("Temperature values must be strictly increasing");
    }
  }
}

bool EquationOfState::is_valid() const 
{
  try {
    validate();
    return true;
  } catch (...) {
    return false;
  }
}

void EquationOfState::write(const std::string& filename) const {
    std::cout << "Writing equation of state to disk... ";

    std::ofstream file(filename);
    file << "T,ps,pb,es,eb\n";

    for (size_t i = 0; i < T_vals.size(); ++i) {
        file << T_vals[i] << "," << ps_vals[i] << "," << pb_vals[i] << "," << es_vals[i] << "," << eb_vals[i] << "\n";
    }
    file.close();

    std::cout << "Saved to " << filename << "!\n";

    return;
}

/********************************** PTParams_Veff *********************************/
// NOTE: PTParams_Veff only uses alN to check if hydrodynamic mode agrees with Bag model!

// New primary constructor
PTParams_Veff::PTParams_Veff(double vw, double alN, double TN, const EquationOfState& eos_data)
    : PTParams_Veff(vw, alN, TN, dflt_PTParams::beta, dflt_PTParams::Rs, dflt_PTParams::nuc_type, default_universe(), eos_data) {}

PTParams_Veff::PTParams_Veff(double vw, double alN, double TN, double beta, double Rs, const std::string nuc_type, const Universe& un, const EquationOfState& eos_data)
    : PTParams(vw, alN, TN, beta, Rs, nuc_type, un) {
    
    std::cout << "Storing phase transition parameters. Note that PTParams_Veff must be given alN defined for mu-nu model to identify if hydrodynamic mode differs to simplified EoS!!\n";
    initialize_from_eos_data(eos_data);
}

// Backward compatibility constructors
PTParams_Veff::PTParams_Veff(double vw, double alN, double TN, const std::string& veff_eos_filename)
    : PTParams_Veff(vw, alN, TN, EquationOfState::from_file(veff_eos_filename)) {}

PTParams_Veff::PTParams_Veff(double vw, double alN, double TN, double beta, double Rs, const std::string nuc_type, const Universe& un, const std::string& veff_eos_filename)
    : PTParams_Veff(vw, alN, TN, beta, Rs, nuc_type, un, EquationOfState::from_file(veff_eos_filename)) {}

// Private initialization method
void PTParams_Veff::initialize_from_eos_data(const EquationOfState& eos_data) {
  const double TN_inv = 1.0 / TN_;

  // Convert to T/TN and store
  veff_TTN_vals_.reserve(eos_data.size());
  for (const auto& T : eos_data.T_vals) {
      veff_TTN_vals_.push_back(T * TN_inv);
  }

  veff_ps_vals_ = eos_data.ps_vals;
  veff_pb_vals_ = eos_data.pb_vals;
  veff_es_vals_ = eos_data.es_vals;
  veff_eb_vals_ = eos_data.eb_vals;

  // Calculate w = e + p
  const auto n = eos_data.size();
  veff_ws_vals_.reserve(n);
  veff_wb_vals_.reserve(n);
  for (size_t i = 0; i < n; ++i) {
      veff_ws_vals_.push_back(veff_es_vals_[i] + veff_ps_vals_[i]);
      veff_wb_vals_.push_back(veff_eb_vals_[i] + veff_pb_vals_[i]);
  }

  // Construct interpolating functions p(T/TN), e(T/TN) in s/b phases
  // Smooth spline needed here to remove numerical noise in eos
  alglib::real_1d_array veff_TTN_array, veff_ps_array, veff_es_array, veff_ws_array, veff_pb_array, veff_eb_array, veff_wb_array;
  alglib::spline1dfitreport rep;
  // const double smooth_fac = 1e-2;
  const double smooth_fac = 1e-2;
  const int basis_size = 50;

  veff_TTN_array.setcontent(n, veff_TTN_vals_.data());
  
  veff_ps_array.setcontent(n, veff_ps_vals_.data());
  veff_es_array.setcontent(n, veff_es_vals_.data());
  veff_ws_array.setcontent(n, veff_ws_vals_.data());
  alglib::spline1dfit(veff_TTN_array, veff_ps_array, basis_size, smooth_fac, veff_ps_interp_, rep);
  alglib::spline1dfit(veff_TTN_array, veff_es_array, basis_size, smooth_fac, veff_es_interp_, rep);
  alglib::spline1dfit(veff_TTN_array, veff_ws_array, basis_size, smooth_fac, veff_ws_interp_, rep);

  veff_pb_array.setcontent(n, veff_pb_vals_.data());
  veff_eb_array.setcontent(n, veff_eb_vals_.data());
  veff_wb_array.setcontent(n, veff_wb_vals_.data());
  alglib::spline1dfit(veff_TTN_array, veff_pb_array, basis_size, smooth_fac, veff_pb_interp_, rep);
  alglib::spline1dfit(veff_TTN_array, veff_eb_array, basis_size, smooth_fac, veff_eb_interp_, rep);
  alglib::spline1dfit(veff_TTN_array, veff_wb_array, basis_size, smooth_fac, veff_wb_interp_, rep);

  // Calculate thermodynamic quantities at nucleation (T/TN = 1)
  pN_ = alglib::spline1dcalc(veff_ps_interp_, 1.0);
  if (pN_ <= 0.0) {
    throw std::invalid_argument("Unphysical nucleation pressure. Must have pN > 0.");
  }

  eN_ = alglib::spline1dcalc(veff_es_interp_, 1.0);
  if (eN_ <= 0.0) {
    throw std::invalid_argument("Unphysical nucleation energy density. Must have eN > 0.");
  }

  wN_ = alglib::spline1dcalc(veff_ws_interp_, 1.0);
  if (wN_ <= 0.0) {
    throw std::invalid_argument("Unphysical nucleation enthalpy. Must have wN > 0.");
  }

  // Calculate sound speeds
  double s_unused, dps, des, dpb, deb, dps2_unused, des2_unused, dpb2_unused, deb2_unused;
  cpsq_vals_.reserve(n);
  cmsq_vals_.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    const auto TTN = veff_TTN_vals_[i];

    alglib::spline1ddiff(veff_ps_interp_, TTN, s_unused, dps, dps2_unused);
    alglib::spline1ddiff(veff_es_interp_, TTN, s_unused, des, des2_unused);
    cpsq_vals_.push_back(dps / des);

    alglib::spline1ddiff(veff_pb_interp_, TTN, s_unused, dpb, dpb2_unused);
    alglib::spline1ddiff(veff_eb_interp_, TTN, s_unused, deb, deb2_unused);
    cmsq_vals_.push_back(dpb / deb);
  }

  // Fit sound speed splines
  alglib::real_1d_array cpsq_array, cmsq_array;

  cpsq_array.setcontent(n, cpsq_vals_.data());
  // alglib::spline1dfit(veff_TTN_array, cpsq_array, basis_size, smooth_fac, cpsq_fit_, rep);
  alglib::spline1dbuildcubic(veff_TTN_array, cpsq_array, cpsq_fit_);

  cmsq_array.setcontent(n, cmsq_vals_.data());
  // alglib::spline1dfit(veff_TTN_array, cmsq_array, basis_size, smooth_fac, cmsq_fit_, rep);
  alglib::spline1dbuildcubic(veff_TTN_array, cmsq_array, cmsq_fit_);

  wNeN_rat_ = wN_ / eN_;

  // check for normalisation issue in eos (adjust tolerance as needed)
  // impacts prefactor for gw spectrum - only changes max amplitude of spectrum
  if (wNeN_rat_ > 2.0) {
    std::cerr << "Warning: Equation of state normalisation issue. wN/eN=" << wNeN_rat_ << " is abnormally large! "
              << "Using wN/eN = 1 + cpsq approximation instead!\n";
    wNeN_rat_ = 1.0 + alglib::spline1dcalc(cpsq_fit_, 1.0); // munu approx
  }

  // calculate alN_bag and alN_munu
  // const auto theta_s_bag = es_val(1.0) - 3.0 * ps_val(1.0);
  // const auto theta_b_bag = eb_val(1.0) - 3.0 * pb_val(1.0);
  // const auto alN_bag = (theta_s_bag - theta_b_bag) / (3.0 * ws_val(1.0));

  // const auto theta_s = es_val(1.0) - ps_val(1.0) / csq_b(1.0);
  // const auto theta_b = eb_val(1.0) - pb_val(1.0) / csq_b(1.0);
  // const auto alN_munu = (theta_s - theta_b) / (3.0 * ws_val(1.0));
  // std::cout << "alN_bag=" << alN_bag << ", alN_munu=" << alN_munu << "\n";
}

// Public:
#ifdef ENABLE_MATPLOTLIB
void PTParams_Veff::plot_thermo(const std::string& filename) const {
    std::map<std::string, std::string> bp;
    bp["label"] = "broken phase";

    std::map<std::string, std::string> sp;
    sp["label"] = "symmetric phase";


    plt::figure_size(1600, 1600);

    plt::subplot2grid(2, 2, 0, 0);
    plt::plot(veff_TTN_vals_, veff_es_vals_, sp);
    plt::plot(veff_TTN_vals_, veff_eb_vals_, bp);
    plt::xlabel("T/TN");
    plt::ylabel("energy density");
    plt::grid(true);
    plt::legend();

    plt::subplot2grid(2, 2, 0, 1);
    plt::plot(veff_TTN_vals_, veff_ps_vals_, sp);
    plt::plot(veff_TTN_vals_, veff_pb_vals_, bp);
    plt::xlabel("T/TN");
    plt::ylabel("pressure");
    plt::grid(true);
    plt::legend();

    plt::subplot2grid(2, 2, 1, 0);
    plt::plot(veff_TTN_vals_, veff_ws_vals_, sp);
    plt::plot(veff_TTN_vals_, veff_wb_vals_, bp);
    plt::xlabel("T/TN");
    plt::ylabel("enthalpy");
    plt::grid(true);
    plt::legend();

    plt::save(filename);

    return;
}

void PTParams_Veff::plot_thermo2(const std::string& filename) const {
    const auto n = veff_TTN_vals_.size();
    std::vector<double> es_spline_vals(n), eb_spline_vals(n), ps_spline_vals(n), pb_spline_vals(n), ws_spline_vals(n), wb_spline_vals(n),
                        veff_es_vals_norm(n), veff_eb_vals_norm(n), veff_ps_vals_norm(n), veff_pb_vals_norm(n), veff_ws_vals_norm(n), veff_wb_vals_norm(n);

    // normalisation constants
    const auto es_TN = alglib::spline1dcalc(veff_es_interp_, 1.0);
    const auto eb_TN = alglib::spline1dcalc(veff_es_interp_, 1.0);
    const auto ps_TN = alglib::spline1dcalc(veff_es_interp_, 1.0);
    const auto pb_TN = alglib::spline1dcalc(veff_es_interp_, 1.0);

    for (int i = 0; i < n; i++) {
      const auto TTN = veff_TTN_vals_[i];

      es_spline_vals[i] = alglib::spline1dcalc(veff_es_interp_, TTN) / es_TN;
      eb_spline_vals[i] = alglib::spline1dcalc(veff_eb_interp_, TTN) / eb_TN;

      ps_spline_vals[i] = alglib::spline1dcalc(veff_ps_interp_, TTN) / ps_TN;
      pb_spline_vals[i] = alglib::spline1dcalc(veff_pb_interp_, TTN) / pb_TN;

      ws_spline_vals[i] = alglib::spline1dcalc(veff_ws_interp_, TTN) / (es_TN + ps_TN);
      wb_spline_vals[i] = alglib::spline1dcalc(veff_wb_interp_, TTN) / (eb_TN + pb_TN);

      veff_es_vals_norm[i] = veff_es_vals_[i] / es_TN;
      veff_eb_vals_norm[i] = veff_eb_vals_[i] / eb_TN;
      veff_ps_vals_norm[i] = veff_ps_vals_[i] / ps_TN;
      veff_pb_vals_norm[i] = veff_pb_vals_[i] / pb_TN;
      veff_ws_vals_norm[i] = veff_ws_vals_[i] / (es_TN + ps_TN);
      veff_wb_vals_norm[i] = veff_wb_vals_[i] / (eb_TN + pb_TN);
    }  

    std::map<std::string, std::string> data;
    data["label"] = "data";

    std::map<std::string, std::string> spline;
    spline["label"] = "spline";


    plt::figure_size(2400, 1600);

    plt::subplot2grid(3, 2, 0, 0);
    // plt::plot(veff_TTN_vals_, veff_es_vals_, data);
    plt::plot(veff_TTN_vals_, veff_es_vals_norm, data);
    plt::plot(veff_TTN_vals_, es_spline_vals, spline);
    plt::xlabel("T/TN");
    plt::ylabel("es/es(TN)");
    plt::grid(true);
    plt::legend();

    plt::subplot2grid(3, 2, 0, 1);
    // plt::plot(veff_TTN_vals_, veff_eb_vals_, data);
    plt::plot(veff_TTN_vals_, veff_eb_vals_norm, data);
    plt::plot(veff_TTN_vals_, eb_spline_vals, spline);
    plt::xlabel("T/TN");
    plt::ylabel("eb/eb(TN)");
    plt::grid(true);
    plt::legend();

    plt::subplot2grid(3, 2, 1, 0);
    // plt::plot(veff_TTN_vals_, veff_ps_vals_, data);
    plt::plot(veff_TTN_vals_, veff_ps_vals_norm, data);
    plt::plot(veff_TTN_vals_, ps_spline_vals, spline);
    plt::xlabel("T/TN");
    plt::ylabel("ps/ps(TN)");
    plt::grid(true);
    plt::legend();

    plt::subplot2grid(3, 2, 1, 1);
    // plt::plot(veff_TTN_vals_, veff_pb_vals_, data);
    plt::plot(veff_TTN_vals_, veff_pb_vals_norm, data);
    plt::plot(veff_TTN_vals_, pb_spline_vals, spline);
    plt::xlabel("T/TN");
    plt::ylabel("pb/pb(TN)");
    plt::grid(true);
    plt::legend();

    plt::subplot2grid(3, 2, 2, 0);
    // plt::plot(veff_TTN_vals_, veff_ws_vals_, data);
    plt::plot(veff_TTN_vals_, veff_ws_vals_norm, data);
    plt::plot(veff_TTN_vals_, ws_spline_vals, spline);
    plt::xlabel("T/TN");
    plt::ylabel("ws/ws(TN)");
    plt::grid(true);
    plt::legend();

    plt::subplot2grid(3, 2, 2, 1);
    // plt::plot(veff_TTN_vals_, veff_wb_vals_, data);
    plt::plot(veff_TTN_vals_, veff_wb_vals_norm, data);
    plt::plot(veff_TTN_vals_, wb_spline_vals, spline);
    plt::xlabel("T/TN");
    plt::ylabel("wb/wb(TN)");
    plt::grid(true);
    plt::legend();

    plt::save(filename);

    return;
}

void PTParams_Veff::plot_csq(const std::string& filename) const {
  const auto n = veff_TTN_vals_.size();
  std::vector<double> cpsq_spline_vals(n), cmsq_spline_vals(n);
  for (int i = 0; i < n; i++) {
    const auto TTN = veff_TTN_vals_[i];
    cpsq_spline_vals[i] = alglib::spline1dcalc(cpsq_fit_, TTN);
    cmsq_spline_vals[i] = alglib::spline1dcalc(cmsq_fit_, TTN);
  }  

  std::map<std::string, std::string> data;
  data["label"] = "data";

  std::map<std::string, std::string> spline;
  spline["label"] = "spline";

  plt::figure_size(1600, 800);

  plt::subplot2grid(1, 2, 0, 0);
  plt::plot(veff_TTN_vals_, cpsq_vals_, data);
  plt::plot(veff_TTN_vals_, cpsq_spline_vals, spline);
  plt::xlabel("T/TN");
  plt::ylabel("cpsq");
  plt::grid(true);
  plt::legend();

  plt::subplot2grid(1, 2, 0, 1);
  plt::plot(veff_TTN_vals_, cmsq_vals_, data);
  plt::plot(veff_TTN_vals_, cmsq_spline_vals, spline);
  plt::xlabel("T/TN");
  plt::ylabel("cmsq");
  plt::grid(true);
  plt::legend();

  plt::save(filename);
}
#endif

void PTParams_Veff::print() const {
  PTParams::print();
  std::cout << "*************************************************\n";
}


} // namespace PhaseTransition