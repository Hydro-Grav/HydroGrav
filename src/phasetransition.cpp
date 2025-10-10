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
    : Universe(dflt_universe::T0, dflt_universe::Ts, dflt_universe::g0, dflt_universe::gs, dflt_universe::H0) {}

Universe::Universe(double Ts, double gs)
    : Universe(dflt_universe::T0, Ts, dflt_universe::g0, gs, dflt_universe::H0) {}

Universe::Universe(double T0, double Ts, double g0, double gs, double H0)
    : T0_(T0), 
      Ts_(Ts), 
      g0_(g0), 
      gs_(gs), 
      H0_(H0), 
      Hs_(std::pow(4.0 * std::pow(M_PI, 3) * gs * std::pow(Ts, 4) / (45.0 * mP * mP), 0.5)) {}

std::ostream& operator<<(std::ostream& os, const Universe& un) {
    os << "************** Universe parameters **************\n"
       << std::left
       << std::setw(20) << " " << std::setw(13) << "Today" << "Start of PT\n"
       << std::setw(20) << " " << std::setw(13) << "-----" << "-----------\n"
       << std::setw(20) << "Temperature:" << "T0=" << std::setw(10) << un.T0() << "Ts=" << un.Ts() << "\n"
       << std::setw(20) << "Hubble constant:" << "H0=" << std::setw(10) << un.H0() << "Hs=" << un.Hs() << "\n"
       << std::setw(20) << "Number of DoF:" << "g0=" << std::setw(10) << un.g0() << "gs=" << un.gs() << "\n"
       << "*************************************************\n";
       
    return os;
}

// some way to combine this with PTParams print()?
void Universe::print() const {
    std::cout << *this;
}

const Universe& default_universe() {
    static Universe u;;
    return u;
}

/************************************ PTParams ************************************/
PTParams::PTParams(double vw, double alN, double TN, double beta, double dtau, const char* nuc_type, const Universe& un)
    : un_(un),
      vw_(vw),
      alN_(alN),
      TN_(TN),
      wNeN_rat_(std::numeric_limits<double>::quiet_NaN()),
      beta_(beta),
      Rs_(),
      tau_s_(),
      tau_fin_(),
      dtau_(dtau),
      nuc_type_(nuc_type) {

      // check valid vw
      if (vw_ < 0.0 ) {
        std::cout << "Warning: vw < 0. Taking |vw| as input instead.";
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

      // check valid sound wave duration
      if (dtau_ < 0.0) {
        std::cout << "Warning: dtau < 0. Taking |dtau| as input instead.";
        dtau_ = std::abs(dtau);
      }

      // check valid bubble nucleation type
      const char* allowed_nuc[] = {"exp", "sim"};
      const auto m = sizeof(allowed_nuc) / sizeof(allowed_nuc[0]);

      if (is_valid_model(nuc_type, allowed_nuc, m)) {
        nuc_type_ = nuc_type;
      } else {
        std::cout << "Warning: Invalid model '" << nuc_type << "' for bubble nucleation. Using default nucleation type (" << dflt_PTParams::nuc_type << ")\n";
        nuc_type_ = dflt_PTParams::nuc_type;
      }

      Rs_ = std::pow(8 * M_PI, 1. / 3.) * vw_ / beta_;

      // define duration of sound waves
      tau_s_ = 1.0 / un_.Hs();
      tau_fin_ = tau_s_ + dtau_;
      
      // Note on Hs_conformal:
      // Usual definition of conformal Hubble rate is Hs_conformal = as*Hs, but here we include a redshift factor 1/a0 (i.e. so this is Hs_conformal as measured today)
      // tau_s should be defined as 1.0 / Hs_conformal (don't understand why!), but this breaks calculation
      // since tau_s is too large -> Si and Ci integrals are const for such large values so dlt=0 -> OmegaGW=0

      // const auto asa0_rat = std::pow(un_.g0() / un_.gs(), 1./3.) * un_.T0() / un_.Ts(); // a_* / a_0
      // const auto Hs_conformal = un_.Hs() * asa0_rat; // conformal Hubble rate at PT as measured today
      // tau_s_ = 1.0 / Hs_conformal;
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
            << std::setw(35) << "PT duration" << "dtau=" << dtau_ << "\n"
            << std::setw(35) << "Mean bubble separation:" << "Rs=" << Rs_ << "\n";
}

// Private:
bool PTParams::is_valid_model(const char* model, const char* allowed_models[], const int n) const {
  for (int i = 0; i < n; i++) {
    if (std::strcmp(model, allowed_models[i]) == 0) {
      return true;
    }
  }
  return false;
}

/********************************** PTParams_Bag **********************************/
PTParams_Bag::PTParams_Bag(double vw, double alN) // Bag model
    : PTParams_Bag(vw, alN, dflt_PTParams::TN, dflt_PTParams::beta, dflt_PTParams::dtau, dflt_PTParams::nuc_type, default_universe(), dflt_PTParams::cpsq, dflt_PTParams::cmsq) {}

PTParams_Bag::PTParams_Bag(double vw, double alN, double TN, double cpsq, double cmsq) // mu nu model
    : PTParams_Bag(vw, alN, dflt_PTParams::TN, dflt_PTParams::beta, dflt_PTParams::dtau, dflt_PTParams::nuc_type, default_universe(), cpsq, cmsq) {}

PTParams_Bag::PTParams_Bag(double vw, double alN, double TN, double beta, double dtau, const char* nuc_type, const Universe& un, double cpsq, double cmsq)
    : PTParams(vw, alN, TN, beta, dtau, nuc_type, un),
      cpsq_(cpsq),
      cmsq_(cmsq) {

      // check valid speed of sound
      if (!is_valid_csq(cpsq_)) {
        throw std::invalid_argument("Unphysical speed of sound in symmetric phase passed into PTParams_Bag. Must have 0 < cpsq < 1.");
      }
      if (!is_valid_csq(cmsq_)) {
        throw std::invalid_argument("Unphysical speed of sound in broken phase passed into PTParams_Bag. Must have 0 < cmsq < 1.");
      }

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

/********************************** PTParams_Veff *********************************/
PTParams_Veff::PTParams_Veff(double vw, double alN, double TN, const std::string& veff_eos_filename)
    : PTParams_Veff(vw, alN, TN, dflt_PTParams::beta, dflt_PTParams::dtau, dflt_PTParams::nuc_type, default_universe(), veff_eos_filename) {}

// TO DO: remove cpsq and cmsq from inputs when cs(T) implemented
PTParams_Veff::PTParams_Veff(double vw, double alN, double TN, double beta, double dtau, const char* nuc_type, const Universe& un, const std::string& veff_eos_filename)
    : PTParams(vw, alN, TN, beta, dtau, nuc_type, un),
      cpsq_(1.0 / 3.0),
      cmsq_(1.0 / 3.0) {

      if (veff_eos_filename.empty()) {
        std::cout << "Warning: Equation of state from effective potential not found. Using Bag equation of state instead.\n";
        throw std::invalid_argument("Must provide veff eos filename to PTParams_Veff.");
        // call bag instead here? not sure if possible with current structure
      }

      std::cout << "Reading in generic equation of state from file: " << veff_eos_filename << "\n"
                << "Note: File must be formated as T, pb, ps, eb, es (comma separated) with header line\n";
      std::ifstream file(veff_eos_filename);
      if (!file) {
          throw std::runtime_error("Could not open file " + veff_eos_filename);
      }

      std::string line;
      std::getline(file, line); // Skip header

      const auto TN_inv = 1.0 / TN_;

      // more efficient way to read in?
      while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::array<double,5> values;
        std::string token;

        for (auto& val : values) {
            if (!std::getline(ss, token, ',')) {
                throw std::runtime_error("Malformed line in " + veff_eos_filename + ": " + line);
            }
            val = std::stod(token);
        }

        veff_TTN_vals_.push_back(values[0] * TN_inv); // T/TN

        const auto pb_val = values[1];
        const auto ps_val = values[2];
        const auto eb_val = values[3];
        const auto es_val = values[4];

        veff_pb_vals_.push_back(pb_val);
        veff_ps_vals_.push_back(ps_val);
        veff_eb_vals_.push_back(eb_val);
        veff_es_vals_.push_back(es_val);

        veff_ws_vals_.push_back(es_val + ps_val); // w=e+p
        veff_wb_vals_.push_back(eb_val + pb_val);
      }

      // construct interpolating functions p(T/TN), e(T/TN) in s/b phases
      const auto nT = veff_TTN_vals_.size();
      alglib::real_1d_array veff_TTN_array, veff_ps_array, veff_es_array, veff_ws_array, veff_pb_array, veff_eb_array, veff_wb_array;
      veff_TTN_array.setcontent(nT, veff_TTN_vals_.data());

      veff_ps_array.setcontent(nT, veff_ps_vals_.data()); // symmetric phase
      veff_es_array.setcontent(nT, veff_es_vals_.data());
      veff_ws_array.setcontent(nT, veff_ws_vals_.data());
      alglib::spline1dbuildcubic(veff_TTN_array, veff_ps_array, veff_ps_interp_);
      alglib::spline1dbuildcubic(veff_TTN_array, veff_es_array, veff_es_interp_);
      alglib::spline1dbuildcubic(veff_TTN_array, veff_ws_array, veff_ws_interp_);

      veff_pb_array.setcontent(nT, veff_pb_vals_.data()); // broken phase
      veff_eb_array.setcontent(nT, veff_eb_vals_.data());
      veff_wb_array.setcontent(nT, veff_wb_vals_.data());
      alglib::spline1dbuildcubic(veff_TTN_array, veff_pb_array, veff_pb_interp_);
      alglib::spline1dbuildcubic(veff_TTN_array, veff_eb_array, veff_eb_interp_);
      alglib::spline1dbuildcubic(veff_TTN_array, veff_wb_array, veff_wb_interp_);

      // define thermo quantities at nucleation temperature (T/TN = 1)
      pN_ = alglib::spline1dcalc(veff_ps_interp_, 1.0); // pN = p(T/TN=1)
      if (pN_ <= 0.0) {
        throw std::invalid_argument("Unphysical nucleation pressure passed into PTParams. Must have eN > 0.");
      }

      eN_ = alglib::spline1dcalc(veff_es_interp_, 1.0); // eN = e(T/TN=1)
      if (eN_ <= 0.0) {
        throw std::invalid_argument("Unphysical nucleation energy density passed into PTParams. Must have eN > 0.");
      }

      wN_ = alglib::spline1dcalc(veff_ws_interp_, 1.0); // wN = w(T/TN=1)
      if (wN_ <= 0.0) {
        throw std::invalid_argument("Unphysical nucleation enthalpy passed into PTParams. Must have wN > 0.");
      }

      wNeN_rat_ = wN_ / eN_;

      // sound speed in fluid
      double s_unused, ds2_unused;
      double dpsdT, desdT, dpbdT, debdT;

      for (const auto& TTN : veff_TTN_vals_) {
        alglib::spline1ddiff(veff_ps_interp_, TTN, s_unused, dpsdT, ds2_unused);
        alglib::spline1ddiff(veff_es_interp_, TTN, s_unused, desdT, ds2_unused);
        cpsq_vals_.push_back(dpsdT / desdT);

        alglib::spline1ddiff(veff_pb_interp_, TTN, s_unused, dpbdT, ds2_unused);
        alglib::spline1ddiff(veff_eb_interp_, TTN, s_unused, debdT, ds2_unused);
        cmsq_vals_.push_back(dpbdT / debdT);
      }

      alglib::real_1d_array cpsq_array, cmsq_array;

      cpsq_array.setcontent(nT, cpsq_vals_.data());
      alglib::spline1dbuildcubic(veff_TTN_array, cpsq_array, cpsq_interp_);

      cmsq_array.setcontent(nT, cmsq_vals_.data());
      alglib::spline1dbuildcubic(veff_TTN_array, cmsq_array, cmsq_interp_);

      std::cout << "Equation of state read successfully!\n";
    }

// Public:
#ifdef ENABLE_MATPLOTLIB
void PTParams_Veff::plot_thermo(const std::string& filename) const {
    plt::figure_size(2400, 1600);

    plt::subplot2grid(3, 2, 0, 0);
    plt::plot(veff_TTN_vals_, veff_es_vals_);
    plt::xlabel("T/TN");
    plt::ylabel("es(T/TN)");
    plt::grid(true);

    plt::subplot2grid(3, 2, 0, 1);
    plt::plot(veff_TTN_vals_, veff_eb_vals_);
    plt::xlabel("T/TN");
    plt::ylabel("eb(T/TN)");
    plt::grid(true);

    plt::subplot2grid(3, 2, 1, 0);
    plt::plot(veff_TTN_vals_, veff_ps_vals_);
    plt::xlabel("T/TN");
    plt::ylabel("ps(T/TN)");
    plt::grid(true);

    plt::subplot2grid(3, 2, 1, 1);
    plt::plot(veff_TTN_vals_, veff_pb_vals_);
    plt::xlabel("T/TN");
    plt::ylabel("pb(T/TN)");
    plt::grid(true);

    plt::subplot2grid(3, 2, 2, 0);
    plt::plot(veff_TTN_vals_, veff_ws_vals_);
    plt::xlabel("T/TN");
    plt::ylabel("ws(T/TN)");
    plt::grid(true);

    plt::subplot2grid(3, 2, 2, 1);
    plt::plot(veff_TTN_vals_, veff_wb_vals_);
    plt::xlabel("T/TN");
    plt::ylabel("wb(T/TN)");
    plt::grid(true);

    plt::save(filename);

    return;
}

void PTParams_Veff::plot_csq(const std::string& filename) const {
  plt::figure_size(1600, 800);

  plt::subplot2grid(1, 2, 0, 0);
  plt::plot(veff_TTN_vals_, cpsq_vals_);
  plt::xlabel("T/TN");
  plt::ylabel("cpsq(T/TN)");
  plt::grid(true);

  plt::subplot2grid(1, 2, 0, 1);
  plt::plot(veff_TTN_vals_, cmsq_vals_);
  plt::xlabel("T/TN");
  plt::ylabel("cmsq(T/TN)");
  plt::grid(true);

  plt::save(filename);
}
#endif

void PTParams_Veff::print() const {
  PTParams::print();
  std::cout << "*************************************************\n";
}


// double PTParams_Veff::cpsq(double TTN) const {
//   if (TTN < TTN_min() || TTN > TTN_max()) {
//     throw std::out_of_range("Temperature T/TN out of range for Veff eos. Must have " + std::to_string(TTN_min()) + " < T/TN < " + std::to_string(TTN_max()));
//   }
//   return alglib::spline1dcalc(cpsq_interp_, TTN);
// }

// double PTParams_Veff::cmsq(double TTN) const {
//   if (TTN < TTN_min() || TTN > TTN_max()) {
//     throw std::out_of_range("Temperature T/TN out of range for Veff eos. Must have " + std::to_string(TTN_min()) + " < T/TN < " + std::to_string(TTN_max()));
//   }
//   return alglib::spline1dcalc(cmsq_interp_, TTN);
// }




} // namespace PhaseTransition