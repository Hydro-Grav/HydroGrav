// PhaseTransition.cpp
#include <string>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <unordered_set>
#include <cassert>

#include "phasetransition.hpp"

/*
TO DO:
- update default vals for beta, wp, wm in PTParams
- check if UF_trans is the same for vp and vm
- only inputs to PTParams should be model, alpha, vw and nucleation type i think
*/

namespace PhaseTransition {

// Universe
Universe::Universe()
    : Universe(dflt_universe::T0, dflt_universe::Ts, dflt_universe::H0, dflt_universe::Hs, dflt_universe::g0, dflt_universe::gs) {}

Universe::Universe(double T0, double Ts, double H0, double Hs, double g0, double gs)
    : T0_(T0),
      Ts_(Ts), 
      H0_(H0), 
      Hs_(Hs), 
      g0_(g0), 
      gs_(gs) {}

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

// PTParams
// make new ctor for reading in Veff since we calculate all the params
// alpha only exists for Bag model, maybe change how it is input/stored in PTParams
PTParams::PTParams()
    : PTParams(dflt_PTParams::vw, dflt_PTParams::alpha, dflt_PTParams::beta, dflt_PTParams::dtau, dflt_PTParams::nuc_type, default_universe()) {}

PTParams::PTParams(double vw, double alpha, double beta, double dtau, const char* nuc_type, const Universe& un)
    : universe_(un),
      vw_(vw),
      alN_(alpha),
      beta_(beta),
      Rs_(),
      tau_s_(),
      tau_fin_(),
      dtau_(dtau),
      cpsq_(),
      cmsq_(),
      nuc_type_()
    {
      // check vw
      if (vw < 0.0 ) {
        std::cout << "Warning: vw < 0. Taking |vw| as input instead.";
        vw_ = std::abs(vw);
      } else if (vw == 0.0 || vw >= 1.0) {
        throw std::invalid_argument("Unphysical wall velocity passed into PTParams. Must have 0 < vw < 1.");
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

      // check valid PT duration
      if (dtau < 0.0) {
        std::cout << "Warning: Invalid PT duration (dtau < 0). Using default value dtau=" << dflt_PTParams::dtau << "\n";
        dtau_ = dflt_PTParams::dtau;
      }
      tau_s_ = 1.0 / universe_.Hs();
      tau_fin_ = tau_s_ + dtau_;

      /***** calc model-dependent quantities *****/
      // if (std::string(model) == "bag") {
        cpsq_ = 1.0 / 3.0;
        cmsq_ = cpsq_;
      // } else {
      //     throw std::invalid_argument("Only Bag model has been implemented so far");
      // }

      Rs_ = std::pow(8 * M_PI, 1. / 3.) * vw_ / beta_;
      /*******************************************/
    }

std::ostream& operator<<(std::ostream& os, const PTParams& params) {
    os << "********** Phase Transition parameters **********\n"
       << std::left
      //  << std::setw(35) << "Equation of state:" << params.model_ << "\n"
       << std::setw(35) << "Nucleation type:" << params.nuc_type_ << "\n"
       << std::setw(35) << "Wall velocity:" << "vw=" << params.vw_ << "\n"
       << std::setw(35) << "PT strength parameter:" << "alN=" << params.alN_ << "\n"
       << std::setw(35) << "Transition rate parameter:" << "beta=" << params.beta_ << "\n"
       << std::setw(35) << "PT duration" << "dtau=" << params.dtau_ << "\n"
       << std::setw(35) << "Speed of sound (broken phase):" << "cpsq=" << params.cpsq_ << "\n"
       << std::setw(35) << "Speed of sound (old phase):" << "cmsq=" << params.cmsq_ << "\n"
       << std::setw(35) << "Mean bubble separation:" << "Rs=" << params.Rs_ << "\n"
       << "*************************************************\n";
       
    return os;
}

void PTParams::print() const {
    std::cout << *this;
}

// Private
bool PTParams::is_valid_model(const char* model, const char* allowed_models[], const int n) {
  for (int i = 0; i < n; i++) {
    if (std::strcmp(model, allowed_models[i]) == 0) {
      return true;
    }
  }
  return false;
}

/********************/

} // namespace PhaseTransition