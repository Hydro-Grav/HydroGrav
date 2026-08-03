#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include <chrono>
#include <memory>

#include "phasetransition.hpp"
#include "ssm.hpp"
#include "profile.hpp"

#ifdef ENABLE_MATPLOTLIB
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

class benchmark_point {
    public:
        benchmark_point(double vw, double Ts, double alN_bag, double alN_munu, double betaHs, double Hs, double cpsq, double cmsq, double gs, const char* nuc_type, const std::string& id, const std::string& dir)
        : benchmark_point(vw, Ts, alN_bag, alN_munu, betaHs, Hs, PhaseTransition::Rs_approx(vw, betaHs * Hs), cpsq, cmsq, gs,nuc_type, id, dir) {}
        benchmark_point(double vw, double Ts, double alN_bag, double alN_munu, double betaHs, double Hs, double Rs, double cpsq, double cmsq, double gs, const char* nuc_type, const std::string& id, const std::string& dir)
        : vw_(vw), Ts_(Ts), alN_bag_(alN_bag), alN_munu_(alN_munu), betaHs_(betaHs), Hs_(Hs), Rs_(Rs), cpsq_(cpsq), cmsq_(cmsq), gs_(gs), nuc_type_(nuc_type), id_(id), dir_(dir) {
            beta_ = betaHs_ * Hs_;
        }
        
        double vw() const { return vw_; }
        double Ts() const { return Ts_; }
        double alN_bag() const { return alN_bag_; }
        double alN_munu() const { return alN_munu_; }
        double betaHs() const { return betaHs_; }
        double Hs() const { return Hs_; }
        double beta() const { return beta_; }
        double cpsq() const { return cpsq_; }
        double cmsq() const { return cmsq_; }
        double gs() const { return gs_; }
        double Rs() const { return Rs_;}
        const char* nuc_type() const { return nuc_type_; }
        std::string name() const { return id_; }
        std::string dir() const { return dir_; }

        void print() const {
            std::cout << "********** Phase Transition parameters **********\n"
                      << std::left
                      << std::setw(35) << "vw ="        << vw_        << "\n"
                      << std::setw(35) << "Ts ="        << Ts_        << "\n"
                      << std::setw(35) << "alN_bag ="   << alN_bag_   << "\n"
                      << std::setw(35) << "alN_munu ="  << alN_munu_  << "\n"
                      << std::setw(35) << "betaHs ="    << betaHs_    << "\n"
                      << std::setw(35) << "Hs ="        << Hs_        << "\n"
                      << std::setw(35) << "beta ="      << beta_      << "\n"
                      << std::setw(35) << "cpsq ="      << cpsq_      << "\n"
                      << std::setw(35) << "cmsq ="      << cmsq_      << "\n"
                      << std::setw(35) << "gs ="        << gs_        << "\n"
                      << std::setw(35) << "Rs ="        << Rs_        << "\n"
                      << std::setw(35) << "nuc_type ="  << nuc_type_  << "\n"
                      << std::setw(35) << "dir ="       << dir_       << "\n"
                      << "*************************************************\n";
        }

        std::unique_ptr<PhaseTransition::PTParams> get_PTParams_Bag() const {
            const PhaseTransition::Universe un(Ts_, gs_, Hs_);
            return std::make_unique<PhaseTransition::PTParams_Bag>(vw_, alN_bag_, Ts_, beta_, Rs_, nuc_type_, un);
        }

        std::unique_ptr<PhaseTransition::PTParams> get_PTParams_munu() const {
            const PhaseTransition::Universe un(Ts_, gs_, Hs_);
            return std::make_unique<PhaseTransition::PTParams_Bag>(vw_, alN_munu_, Ts_, beta_, Rs_, nuc_type_, un, cpsq_, cmsq_);
        }

        std::unique_ptr<PhaseTransition::PTParams> get_PTParams_Veff() const {
            const PhaseTransition::Universe un(Ts_, gs_, Hs_);
            return std::make_unique<PhaseTransition::PTParams_Veff>(vw_, alN_munu_, Ts_, beta_, Rs_, nuc_type_, un, dir_ + "eos.csv");
        }

    private:
        const double vw_, Ts_, alN_bag_, alN_munu_, betaHs_, Hs_, Rs_, cpsq_, cmsq_, gs_;
        const char* nuc_type_;
        const std::string id_, dir_;
        double beta_;
};

int main() {
    /*
        This example script can be used to generate the fluid profiles 
        and gravitaitonal wave spectra that appear in Figures 4 & 5 of 
        arXiv:2606.27775
    */

    const auto ti = std::chrono::high_resolution_clock::now(); // clock

    const auto gs = 106.75;
    const auto nuc_type = "exp";

    // deflagration bp
    benchmark_point bp_def(
        0.496876, // vw
        39.7773, // Ts
        0.073497, // alN_bag
        0.0717017, // alN_munu
        4077.4825558235066, // beta/Hs
        2.33459e-15, // Hs
        0.568415 * 0.568415, // cpsq
        0.588193 * 0.588193, // cmsq
        gs,
        nuc_type,
        "bp_def",
        "bps/eos_88.752800_0.826178.csv"
    );

    // detonation bp
    benchmark_point bp_det(
        0.718823, // vw
        64.3442, // Ts
        0.0456357, // alN_bag
        0.0465414, // alN_munu
        2.39209e-12 / 6.07892e-15, // beta/Hs
        6.07892e-15, // Hs
        0.572472 * 0.572472, // cpsq
        0.567136 * 0.567136, // cmsq
        gs,
        nuc_type,
        "bp_det",
        "bps/eos_118.884000_0.949788.csv"
    );

    // hybrid bp
    benchmark_point bp_hyb(
        0.659366,
        62.6341,
        0.0408514,
        0.0416136,
        761.0359024488668,
        5.74919e-15,
        0.572315 * 0.572315,
        0.56722 * 0.56722,
        gs,
        nuc_type,
        "bp_hyb",
        "bps/eos_106.490000_0.880577.csv"
    );

    const auto bp = bp_def;

    const auto veff_file = bp.dir();
    const auto vw = bp.vw();
    const auto Ts = bp.Ts();
    const auto TN = Ts;
    const auto alN_bag = bp.alN_bag();
    const auto alN_munu = bp.alN_munu();
    const auto beta = bp.beta();
    const auto Hs = bp.Hs();
    const auto Rs = bp.Rs();
    const auto cpsq = bp.cpsq();
    const auto cmsq = bp.cmsq();

    const PhaseTransition::Universe un(Ts, gs, Hs);
    const auto kRs_vals = logspace(-3.0, 3.0, 100);

    const PhaseTransition::PTParams_Bag params_bag(vw, alN_bag, TN, beta, Rs, nuc_type, un, 1.0 / 3.0, 1.0 / 3.0);
    const auto OmegaGW_bag = Spectrum::GWSpec(kRs_vals, params_bag);
    OmegaGW_bag.write("gw_bag_" + OmegaGW_bag.profile().mode_str() + ".csv");
    OmegaGW_bag.profile().write("profile_bag_" + OmegaGW_bag.profile().mode_str() + ".csv");


    const PhaseTransition::PTParams_Bag params_munu(vw, alN_munu, TN, beta, Rs, nuc_type, un, cpsq, cmsq);
    const auto OmegaGW_munu = Spectrum::GWSpec(kRs_vals, params_munu);
    OmegaGW_munu.write("gw_munu_" + OmegaGW_munu.profile().mode_str() + ".csv");
    OmegaGW_munu.profile().write("profile_munu_" + OmegaGW_munu.profile().mode_str() + ".csv");

    const PhaseTransition::PTParams_Veff params_veff(vw, alN_munu, TN, beta, Rs, nuc_type, un, veff_file);
    const auto OmegaGW_veff = Spectrum::GWSpec(kRs_vals, params_veff);
    OmegaGW_veff.write("gw_veff_" + OmegaGW_veff.profile().mode_str() + ".csv");
    OmegaGW_veff.profile().write("profile_veff_" + OmegaGW_veff.profile().mode_str() + ".csv");

    #ifdef ENABLE_MATPLOTLIB
    const std::string filename_fp = "fp_" + bp.name() + ".png";
    const std::string filename_gw = "gw_" + bp.name() + ".png";
    Hydrodynamics::plot_profiles(OmegaGW_bag.profile(), OmegaGW_munu.profile(), OmegaGW_veff.profile(), filename_fp);
    Spectrum::plot_spectra(OmegaGW_bag, OmegaGW_munu, OmegaGW_veff, filename_gw);
    #endif

    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;

    return 0;
}
