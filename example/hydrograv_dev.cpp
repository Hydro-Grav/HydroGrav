// main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include <cassert>
#include <chrono>
// #include <gperftools/profiler.h>
#include <omp.h>
#include <streambuf>
#include <random>
#include <memory>

// modify include list when testing of program finished - currently includes everything
#include "phasetransition.hpp"
#include "ssm.hpp"
#include "profile.hpp"
#include "physics.hpp"
#include "maths_ops.hpp"
#include "constants.hpp"
#include "snr.hpp"

#include "ap.h"
#include "interpolation.h"
#include "specialfunctions.h"

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

// Fluid profile
void example_FluidProfile(const benchmark_point& bp) {
    const auto dir = bp.dir();
    const auto veff_file = dir;

    const auto vw = bp.vw();
    const auto Ts = bp.Ts();
    const auto TN = Ts;
    const auto alN_bag = bp.alN_bag();
    const auto alN_munu = bp.alN_munu();
    const auto beta = bp.beta();
    const auto Rs = bp.Rs();
    const auto Hs = bp.Hs();
    const auto cpsq = bp.cpsq();
    const auto cmsq = bp.cmsq();
    const auto gs = bp.gs();
    const auto nuc_type = bp.nuc_type();

    const PhaseTransition::Universe un(Ts, gs, Hs);

    const PhaseTransition::PTParams_Bag params_bag(vw, alN_bag, TN, beta, Rs, nuc_type, un, 1.0 / 3.0, 1.0 / 3.0);
    const PhaseTransition::PTParams_Bag params_munu(vw, alN_munu, TN, beta, Rs, nuc_type, un, cpsq, cmsq);
    const PhaseTransition::PTParams_Veff params_veff(vw, alN_munu, TN, beta, Rs, nuc_type, un, veff_file);

    un.print();
    params_veff.print();

    const Hydrodynamics::FluidProfile profile_bag(params_bag);
    // profile_bag.plot("profile_bag.png");
    // profile_bag.write("prof_bag.csv");

    const Hydrodynamics::FluidProfile profile_munu(params_munu);
    // profile_munu.plot("fp_" + bp.name() + "_new.png");
    // profile_munu.write("prof_munu.csv");

    const Hydrodynamics::FluidProfile profile_veff(params_veff); // veff
    // profile_veff.plot("profile_veff.png");
    // profile_veff.write("prof_veff.csv");

    // std::cout << "Bag: cpsq=" << params_bag.cpsq() << ", cmsq=" << params_bag.cmsq() << "\n"
    //           << "mu nu: cpsq=" << params_munu.cpsq() << ", cmsq=" << params_munu.cmsq() << "\n"
    //           << "Veff: cpsq=" << params_veff.csq_s(1.0) << ", cmsq=" << params_veff.csq_b(1.0) << "\n";

    // params_veff.plot_thermo();
    // params_veff.plot_csq();

    #ifdef ENABLE_MATPLOTLIB
    const auto filename = bp.dir() + "profile_" + profile_bag.mode_str() + ".png";
    // const auto filename = "fp.png";
    Hydrodynamics::plot_profiles(profile_bag, profile_munu, profile_veff, filename);
    #endif
}

// Gravitational wave power spectrum
void example_GW_Spec(const benchmark_point& bp) {
    const auto dir = bp.dir();
    // const auto veff_file = dir + "eos.csv";
    // const std::string dir = "parameter_scan/eos_scan/";
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
    const auto gs = bp.gs();
    const auto nuc_type = bp.nuc_type();

    const auto dtau = 10 * Rs;

    const PhaseTransition::Universe un(Ts, gs, Hs);

    // Define GW spectrum
    const auto kRs_vals = logspace(-3.0, 3.0, 100);

        const PhaseTransition::PTParams_Bag params_bag(vw, alN_bag, TN, beta, Rs, nuc_type, un, 1.0 / 3.0, 1.0 / 3.0);
        // const auto OmegaGW_bag = Spectrum::GWSpec(kRs_vals, params_bag, dtau);
        const auto OmegaGW_bag = Spectrum::GWSpec(kRs_vals, params_bag, dtau);
        OmegaGW_bag.write(dir + "GWSpec_bag_" + OmegaGW_bag.profile().mode_str() + ".csv");
        OmegaGW_bag.profile().write(dir + "profile_bag_" + OmegaGW_bag.profile().mode_str() + ".csv");

        const PhaseTransition::PTParams_Bag params_munu(vw, alN_munu, TN, beta, Rs, nuc_type, un, cpsq, cmsq);
        const auto OmegaGW_munu = Spectrum::GWSpec(kRs_vals, params_munu, dtau);
        OmegaGW_munu.write(dir + "GWSpec_munu_" + OmegaGW_munu.profile().mode_str() + ".csv");
        OmegaGW_munu.profile().write(dir + "profile_munu_" + OmegaGW_munu.profile().mode_str() + ".csv");

        const PhaseTransition::PTParams_Veff params_veff(vw, alN_munu, TN, beta, Rs, nuc_type, un, veff_file);
        const auto OmegaGW_veff = Spectrum::GWSpec(kRs_vals, params_veff, dtau);
        OmegaGW_veff.write(dir + "GWSpec_veff_" + OmegaGW_veff.profile().mode_str() + ".csv");
        OmegaGW_veff.profile().write(dir + "profile_veff_" + OmegaGW_veff.profile().mode_str() + ".csv");

    #ifdef ENABLE_MATPLOTLIB
    // const std::string filename_fp = dir + "fp_" + OmegaGW_veff.profile().mode_str() + ".png";
    // const std::string filename_gw = dir + "gw_" + OmegaGW_veff.profile().mode_str() + ".png";
    const std::string filename_fp = "fp_" + bp.name() + ".png";
    const std::string filename_gw = "gw_" + bp.name() + ".png";
    Hydrodynamics::plot_profiles(OmegaGW_bag.profile(), OmegaGW_munu.profile(), OmegaGW_veff.profile(), filename_fp);
    Spectrum::plot_spectra(OmegaGW_bag, OmegaGW_munu, OmegaGW_veff, filename_gw);
    #endif

    return;
}

// Tests parameter space (vw, alN) for fluid profile calculation
void test_FluidProfile(const size_t n, const std::string& filename = "fluid_profile_test.csv") {
    std::cout << "Running fluid profiles tests for (vw, alN) parameter space...\n";

    // dflt vals
    const auto TN = PhaseTransition::dflt_PTParams::TN;
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;
    // const PhaseTransition::Universe un();

    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    
    std::uniform_real_distribution<double> vw_distr(0.001, 0.999);
    std::uniform_real_distribution<double> log_alN_distr(-6.0, 0.0);

    const std::array<std::string, 3> unphysical_exception = {"alN too small for shock!", "alpha_+ too small for shock", "alpha_+ too large for shock"};
    int pass_count = 0;
    int unphysical_count = 0;

    std::ofstream file(filename);
    file << "vw,alN,mode,note\n";

    // Suppress console output during testing
    std::streambuf* original_cout_buffer = std::cout.rdbuf();
    std::cout.rdbuf(nullptr);

    for (size_t i = 0; i < n; ++i) {
        const auto vw = vw_distr(gen);
        const auto alN = std::pow(10.0, log_alN_distr(gen));
        const auto Rs = PhaseTransition::Rs_approx(vw, beta);

        const PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, nuc_type, PhaseTransition::default_universe());
        // const PhaseTransition::PTParams_Veff params(vw, alN, TN, beta, Rs, nuc_type, PhaseTransition::default_universe(), "benchmark_pts/SS/BP0/eos.csv");
        
        file << vw << "," << alN << ",";

        try {
            const Hydrodynamics::FluidProfile profile(params);
            pass_count++;
            file << profile.mode_str() << "," << "\n";
        } catch (const std::exception& e) {
            // flags unphysical parameter choices
            if (e.what() == unphysical_exception[0] || e.what() == unphysical_exception[1] || e.what() == unphysical_exception[2]) {
                file << "unphysical" << "," << "\n";
                unphysical_count++;
            } else {
                file << "fail" << "," << e.what() << "\n";
            }
        }
    }
    
    // Restore console output
    std::cout.rdbuf(original_cout_buffer);

    file.close();
    std::cout << "Fluid profile test complete: " << pass_count << "/" << n << " cases passed (" << unphysical_count << " unphysical).\n";
    std::cout << "Results saved to '" << filename << "'.\n";

    return;
}

// Tests parameter space (vw, alN) for fluid profile calculation
void test_GWSpec(const std::string& filename = "GWSpec_test.csv") {
    std::cout << "Running GW spectrum tests for (vw, alN) parameter space...\n";

    const int n = 5000;
    // const int n = 1;

    // fixed params
    const auto gs = PhaseTransition::dflt_universe::gs;
    const auto Ts = PhaseTransition::dflt_universe::Ts; // vary this param too?
    const auto Hs = PhaseTransition::dflt_universe::Hs;
    const PhaseTransition::Universe un(Ts, gs, Hs);

    const auto cpsq = 1.0 / 3.0;
    const auto cmsq = cpsq;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;
    const auto TN = Ts;

    const auto kRs_vals = logspace(-3.0, 3.0, 100);

    // parameter range
    std::uniform_real_distribution<double> vw_distr(0.001, 0.999);
    std::uniform_real_distribution<double> log_alN_distr(-6, 0.0);
    std::uniform_real_distribution<double> log_betaH_distr(1.75, 3.75);

    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator

    const std::array<std::string, 3> unphysical_exception = {"alN too small for shock!", "alpha_+ too small for shock", "alpha_+ too large for shock"};
    const std::string profile_failed = "Fluid profile construction failed, aborting GW calculation!";
    int pass_count = 0;
    int unphysical_count = 0;

    std::ofstream file(filename);
    file << "vw,alN,beta,mode\n";

    // Suppress console output during testing
    std::streambuf* original_cout_buffer = std::cout.rdbuf();
    std::cout.rdbuf(nullptr);

    for (size_t i = 0; i < n; ++i) {
        const auto vw = vw_distr(gen);
        const auto alN = std::pow(10.0, log_alN_distr(gen));
        const auto beta = Hs * std::pow(10.0, log_betaH_distr(gen));
        const auto Rs = PhaseTransition::Rs_approx(vw, beta);
        const auto dtau = 10*Rs;
        
        const PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, nuc_type, un, cpsq, cmsq);
        // const PhaseTransition::PTParams_Veff params(vw, alN, TN, beta, Rs, nuc_type, un, "thermo.csv");

        file << vw << "," << alN << "," << beta << ",";

        // attempt to construct fluid profile
        try {
            const Hydrodynamics::FluidProfile profile(params);
        } catch (const std::exception& e) {
            if (e.what() == unphysical_exception[0] || e.what() == unphysical_exception[1] || e.what() == unphysical_exception[2]) {
                file << "unphysical\n";
                unphysical_count++;
            } else {
                file << "profile failed\n";
            }
            continue;
        }

        // attempt to construct GW spectrum
        try {         
            const auto OmegaGW = Spectrum::GWSpec(kRs_vals, params, dtau);
            pass_count++;
            file << OmegaGW.profile().mode_str() << "\n";
        } catch (const std::exception& e) {
            file << "GWSpec failed\n";
        }
    }
    
    // Restore console output
    std::cout.rdbuf(original_cout_buffer);

    file.close();
    std::cout << "Fluid profile test complete: " << pass_count << "/" << n << " cases passed (" << unphysical_count << " unphysical).\n";
    std::cout << "Results saved to '" << filename << "'.\n";

    return;
}

void gw_param_scan(const std::vector<benchmark_point>& bp_list, const std::string& filename="gw_param_scan.csv") {
    // need to get lagrangian parameters for each benchmark point!!
    std::ofstream file(filename);
    // file << "BP,W_bag,W_munu\n";
    file << "la_s,m_s,la_hs,vw,Ts,TN,alN_bag,alN_munu,beta,Hs,Rs,cpsq,cmsq,gs,nuc_type,"
         << "snr_bag,snr_munu,snr_veff,df_peak_bag,df_peak_munu,dOmega_peak_bag,dOmega_peak_munu,"
         << "L1_bag,L1_munu,L2_bag,L2_munu,W1_bag,W1_munu\n";

    const auto kRs_vals = logspace(-3.0, 3.0, 200);

    for (const auto& bp : bp_list) {
        // define PT parameters
        const auto dir = bp.dir();
        const auto veff_file = dir + "eos.csv";

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
        const auto gs = bp.gs();
        const auto nuc_type = bp.nuc_type();

        const auto dtau = 10*Rs;

        const PhaseTransition::Universe un(Ts, gs, Hs);
        const auto Tyears = 4.0; // No. of LISA observation years

        // bag spectrum
        const PhaseTransition::PTParams_Bag params_bag(vw, alN_bag, TN, beta, Rs, nuc_type, un, 1.0 / 3.0, 1.0 / 3.0);
        const auto OmegaGW_bag = Spectrum::GWSpec(kRs_vals, params_bag, dtau);
        OmegaGW_bag.write(dir + "GWSpec_bag_" + OmegaGW_bag.profile().mode_str() + ".csv");
        const auto snr_bag = LISA_snr(OmegaGW_bag.freq(), OmegaGW_bag.P(), Tyears);

        // mu nu spectrum
        const PhaseTransition::PTParams_Bag params_munu(vw, alN_munu, TN, beta, Rs, nuc_type, un, cpsq, cmsq);
        const auto OmegaGW_munu = Spectrum::GWSpec(kRs_vals, params_munu, dtau);
        OmegaGW_munu.write(dir + "GWSpec_munu_" + OmegaGW_munu.profile().mode_str() + ".csv");
        const auto snr_munu = LISA_snr(OmegaGW_munu.freq(), OmegaGW_munu.P(), Tyears);

        // veff spectrum
        const PhaseTransition::PTParams_Veff params_veff(vw, alN_munu, TN, beta, Rs, nuc_type, un, veff_file);
        const auto OmegaGW_veff = Spectrum::GWSpec(kRs_vals, params_veff, dtau);
        OmegaGW_veff.write(dir + "GWSpec_veff_" + OmegaGW_veff.profile().mode_str() + ".csv");
        const auto snr_veff = LISA_snr(OmegaGW_veff.freq(), OmegaGW_veff.P(), Tyears);

        // calculate peak freq/amplitude differences
        const auto [f_peak_bag, Omega_peak_bag] = OmegaGW_bag.peak_vals();
        const auto [f_peak_munu, Omega_peak_munu] = OmegaGW_munu.peak_vals();
        const auto [f_peak_veff, Omega_peak_veff] = OmegaGW_veff.peak_vals();

        const auto df_peak_bag = std::abs(std::log(f_peak_bag / f_peak_veff));
        const auto dOmega_peak_bag = std::abs(std::log(Omega_peak_bag / Omega_peak_veff));
        const auto df_peak_munu = std::abs(std::log(f_peak_munu / f_peak_veff));
        const auto dOmega_peak_munu = std::abs(std::log(Omega_peak_munu / Omega_peak_veff));

        // log frequencies
        std::vector<double> log_freqs(OmegaGW_bag.freq().size());
        for (size_t i = 0; i < OmegaGW_bag.freq().size(); ++i) {
            if (OmegaGW_bag.freq()[i] != OmegaGW_veff.freq()[i] || OmegaGW_munu.freq()[i] != OmegaGW_veff.freq()[i]) {
                throw std::runtime_error("Frequency arrays do not match!");
            }
            log_freqs[i] = std::log(OmegaGW_bag.freq()[i]);
        }

        // Calculate L1 & L2 norms
        const auto L1_bag = L1_norm(log_freqs, OmegaGW_bag.P(), OmegaGW_veff.P());
        const auto L1_munu = L1_norm(log_freqs, OmegaGW_munu.P(), OmegaGW_veff.P());

        const auto L2_bag = L2_norm(log_freqs, OmegaGW_bag.P(), OmegaGW_veff.P());
        const auto L2_munu = L2_norm(log_freqs, OmegaGW_munu.P(), OmegaGW_veff.P());

        // Calculate Wasserstein distance between spectra
        const auto W1_bag = wasserstein_distance_1d(log_freqs, OmegaGW_bag.P(), log_freqs, OmegaGW_veff.P());
        const auto W1_munu = wasserstein_distance_1d(log_freqs, OmegaGW_munu.P(), log_freqs, OmegaGW_veff.P());

        // file << bp.name() << "," << W1_bag_veff << "," << W1_munu_veff << "\n";
        file << "1,1,1" << "," << vw << "," << Ts << "," << TN << "," << alN_bag << "," << alN_munu
             << "," << beta << "," << Hs << "," << Rs << "," << cpsq << "," << cmsq << "," << gs 
             << "," << nuc_type << "," << snr_bag << "," << snr_munu << "," << snr_veff << "," << df_peak_bag << "," << df_peak_munu << "," << dOmega_peak_bag 
             << "," << dOmega_peak_munu << "," << L1_bag << "," << L1_munu << "," << L2_bag << "," << L2_munu << "," << W1_bag << "," << W1_munu << "\n";

    }

    file.close();
    std::cout << "GW parameter scan complete. Results saved to '" << filename << "'.\n";
    return;
}

void calc_prefacs(const benchmark_point& bp) {
    const PhaseTransition::Universe un(bp.Ts(), bp.gs(), bp.Hs());
    const auto kRs_vals = logspace(-3.0, 3.0, 100);

    const PhaseTransition::PTParams_Bag params_bag(bp.vw(), bp.alN_munu(), bp.Ts(), bp.beta(), bp.Rs(), bp.nuc_type(), un, 1./3., 1./3.);
    const Hydrodynamics::FluidProfile profile_bag(params_bag);
    const auto Ek_bag = Spectrum::Ekin(kRs_vals, profile_bag);
    const auto Ek_bag_max = Ek_bag.peak_vals().second;
    const auto prefac_bag = Spectrum::gw_prefac(kRs_vals, profile_bag);
    const auto wNeN_ratio_bag = params_bag.wNeN_rat();

    const PhaseTransition::PTParams_Bag params_munu(bp.vw(), bp.alN_munu(), bp.Ts(), bp.beta(), bp.Rs(), bp.nuc_type(), un, bp.cpsq(), bp.cmsq());
    const Hydrodynamics::FluidProfile profile_munu(params_munu);
    const auto Ek_munu = Spectrum::Ekin(kRs_vals, profile_munu);
    const auto Ek_munu_max = Ek_munu.peak_vals().second;
    const auto prefac_munu = Spectrum::gw_prefac(kRs_vals, profile_munu);
    const auto wNeN_ratio_munu = params_munu.wNeN_rat();
    
    const PhaseTransition::PTParams_Veff params_veff(bp.vw(), bp.alN_munu(), bp.Ts(), bp.beta(), bp.Rs(), bp.nuc_type(), un, bp.dir());
    const Hydrodynamics::FluidProfile profile_veff(params_veff);
    const auto Ek_veff = Spectrum::Ekin(kRs_vals, profile_veff);
    const auto Ek_veff_max = Ek_veff.peak_vals().second;
    const auto prefac_veff = Spectrum::gw_prefac(kRs_vals, profile_veff);
    const auto wNeN_ratio_veff = params_veff.wNeN_rat();

    std::cout << bp.name() << ":\n";

    std::cout << "Bag: wN/eN=" << wNeN_ratio_bag << ", Ek_max=" << Ek_bag_max << ", prefac=" << prefac_bag << "\n"
              << "mu nu: wN/eN=" << wNeN_ratio_munu << ", Ek_max=" << Ek_munu_max << ", prefac=" << prefac_munu << "\n"
              << "Veff: wN/eN=" << wNeN_ratio_veff << ", Ek_max=" << Ek_veff_max << ", prefac=" << prefac_veff << "\n\n";

    return;
}

int main() {
    /************************ CLOCK / PROFILER *************************/
    // ProfilerStart("profile.out");
    const auto ti = std::chrono::high_resolution_clock::now();
    /******************************************************************/

    const auto gs = 106.75;
    const auto nuc_type = "exp";

    benchmark_point BP0_def(
        0.4, // vw (def)
        53.370765185008004,  // Ts
        0.11384915003991744, // alN_bag
        0.11384915003991744, // alN_munu
        953.267, // beta/Hs
        3.99871e-15, // Hs
        1.0 / 3.0, // cpsq
        1.0 / 3.0, // cmsq
        gs,
        nuc_type,
        "BP0_def",
        "benchmark_pts/SS/BP0/eos.csv"
    );

    benchmark_point BP0_hyb(
        0.6, // vw (hyb)
        53.370765185008004,  // Ts
        0.11384915003991744, // alN_bag
        0.11384915003991744, // alN_munu
        953.267, // beta/Hs
        3.99871e-15, // Hs
        1.0 / 3.0, // cpsq
        1.0 / 3.0, // cmsq
        gs,
        nuc_type,
        "BP0_hyb",
        "benchmark_pts/SS/BP0/eos.csv"
    );

    benchmark_point BP0_det(
        0.8, // vw (det)
        53.370765185008004,  // Ts
        0.11384915003991744, // alN_bag
        0.11384915003991744, // alN_munu
        953.267, // beta/Hs
        3.99871e-15, // Hs
        1.0 / 3.0, // cpsq
        1.0 / 3.0, // cmsq
        gs,
        nuc_type,
        "BP0_det",
        "benchmark_pts/SS/BP0/eos.csv"
    );

    benchmark_point BP1_hyb(
        0.588525, // vw
        98.1547,  // Ts
        // 0.00821205, // alN_bag (Will)
        // 0.00821205, // alN_munu (Will)
        0.00821094, // alN_bag
        0.00827668, // alN_munu
        1106.16, // betaHs
        1.3763e-14, // Hs
        0.574624 * 0.574624, // cpsq
        0.569767 * 0.569767, // cmsq
        gs,
        nuc_type,
        "BP1_hyb",
        "benchmark_pts/SS/BP1/3deft/eos.csv"
    );

    benchmark_point BP1_det(
        0.729893, // vw
        98.1547,  // Ts
        // 0.00821205, // alN_bag (Will)
        // 0.00821205, // alN_munu (Will)
        0.00821094, // alN_bag
        0.00827668, // alN_munu
        1106.16, // betaHs
        1.3763e-14, // Hs
        0.574624 * 0.574624, // cpsq
        0.569767 * 0.569767, // cmsq
        gs,
        nuc_type,
        "BP1_det",
        "benchmark_pts/SS/BP1/3deft/eos.csv"
    );

    benchmark_point BP2_hyb(
        0.594842, // vw (hybrid)
        114.579,  // Ts
        // 0.00386107, // alN_bag (Will)
        // 0.00386107, // alN_munu (Will)
        0.00386064, // alN_bag
        0.00387911, // alN_munu
        1445.85, // betaHs
        1.86566e-14, // Hs
        0.575068 * 0.575068, // cpsq
        0.570803 * 0.570803, // cmsq
        gs,
        nuc_type,
        "BP2_hyb",
        "benchmark_pts/SS/BP2/3deft/eos.csv"
    );

    benchmark_point BP2_det(
        0.754187, // vw (detonation)
        114.579,  // Ts
        // 0.00386107, // alN_bag (Will)
        // 0.00386107, // alN_munu (Will)
        0.00386064, // alN_bag
        0.00387911, // alN_munu
        1445.85, // betaHs
        1.86566e-14, // Hs
        0.575068 * 0.575068, // cpsq
        0.570803 * 0.570803, // cmsq
        gs,
        nuc_type,
        "BP2_det",
        "benchmark_pts/SS/BP2/3deft/eos.csv"
    );

    benchmark_point BP2_det_new(
        0.754187, // vw (detonation)
        114.579,  // Ts
        // 0.00386107, // alN_bag (Will)
        // 0.00386107, // alN_munu (Will)
        0.00386064, // alN_bag
        0.00387911, // alN_munu
        1445.85, // betaHs
        1.86566e-14, // Hs
        0.575068 * 0.575068, // cpsq
        0.570803 * 0.570803, // cmsq
        gs,
        nuc_type,
        "BP2_det_new",
        "benchmark_pts/SS/BP2_new/eos.csv"
    );

    benchmark_point BP3(
        0.605811, // vw (hybrid)
        90.7257,  // Ts
        // 0.0130246, // alN_bag (Will)
        // 0.0130246, // alN_munu (Will)
        0.0130244, // alN_bag
        0.0131376, // alN_munu
        1389.28, // betaHs
        1.18307e-14, // Hs
        0.331032, // cpsq
        0.324375, // cmsq
        gs,
        nuc_type,
        "BP3",
        "benchmark_pts/SS/BP3/eos.csv"
    );

    // NOTE: TmTN < 1 for veff
    benchmark_point BP4(
        0.675122, // vw (hybrid)
        52.9772,  // Ts
        // 0.0972391, // alN_bag (Will)
        // 0.0972391, // alN_munu (Will)
        0.0972393, // alN_bag
        0.103742, // alN_munu
        1231.05, // betaHs
        4.34679e-15, // Hs
        0.56705 * 0.56705, // cpsq
        0.539046 * 0.539046, // cmsq
        gs,
        nuc_type,
        "BP4",
        "benchmark_pts/SS/BP4/eos.csv"
    );

    // NOTE: TmTN < 1 for veff
    benchmark_point BP5(
        0.626002, // vw (hybrid)
        76.2128,  // Ts
        // 0.0283037, // alN_bag (Will)
        // 0.0283037, // alN_munu (Will)
        0.0283029, // alN_bag
        0.0289454, // alN_munu
        941.912, // betaHs
        8.49472e-15, // Hs
        0.572982 * 0.572982, // cpsq
        0.55958 * 0.55958, // cmsq
        gs,
        nuc_type,
        "BP5",
        "benchmark_pts/SS/BP5/eos.csv"
    );

    benchmark_point BP6(
        0.765762, // vw (deflagration)
        47.417670, // Ts
        // 0.107000, // alN_bag (Will)
        // 0.107000, // alN_munu (Will)
        0.107001, // alN_bag
        0.10778, // alN_munu
        604.582851, // betaHs
        3.397981e-15, // Hs
        0.570803 * 0.570803, // cpsq
        0.572998 * 0.572998, // cmsq
        gs,
        nuc_type,
        "BP6",
        "benchmark_pts/SS/BP6/3deft/eos.csv"
    );

    // bag gets mode wrong
    benchmark_point BP7(
        0.568989, // vw (deflagration)
        96.660324, // Ts
        0.009315, // alN_bag (Will)
        0.009315, // alN_munu (Will)
        2084.903549, // betaHs
        8.485469e-15, // Hs
        0.574964 * 0.574964, // cpsq
        0.566480 * 0.566480, // cmsq
        gs,
        nuc_type,
        "BP7",
        "benchmark_pts/SS/BP7/eos.csv"
    );

    /*
    BP1_hyb: alN_bag=0.00821094, alN_munu=0.00827668
    BP2_hyb: alN_bag=0.00386064, alN_munu=0.00387911
    BP3: alN_bag=0.0130244, alN_munu=0.0131376
    BP4: alN_bag=0.0972393, alN_munu=0.103742
    BP5: alN_bag=0.0283029, alN_munu=0.0289454
    BP6: alN_bag=0.107001, alN_munu=0.10778
    */

    benchmark_point scan_BP0(
        0.589995, // vw_def
        // 0.763105, // vw_det
        99.5367, // Ts
        0.00783581, // alN_bag
        0.00789714, // alN_munu
        1020.2448610286425, // beta/Hs
        1.41468e-14, // Hs
        0.574665 * 0.574665, // cpsq
        0.569859 * 0.569859, // cmsq
        gs,
        nuc_type,
        "scan_BP0",
        "parameter_scan/eos_scan/eos_127.860000_0.963933.csv"
    );

    benchmark_point scan_BP1(
        0.589344, // vw_def
        // 0.741504, // vw_det
        86.6057, // Ts
        0.0126859, // alN_bag
        0.0128179, // alN_munu
        1241.5227870085725, // beta/Hs
        1.07671e-14, // Hs
        0.574126 * 0.574126, // cpsq
        0.568808 * 0.568808, // cmsq
        gs,
        nuc_type,
        "scan_BP1",
        "parameter_scan/eos_scan/eos_111.992000_0.877620.csv"
    );

    benchmark_point scan_BP2(
        0.573406, // vw_def
        // 0.73521, // vw_det
        118.632, // Ts
        0.00295754, // alN_bag (will)
        0.002968, // alN_munu (will)
        2350.7934124242875, // beta/Hs
        1.9977e-14, // Hs
        0.575253 * 0.575253, // cpsq
        0.571146 * 0.571146, // cmsq
        gs,
        nuc_type,
        "scan_BP2",
        "parameter_scan/eos_scan/eos_128.431000_0.929978.csv"
    );

    benchmark_point scan_BP3(
        0.573272, // vw_def
        // 0.68414, // vw_det
        96.8226, // Ts
        0.00797781, // alN_bag
        0.00803666, // alN_munu
        1558.4086853285744, // beta/Hs
        1.33927e-14, // Hs
        0.57461 * 0.57461, // cpsq
        0.569731 * 0.569731, // cmsq
        gs,
        nuc_type,
        "scan_BP3",
        "parameter_scan/eos_scan/eos_113.929000_0.870927.csv"
    );

    benchmark_point scan_BP4(
        0.734609, // vw
        93.1892, // Ts
        0.00989212, // alN_bag
        0.00998051, // alN_munu
        1189.7451988446653, // beta/Hs
        1.24293e-14, // Hs
        0.574436 * 0.574436, // cpsq
        0.569383 * 0.569383, // cmsq
        gs,
        nuc_type,
        "scan_BP4",
        "parameter_scan/eos_scan/eos_118.419000_0.907345.csv"
    );

    benchmark_point scan_BP5(
        0.823004, // vw
        133.023, // Ts
        0.00110454, // alN_bag
        0.00110601, // alN_munu
        4862.431475263219, // beta/Hs
        2.50457e-14, // Hs
        0.576027 * 0.576027, // cpsq
        0.572424 * 0.572424, // cmsq
        gs,
        nuc_type,
        "scan_BP5",
        "parameter_scan/eos_scan/eos_132.270000_0.929296.csv"
    );

    benchmark_point scan_BP6(
        0.744623, // vw
        136.461, // Ts
        0.000786381, // alN_bag
        0.000787103, // alN_munu
        8424.964601467567, // beta/Hs
        2.63429e-14, // Hs
        0.576126 * 0.576126, // cpsq
        0.572348 * 0.572348, // cmsq
        gs,
        nuc_type,
        "scan_BP6",
        "parameter_scan/eos_scan/eos_119.531000_0.813313.csv"
    );

    benchmark_point scan_BP7(
        0.763095, // vw
        124.337, // Ts
        0.00202748, // alN_bag
        0.00203198, // alN_munu
        4444.665492475749, // beta/Hs
        2.19158e-14, // Hs
        0.57544 * 0.57544, // cpsq
        0.57153 * 0.57153, // cmsq
        gs,
        nuc_type,
        "scan_BP7",
        "parameter_scan/eos_scan/eos_117.046000_0.827426.csv"
    );

    benchmark_point scan_BP8(
        0.739346, // vw
        118.515, // Ts
        0.00299711, // alN_bag
        0.00300799, // alN_munu
        2234.1627087021725, // beta/Hs
        1.99387e-14, // Hs
        0.575238 * 0.575238, // cpsq
        0.571121 * 0.571121, // cmsq
        gs,
        nuc_type,
        "scan_BP8",
        "parameter_scan/eos_scan/eos_130.203000_0.944763.csv"
    );

    benchmark_point scan2_BP0(
        0.496876,
        39.7773,
        0.073497,
        0.0717017,
        4077.4825558235066,
        2.33459e-15,
        0.568415 * 0.568415,
        0.588193 * 0.588193,
        gs,
        nuc_type,
        "scan2_BP0",
        // "parameter_scan/eos_scan2/eos_88.752800_0.826178.csv"
        "parameter_scan/eos_88.752800_0.826178.csv"
    );

    benchmark_point scan2_BP1(
        0.236704,
        92.8928,
        0.00637658,
        0.0063967,
        11474.277831872767,
        1.2324e-14,
        0.574448 * 0.574448,
        0.569513 * 0.569513,
        gs,
        nuc_type,
        "scan2_BP1",
        "parameter_scan/eos_scan2/eos_87.282700_0.738430.csv"
    );

    benchmark_point scan2_BP2(
        0.567028,
        133.416,
        0.00105982,
        0.00106118,
        4696.549287667861,
        2.51919e-14,
        0.576355 * 0.576355,
        0.572526 * 0.572526,
        gs,
        nuc_type,
        "scan2_BP2",
        "parameter_scan/eos_scan2/eos_135.091000_0.953736.csv"
    );

    benchmark_point scan2_BP3(
        0.423937,
        116.589,
        0.00276923,
        0.00277493,
        7350.383499170812,
        1.9296e-14,
        0.575313 * 0.575313,
        0.571248 * 0.571248,
        gs,
        nuc_type,
        "scan2_BP3",
        "parameter_scan/eos_scan2/eos_100.163000_0.738276.csv"
    );

    benchmark_point scan2_BP4(
        0.221472,
        62.3914,
        0.0185038,
        0.0186734,
        9020.927842358946,
        5.6432e-15,
        0.57224 * 0.57224,
        0.567299 * 0.567299,
        gs,
        nuc_type,
        "scan2_BP4",
        "parameter_scan/eos_scan2/eos_82.596000_0.778311.csv"
    );

    benchmark_point scan2_BP5(
        0.30996,
        94.7739,
        0.00649404,
        0.00652125,
        6523.018620754187,
        1.28244e-14,
        0.574541 * 0.574541,
        0.569659 * 0.569659,
        gs,
        nuc_type,
        "scan2_BP5",
        "parameter_scan/eos_scan2/eos_91.845800_0.753723.csv"
    );

    benchmark_point scan2_BP6(
        0.202468,
        124.473,
        0.00170953,
        0.0017106,
        30035.88816373747,
        2.19571e-14,
        0.575514 * 0.575514,
        0.571717 * 0.571717,
        gs,
        nuc_type,
        "scan2_BP6",
        "parameter_scan/eos_scan2/eos_92.130600_0.668430.csv"
    );

    benchmark_point scan2_BP7(
        0.829277,
        38.6951,
        0.218745,
        0.210724,
        458.3930121295817,
        2.33891e-15,
        0.568339 * 0.568339,
        0.592476 * 0.592476,
        gs,
        nuc_type,
        "scan2_BP7",
        "parameter_scan/eos_scan2/eos_106.197000_0.899633.csv"
    );

    benchmark_point scan2_BP8(
        0.521347,
        120.231,
        0.00250459,
        0.00251085,
        4418.484223460745,
        2.05083e-14,
        0.575373 * 0.575373,
        0.571374 * 0.571374,
        gs,
        nuc_type,
        "scan2_BP8",
        "parameter_scan/eos_scan2/eos_112.013000_0.802738.csv"
    );

    benchmark_point scan2_BP9(
        0.313029,
        94.9333,
        0.00639526,
        0.00642104,
        7039.902381396506,
        1.28664e-14,
        0.574548 * 0.574548,
        0.569676 * 0.569676,
        gs,
        nuc_type,
        "scan2_BP9",
        "parameter_scan/eos_scan2/eos_91.389400_0.751282.csv"
    );

    benchmark_point scan2_BP10(
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
        "scan2_BP10",
        // "parameter_scan/eos_scan2/eos_106.490000_0.880577.csv"
        "parameter_scan/eos_106.490000_0.880577.csv"
    );

    benchmark_point scan2_BP11(
        0.640732,
        83.0458,
        0.0172304,
        0.0174563,
        678.6009460848416,
        9.93357e-15,
        0.573887 * 0.573887,
        0.568395 * 0.568395,
        gs,
        nuc_type,
        "scan2_BP11",
        // "parameter_scan/eos_scan2/eos_124.019000_0.961797.csv"
        "parameter_scan/eos_124.019000_0.961797.csv"
    );

    benchmark_point scan2_BP12(
        0.679833,
        60.4088,
        0.04728,
        0.048188,
        664.3364981890268,
        5.36728e-15,
        0.572074 * 0.572074,
        0.567348 * 0.567348,
        gs,
        nuc_type,
        "scan2_BP12",
        "parameter_scan/eos_scan2/eos_107.167000_0.886585.csv"
    );

    benchmark_point scan2_BP13(
        0.579427,
        124.247,
        0.00217401,
        0.00218047,
        2269.230769230769,
        2.18868e-14,
        0.575278 * 0.575278,
        0.57128 * 0.57128,
        gs,
        nuc_type,
        "scan2_BP13",
        "parameter_scan/eos_scan2/eos_140.291000_1.020820.csv"
    );

    benchmark_point scan2_BP14(
        0.771914,
        50.5187,
        0.10034,
        0.10192,
        370.39982534930135,
        3.84768e-15,
        0.570769 * 0.570769,
        0.570349 * 0.570349,
        gs,
        nuc_type,
        "scan2_BP14",
        "parameter_scan/eos_scan2/eos_111.856000_0.920614.csv"
    );

    benchmark_point scan3_BP0(
        0.424034227585217,
        42.4833159363108,
        0.0524864082888102,
        0.0521341022362761,
        3425.71897510552,
        2.65136823153231e-15,
        0.569099271949071 * 0.569099271949071,
        0.580678905077911 * 0.580678905077911,
        gs,
        nuc_type,
        "scan3_BP0",
        "parameter_scan/eos_scan3/bad_eos/eos_86.463300_0.815641.csv"
    );

    benchmark_point scan3_BP1(
        0.591215456228598,
        43.6364350814534,
        0.0697386824723006,
        0.0695952535344496,
        1861.81925679229,
        2.82256805440367e-15,
        0.569405413616819 * 0.569405413616819,
        0.578317311888028 * 0.578317311888028,
        gs,
        nuc_type,
        "scan3_BP1",
        "parameter_scan/eos_scan3/bad_eos/eos_90.300400_0.828345.csv"
    );

    benchmark_point scan3_BP2(
        0.653061873774833,
        45.4737384453901,
        0.0782957447474071,
        0.0786390701387317,
        1109.54482304556,
        3.081674246438e-15,
        0.569819674068506 * 0.569819674068506,
        0.575315998510338 * 0.575315998510338,
        gs,
        nuc_type,
        "scan3_BP2",
        "parameter_scan/eos_scan3/bad_eos/eos_93.951300_0.840716.csv"
    );

    benchmark_point scan3_BP3(
        0.714907232835199,
        43.5017326427126,
        0.104168292454918,
        0.103874305817539,
        762.744208436858,
        2.84678740435448e-15,
        0.569430234580891 * 0.569430234580891,
        0.578590568293944 * 0.578590568293944,
        gs,
        nuc_type,
        "scan3_BP3",
        "parameter_scan/eos_scan3/bad_eos/eos_97.126100_0.855474.csv"
    );

    benchmark_point scan3_BP4(
        0.714180830765479,
        50.5562498790412,
        0.0779077030974728,
        0.079080634417694,
        484.185141122048,
        3.81570945891244e-15,
        0.570759266512843 * 0.570759266512843,
        0.570324267725689 * 0.570324267725689,
        gs,
        nuc_type,
        "scan3_BP4",
        "parameter_scan/eos_scan3/bad_eos/eos_102.185000_0.871710.csv"
    );

    benchmark_point scan3_BP5(
        0.744648046347456,
        47.9752859785928,
        0.0956335088045258,
        0.0966920287228332,
        435.363580582639,
        3.46056323174485e-15,
        0.570333816555551 * 0.570333816555551,
        0.572383002196562 * 0.572383002196562,
        gs,
        nuc_type,
        "scan3_BP5",
        "parameter_scan/eos_scan3/bad_eos/eos_103.074000_0.878137.csv"
    );

    benchmark_point scan3_BP6(
        0.781296132210045,
        47.1428455093062,
        0.114530743197856,
        0.11560465816251,
        311.885297906255,
        3.36745875576192e-15,
        0.570199039380077 * 0.570199039380077,
        0.573257328971889 * 0.573257328971889,
        gs,
        nuc_type,
        "scan3_BP6",
        "parameter_scan/eos_scan3/bad_eos/eos_107.153000_0.898598.csv"
    );

    benchmark_point scan3_BP7(
        0.548765588044471,
        126.788721271538,
        0.00175159918753548,
        0.00175511448748998,
        3090.90060092291,
        2.2778244074639e-14,
        0.575472173939423 * 0.575472173939423,
        0.571614141750257 * 0.571614141750257,
        gs,
        nuc_type,
        "scan3_BP7",
        "parameter_scan/eos_scan3/bad_eos/eos_119.501000_0.839315.csv"
    );

    benchmark_point scan3_BP8(
        0.55313749583814,
        126.088263677203,
        0.00184340073145435,
        0.0018473660292793,
        2753.10679150686,
        2.25304234185504e-14,
        0.57544931477567 * 0.57544931477567,
        0.571566514056531 * 0.571566514056531,
        gs,
        nuc_type,
        "scan3_BP8",
        "parameter_scan/eos_scan3/bad_eos/eos_121.234000_0.854576.csv"
    );

    benchmark_point scan3_BP9(
        0.562182571359414,
        140.27297767704,
        0.000451703770559099,
        0.000451954212443167,
        7368.52225877582,
        2.78178073618762e-14,
        0.576638720544239 * 0.576638720544239,
        0.572988161588945 * 0.572988161588945,
        gs,
        nuc_type,
        "scan3_BP9",
        "parameter_scan/eos_scan3/bad_eos/eos_122.655000_0.829474.csv"
    );

    benchmark_point scan3_BP10(
        0.55819859823755,
        124.763794685288,
        0.00202014236318521,
        0.00202498946643498,
        2474.71702877477,
        2.20655262818198e-14,
        0.575412049657847 * 0.575412049657847,
        0.571487608970243 * 0.571487608970243,
        gs,
        nuc_type,
        "scan3_BP10",
        "parameter_scan/eos_scan3/bad_eos/eos_123.398000_0.874963.csv"
    );

    benchmark_point scan3_BP11(
        0.560331326707943,
        137.464810537926,
        0.000695102425340838,
        0.000695751416132317,
        5347.25238020574,
        2.67271576109664e-14,
        0.575712912419011 * 0.575712912419011,
        0.57200777409572 * 0.57200777409572,
        gs,
        nuc_type,
        "scan3_BP11",
        "parameter_scan/eos_scan3/bad_eos/eos_125.536000_0.860363.csv"
    );

    benchmark_point scan3_BP12(
        0.561820624700799,
        136.831749909775,
        0.000741970429636302,
        0.00074270989804963,
        4780.98298208269,
        2.64837260343305e-14,
        0.575774595343358 * 0.575774595343358,
        0.572115479512002 * 0.572115479512002,
        gs,
        nuc_type,
        "scan3_BP12",
        "parameter_scan/eos_scan3/bad_eos/eos_128.661000_0.888760.csv"
    );

    benchmark_point scan3_BP13(
        0.580906244923419,
        126.980003343966,
        0.00180937492276206,
        0.00181429273706598,
        1519.48449769515,
        2.28472229809105e-14,
        0.575195639726464 * 0.575195639726464,
        0.571228965647879 * 0.571228965647879,
        gs,
        nuc_type,
        "scan3_BP13",
        "parameter_scan/eos_scan3/bad_eos/eos_148.818000_1.097190.csv"
    );

    benchmark_point scan3_BP14(
        0.57051788160334,
        134.667059569115,
        0.000894747706780823,
        0.000896056374203199,
        2542.90395058529,
        2.5658996207417e-14,
        0.575424298550776 * 0.575424298550776,
        0.571522408659435 * 0.571522408659435,
        gs,
        nuc_type,
        "scan3_BP14",
        "parameter_scan/eos_scan3/bad_eos/eos_150.239000_1.098110.csv"
    );

    benchmark_point scan3_BP15(
        0.581028335518798,
        127.375343815094,
        0.00175790967388128,
        0.00176267755098547,
        1516.28234562835,
        2.29878378934788e-14,
        0.575079946151881 * 0.575079946151881,
        0.571138872336147 * 0.571138872336147,
        gs,
        nuc_type,
        "scan3_BP15",
        "parameter_scan/eos_scan3/bad_eos/eos_150.371000_1.112020.csv"
    );

    benchmark_point scan3_BP16(
        0.612834790344668,
        114.101357387217,
        0.00412426377578814,
        0.00414644653653141,
        683.044309458343,
        1.85061458035062e-14,
        0.574927503070381 * 0.574927503070381,
        0.570620944017831 * 0.570620944017831,
        gs,
        nuc_type,
        "scan3_BP16",
        "parameter_scan/eos_scan3/bad_eos/eos_153.187000_1.162030.csv"
    );

    benchmark_point scan3_BP17(
        0.597688035913066,
        122.938148687558,
        0.0024206435070742,
        0.00242940543861108,
        1011.36088079302,
        2.14349027429355e-14,
        0.575041461142732 * 0.575041461142732,
        0.570966485627116 * 0.570966485627116,
        gs,
        nuc_type,
        "scan3_BP17",
        "parameter_scan/eos_scan3/bad_eos/eos_156.730000_1.184500.csv"
    );

    benchmark_point scan3_BP18(
        0.608827583430373,
        118.736820123508,
        0.00317023280309733,
        0.00318444281719968,
        782.197819931002,
        2.00154715709996e-14,
        0.574978517919574 * 0.574978517919574,
        0.570800995194272 * 0.570800995194272,
        gs,
        nuc_type,
        "scan3_BP18",
        "parameter_scan/eos_scan3/bad_eos/eos_157.389000_1.197710.csv"
    );

    benchmark_point scan3_BP19(
        0.60011975917346,
        122.386620581492,
        0.00251460809469135,
        0.00252405564959401,
        947.197574440905,
        2.12457707364485e-14,
        0.575013475943456 * 0.575013475943456,
        0.570928750245745 * 0.570928750245745,
        gs,
        nuc_type,
        "scan3_BP19",
        "parameter_scan/eos_scan3/bad_eos/eos_158.006000_1.198680.csv"
    );

    benchmark_point scan3_BP20(
        0.58372944603808,
        128.730217899208,
        0.00157068882982505,
        0.00157433556736908,
        1355.65779801245,
        2.34723770061459e-14,
        0.575677851051611 * 0.575677851051611,
        0.571704643833846 * 0.571704643833846,
        gs,
        nuc_type,
        "scan3_BP20",
        "parameter_scan/eos_scan3/bad_eos/eos_158.669000_1.196200.csv"
    );

    benchmark_point scan3_BP21(
        0.603831570353577,
        123.809888969341,
        0.00228385253364145,
        0.0022926182995402,
        693.423786205566,
        2.17350728349524e-14,
        0.574736536836366 * 0.574736536836366,
        0.570581563613331 * 0.570581563613331,
        gs,
        nuc_type,
        "scan3_BP21",
        "parameter_scan/eos_scan3/bad_eos/eos_172.086000_1.352690.csv"
    );

    benchmark_point scan3_BP22(
        0.610153896286683,
        123.066624923648,
        0.00240983782442729,
        0.00241926505671368,
        663.093233187663,
        2.14787125307176e-14,
        0.574724624669337 * 0.574724624669337,
        0.570713419564649 * 0.570713419564649,
        gs,
        nuc_type,
        "scan3_BP22",
        "parameter_scan/eos_scan3/bad_eos/eos_172.123000_1.353910.csv"
    );

    benchmark_point scan4_BP0(
        0.136615, // vw
        126.343, // Ts
        0.00149485, // alN_bag
        0.00149523, // alN_munu
        0.00000000136601 / 2.26138e-14, // beta/Hs
        2.26138e-14, // Hs
        0.575555 * 0.575555, // cpsq
        0.57182 * 0.57182, // cmsq
        gs,
        nuc_type,
        "scan4_BP0",
        "parameter_scan/eos_scan4/bps/eos_89.426800_0.647888.csv"
    );

    benchmark_point scan4_BP1(
        0.144901, // vw
        125.473, // Ts
        0.00156999, // alN_bag
        0.00157043, // alN_munu
        0.00000000128504 / 2.23066e-14, // beta/Hs
        2.23066e-14, // Hs
        0.575537 * 0.575537, // cpsq
        0.57178 * 0.57178, // cmsq
        gs,
        nuc_type,
        "scan4_BP1",
        "parameter_scan/eos_scan4/bps/eos_89.515800_0.651567.csv"
    );

    benchmark_point scan4_BP2(
        0.149231, // vw
        141.24, // Ts
        0.000434583, // alN_bag
        0.000434605, // alN_munu
        1.54959e-10 / 2.82025e-14, // beta/Hs
        2.82025e-14, // Hs
        0.576053 * 0.576053, // cpsq
        0.572553 * 0.572553, // cmsq
        gs,
        nuc_type,
        "scan4_BP2",
        "parameter_scan/eos_scan4/bps/eos_79.952400_0.537056.csv"
    );

    benchmark_point scan4_BP3(
        0.255069, // vw
        139.306, // Ts
        0.000560134, // alN_bag
        0.000560242, // alN_munu
        0.00000000194436 / 2.7442e-14, // beta/Hs
        2.7442e-14, // Hs
        0.576324 * 0.576324, // cpsq
        0.572916 * 0.572916, // cmsq
        gs,
        nuc_type,
        "scan4_BP3",
        "parameter_scan/eos_scan4/bps/eos_86.408700_0.579259.csv"
    );

    benchmark_point scan4_BP4(
        0.796963, // vw
        136.243, // Ts
        0.000786375, // alN_bag
        0.000787255, // alN_munu
        1.52155e-10 / 2.62582e-14, // beta/Hs
        2.62582e-14, // Hs
        0.575658 * 0.575658, // cpsq
        0.57193 * 0.57193, // cmsq
        gs,
        nuc_type,
        "scan4_BP4",
        "parameter_scan/eos_scan4/bps/eos_134.707000_0.944110.csv"
    );

    benchmark_point scan4_BP5(
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
        "scan4_BP5",
        "parameter_scan/eos_scan4/bps/eos_118.884000_0.949788.csv"
    );

    benchmark_point scan4_BP6(
        0.814479, // vw
        134.641, // Ts
        0.000934945, // alN_bag
        0.000936077, // alN_munu
        1.36772e-10 / 2.56511e-14, // beta/Hs
        2.56511e-14, // Hs
        0.576167 * 0.576167, // cpsq
        0.572299 * 0.572299, // cmsq
        gs,
        nuc_type,
        "scan4_BP6",
        "parameter_scan/eos_scan4/bps/eos_133.063000_0.932698.csv"
    );

    benchmark_point scan4_BP7(
        0.523353, // vw
        117.931, // Ts
        0.00283357, // alN_bag
        0.00284157, // alN_munu
        7.85268e-11 / 1.97411e-14, // beta/Hs
        1.97411e-14, // Hs
        0.575317 * 0.575317, // cpsq
        0.571249 * 0.571249, // cmsq
        gs,
        nuc_type,
        "scan4_BP7",
        "parameter_scan/eos_scan4/bps/eos_112.085000_0.809488.csv"
    );

    benchmark_point scan4_BP8(
        0.576963, // vw
        121.954, // Ts
        0.00248277, // alN_bag
        0.00249068, // alN_munu
        4.93057e-11 / 2.10964e-14, // beta/Hs
        2.10964e-14, // Hs
        0.575282 * 0.575282, // cpsq
        0.571248 * 0.571248, // cmsq
        gs,
        nuc_type,
        "scan4_BP8",
        "parameter_scan/eos_scan4/bps/eos_134.370000_0.972618.csv"
    );

    benchmark_point scan4_BP9(
        0.59008, // vw
        87.3169, // Ts
        0.0123841, // alN_bag
        0.0125113, // alN_munu
        1.30721e-11 / 1.10941e-14, // beta/Hs
        1.10941e-14, // Hs
        0.574161 * 0.574161, // cpsq
        0.568868 * 0.568868, // cmsq
        gs,
        nuc_type,
        "scan4_BP9",
        "parameter_scan/eos_scan4/bps/eos_112.856000_0.881694.csv"
    );

    // NOTE: Old BPs will be broken since order of 's' and 'b' in veff eos switched when reading it in -> update veff eos files
    std::vector<benchmark_point> bp_list = {BP0_def, BP0_hyb, BP0_det, BP1_hyb, BP1_det, BP2_hyb, BP2_det, BP2_det_new, BP3, BP4, BP5, BP6, BP7};
    std::vector<benchmark_point> scan_bp_list = {scan_BP0, scan_BP1, scan_BP2, scan_BP3, scan_BP4, scan_BP5, scan_BP6, scan_BP7, scan_BP8};
    std::vector<benchmark_point> scan2_bp_list = {scan2_BP0, scan2_BP1, scan2_BP2, scan2_BP3, scan2_BP4, scan2_BP5, scan2_BP6, scan2_BP7, 
                                                  scan2_BP8, scan2_BP9, scan2_BP10, scan2_BP11, scan2_BP12, scan2_BP13, scan2_BP14};
    std::vector<benchmark_point> scan3_bp_list = {scan3_BP0, scan3_BP1, scan3_BP2, scan3_BP3, scan3_BP4, scan3_BP5, scan3_BP6, scan3_BP7, 
                                                  scan3_BP8, scan3_BP9, scan3_BP10, scan3_BP11, scan3_BP12, scan3_BP13, scan3_BP14, 
                                                  scan3_BP15, scan3_BP16, scan3_BP17, scan3_BP18, scan3_BP19, scan3_BP20, 
                                                  scan3_BP21, scan3_BP22};
    std::vector<benchmark_point> scan4_bp_list = {scan4_BP0, scan4_BP1, scan4_BP2, scan4_BP3, scan4_BP4, scan4_BP5, scan4_BP6, scan4_BP7, scan4_BP8, scan4_BP9};

    // std::vector<std::pair<std::string, std::string>> fail_cases;
    // for (size_t i = 0; i < bp_list.size(); i++) {
    //     const auto bp = bp_list[i];
    //     try {
    //         example_GW_Spec(bp);
    //         // example_FluidProfile(bp);
    //     } catch (std::exception& e) {
    //         std::cout << "Failed on " << bp.name() << "!\n";
    //         fail_cases.push_back({bp.name(), e.what()});
    //     }
    // }

    // std::cout << "Fail Cases:\n";
    // for (size_t i = 0; i < fail_cases.size(); i++) {
    //     std::cout << fail_cases[i].first << ": " << fail_cases[i].second << "\n";
    // }

    // for (size_t i = 0; i < scan2_bp_list.size(); i++) {
        const int i = 9;
        const auto bp = scan4_bp_list[i];
        const PhaseTransition::Universe un(bp.Ts(), bp.gs(), bp.Hs());
        const auto kRs_vals = logspace(-3.0, 3.0, 100);
        const auto dtau = 10*bp.Rs();

        // const PhaseTransition::PTParams_Veff params_veff(bp.vw(), bp.alN_munu(), bp.Ts(), bp.beta(), bp.Rs(), bp.nuc_type(), un, bp.dir());
        // const auto OmegaGW_veff = Spectrum::GWSpec(kRs_vals, params_veff, dtau);
        
        // std::cout << "wNeN_rat=" << params_veff.wNeN_rat() << "\n";
        // OmegaGW_veff.profile().write("fp_" + bp.name() + "_veff.csv");
        // OmegaGW_veff.write("gw_" + bp.name() + "_veff.csv");

        // std::cout << "shock_flag=" << OmegaGW_veff.profile().shock_flag() << "\n";

        // const PhaseTransition::PTParams_Bag params_bag(bp.vw(), bp.alN_bag(), bp.Ts(), bp.beta(), bp.Rs(), bp.nuc_type(), un, 1./3., 1./3.);
        // const auto OmegaGW_bag = Spectrum::GWSpec(kRs_vals, params_bag, dtau);
        
        // OmegaGW_bag.profile().write("fp_" + bp.name() + "_bag.csv");
        // OmegaGW_bag.write("gw_" + bp.name() + "_bag.csv");

        const PhaseTransition::PTParams_Bag params_munu(bp.vw(), bp.alN_munu(), bp.Ts(), bp.beta(), bp.Rs(), bp.nuc_type(), un, bp.cpsq(), bp.cmsq());
        const auto OmegaGW_munu = Spectrum::GWSpec(kRs_vals, params_munu, dtau);
        
        // OmegaGW_munu.profile().write("fp_" + bp.name() + "_munu.csv");
        // OmegaGW_munu.write("gw_" + bp.name() + "_munu.csv");

        // std::vector<double> veff_amp, munu_amp, bag_amp;
        // for (size_t i = 0; i < kRs_vals.size(); i++) {
        //     veff_amp.push_back(h*h*OmegaGW_veff.P()[i]);
        //     munu_amp.push_back(h*h*OmegaGW_munu.P()[i]);
        //     bag_amp.push_back(h*h*OmegaGW_bag.P()[i]);
        // }

        // const auto snr_veff = LISA_snr(OmegaGW_veff.freq(), veff_amp);
        // const auto snr_bag = LISA_snr(OmegaGW_bag.freq(), bag_amp);
        // const auto snr_munu = LISA_snr(OmegaGW_munu.freq(), bag_amp);
        // std::cout << "snr_bag=" << snr_bag << ", snr_munu=" << snr_munu << ", snr_veff=" << snr_veff << "\n";
        

        // #ifdef ENABLE_MATPLOTLIB
        // const std::string filename_fp = "fp_" + bp.name() + ".png";
        // Hydrodynamics::plot_profiles(OmegaGW_bag.profile(), OmegaGW_munu.profile(), OmegaGW_veff.profile(), filename_fp);
        // Spectrum::plot_spectra(OmegaGW_bag, OmegaGW_munu, OmegaGW_veff, "gw_" + bp.name() + ".png");
        // Spectrum::plot_spectra(Ekin_bag, Ekin_munu, Ekin_veff, "ekin_" + bp.name() + ".png");
        // params_veff.plot_thermo("thermo.png");
        // params_veff.plot_csq("csq_" + bp.name() + ".png");
        // #endif
    // }

    // veff solution space stuff
    // const int i = 8;
    // const auto bp = bp_list[i];
    // const PhaseTransition::Universe un(bp.Ts(), bp.gs(), bp.Hs());
    // const PhaseTransition::PTParams_Veff params_veff(bp.vw(), bp.alN_munu(), bp.Ts(), bp.beta(), bp.Rs(), bp.nuc_type(), un, bp.dir() + "eos.csv");
    // Hydrodynamics::veff_solution_space(params_veff);

    /************************ CLOCK / PROFILER *************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;
    // ProfilerStop();
    /******************************************************************/
    return 0;
}
