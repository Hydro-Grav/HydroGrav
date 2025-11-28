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
#include "hydrodynamics.hpp"
#include "phasetransition.hpp"
#include "spectrum.hpp"
#include "profile.hpp"
#include "physics.hpp"
#include "maths_ops.hpp"
#include "constants.hpp"

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
    const auto veff_file = dir + "eos.csv";

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
    params_munu.print();

    // Write fluid profiles to disk
    const Hydrodynamics::FluidProfile profile_bag(params_bag);
    // profile_bag.plot("profile_bag.png");
    // profile_bag.write("prof_bag.csv");
    const auto xi_bag = profile_bag.xi_vals();
    const auto v_bag = profile_bag.v_vals();
    const auto w_bag = profile_bag.w_vals();
    const auto T_bag = profile_bag.T_vals();
    const auto la_bag = profile_bag.la_vals();

    const Hydrodynamics::FluidProfile profile_munu(params_munu);
    // profile_munu.plot("profile_munu.png");
    // profile_munu.write("prof_munu.csv");
    const auto xi_munu = profile_munu.xi_vals();
    const auto v_munu = profile_munu.v_vals();
    const auto w_munu = profile_munu.w_vals();
    const auto T_munu = profile_munu.T_vals();
    const auto la_munu = profile_munu.la_vals();

    const Hydrodynamics::FluidProfile profile_veff(params_veff); // veff
    // profile_veff.plot("profile_veff.png");
    // profile_veff.write("prof_veff.csv");
    const auto xi_veff = profile_veff.xi_vals();
    const auto v_veff = profile_veff.v_vals();
    const auto w_veff = profile_veff.w_vals();
    const auto T_veff = profile_veff.T_vals();
    const auto la_veff = profile_veff.la_vals();

    std::cout << "Bag: cpsq=" << params_bag.cpsq() << ", cmsq=" << params_bag.cmsq() << "\n"
              << "mu nu: cpsq=" << params_munu.cpsq() << ", cmsq=" << params_munu.cmsq() << "\n"
              << "Veff: cpsq=" << params_veff.csq_s(1.0) << ", cmsq=" << params_veff.csq_b(1.0) << "\n";

    // params_veff.plot_thermo();
    // params_veff.plot_csq();

    #ifdef ENABLE_MATPLOTLIB
    const auto filename = bp.dir() + "profile_" + profile_bag.mode_str() + ".png";
    Hydrodynamics::plot_profiles(profile_bag, profile_munu, profile_veff, filename);
    #endif
}

// Gravitational wave power spectrum
void example_GW_Spec(const benchmark_point& bp) {
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

    const PhaseTransition::Universe un(Ts, gs, Hs);

    // Define GW spectrum
    const auto kRs_vals = logspace(1e-3, 1e+3, 100);

    const PhaseTransition::PTParams_Bag params_bag(vw, alN_bag, TN, beta, Rs, nuc_type, un, 1.0 / 3.0, 1.0 / 3.0);
    const auto OmegaGW_bag = Spectrum::GWSpec(kRs_vals, params_bag);
    OmegaGW_bag.write(dir + "GWSpec_bag_" + OmegaGW_bag.profile().mode_str() + ".csv");
    OmegaGW_bag.profile().write(dir + "profile_bag_" + OmegaGW_bag.profile().mode_str() + ".csv");

    const PhaseTransition::PTParams_Bag params_munu(vw, alN_munu, TN, beta, Rs, nuc_type, un, cpsq, cmsq);
    const auto OmegaGW_munu = Spectrum::GWSpec(kRs_vals, params_munu);
    OmegaGW_munu.write(dir + "GWSpec_munu_" + OmegaGW_munu.profile().mode_str() + ".csv");
    OmegaGW_munu.profile().write(dir + "profile_munu_" + OmegaGW_munu.profile().mode_str() + ".csv");

    const PhaseTransition::PTParams_Veff params_veff(vw, alN_munu, TN, beta, Rs, nuc_type, un, veff_file);
    const auto OmegaGW_veff = Spectrum::GWSpec(kRs_vals, params_veff);
    OmegaGW_veff.write(dir + "GWSpec_veff_" + OmegaGW_veff.profile().mode_str() + ".csv");
    OmegaGW_veff.profile().write(dir + "profile_veff_" + OmegaGW_veff.profile().mode_str() + ".csv");

    // #ifdef ENABLE_MATPLOTLIB
    // profile_bag.plot(dir + "profile_bag_" + profile_bag.mode_str() + ".png");
    // profile_munu.plot(dir + "profile_munu_" + profile_munu.mode_str() + ".png");
    // profile_veff.plot(dir + "profile_veff_" + profile_veff.mode_str() + ".png");
    // OmegaGW.plot(dir + ".png");
    // OmegaGW_veff.plot(dir + "_veff.png");
    // #endif

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

    for (int i = 0; i < n; ++i) {
        const auto vw = vw_distr(gen);
        const auto alN = std::pow(10.0, log_alN_distr(gen));
        const auto Rs = PhaseTransition::Rs_approx(vw, beta);

        const PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, nuc_type, PhaseTransition::default_universe());
        // const PhaseTransition::PTParams_Veff params(vw, alN, "thermo.csv");
        
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

    const auto kRs_vals = logspace(1e-3, 1e+3, 100);

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

    for (int i = 0; i < n; ++i) {
        const auto vw = vw_distr(gen);
        const auto alN = std::pow(10.0, log_alN_distr(gen));
        const auto beta = Hs * std::pow(10.0, log_betaH_distr(gen));
        const auto Rs = PhaseTransition::Rs_approx(vw, beta);
        
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
            const auto OmegaGW = Spectrum::GWSpec(kRs_vals, params);
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

void compare_dtau_calc(const benchmark_point& bp) {

    const auto kRs_vals = logspace(1e-3, 1e+3, 100);

    // bag
    auto params_bag_ptr = bp.get_PTParams_Bag();
    const auto OmegaGW_bag = Spectrum::GWSpec(kRs_vals, *params_bag_ptr); // dflt dtau
    OmegaGW_bag.write(bp.dir() + "dtau_comparison/gw_bag_" + OmegaGW_bag.profile().mode_str() + ".csv");

    const auto OmegaGW_bag_dtau = Spectrum::GWSpec(kRs_vals, *params_bag_ptr, true); // calc dtau
    OmegaGW_bag_dtau.write(bp.dir() + "dtau_comparison/gw_bag_dtau_" + OmegaGW_bag_dtau.profile().mode_str() + ".csv");
    
    // munu
    auto params_munu_ptr = bp.get_PTParams_munu();
    const auto OmegaGW_munu = Spectrum::GWSpec(kRs_vals, *params_munu_ptr);
    OmegaGW_munu.write(bp.dir() + "dtau_comparison/gw_munu_" + OmegaGW_munu.profile().mode_str() + ".csv");

    const auto OmegaGW_munu_dtau = Spectrum::GWSpec(kRs_vals, *params_munu_ptr, true);
    OmegaGW_munu_dtau.write(bp.dir() + "dtau_comparison/gw_munu_dtau_" + OmegaGW_munu.profile().mode_str() + ".csv");

    // veff
    auto params_veff_ptr = bp.get_PTParams_Veff();
    const auto OmegaGW_veff = Spectrum::GWSpec(kRs_vals, *params_veff_ptr);
    OmegaGW_veff.write(bp.dir() + "dtau_comparison/gw_veff_" + OmegaGW_veff.profile().mode_str() + ".csv");

    const auto OmegaGW_veff_dtau = Spectrum::GWSpec(kRs_vals, *params_veff_ptr, true);
    OmegaGW_veff_dtau.write(bp.dir() + "dtau_comparison/gw_veff_dtau_" + OmegaGW_veff_dtau.profile().mode_str() + ".csv");

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
        "benchmark_pts/BP0/"
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
        "benchmark_pts/BP0/"
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
        "benchmark_pts/BP0/"
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
        "benchmark_pts/BP1/3deft/"
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
        "benchmark_pts/BP1/3deft/"
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
        "benchmark_pts/BP2/3deft/"
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
        "benchmark_pts/BP2/3deft/"
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
        "benchmark_pts/BP3/"
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
        "benchmark_pts/BP4/"
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
        "benchmark_pts/BP5/"
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
        "benchmark_pts/BP6/3deft/"
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
        "benchmark_pts/BP7/"
    );

    /*
    BP1_hyb: alN_bag=0.00821094, alN_munu=0.00827668
    BP2_hyb: alN_bag=0.00386064, alN_munu=0.00387911
    BP3: alN_bag=0.0130244, alN_munu=0.0131376
    BP4: alN_bag=0.0972393, alN_munu=0.103742
    BP5: alN_bag=0.0283029, alN_munu=0.0289454
    BP6: alN_bag=0.107001, alN_munu=0.10778
    */

    std::vector<benchmark_point> bp_list = {BP0_def, BP0_hyb, BP0_det, BP1_hyb, BP1_det, BP2_hyb, BP2_det, BP3, BP4, BP5, BP6, BP7};

    // example_GW_Spec2(BP5);

    // const int i = 10;
    // for (int i = 0; i < bp_list.size(); i++) {
    //     const auto bp = bp_list[i];
    //     example_GW_Spec(bp);

    //     // compare_dtau_calc(bp);

    //     // Hydrodynamics::plot_profiles(fp_bag, fp_munu, fp_veff, "fp_combined_" + bp.name());
    // }

    // no shock found (xi_sol starts before vw so loop in find_shock_idx always continues)
    const auto vw = 0.620376;
    const auto alN = 0.00412303;

    const auto TN = PhaseTransition::dflt_PTParams::TN;
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto Rs = PhaseTransition::Rs_approx(vw, beta);

    const PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, nuc_type, PhaseTransition::default_universe());
    const Hydrodynamics::FluidProfile profile(params);
    profile.plot("test_prof");

    /************************ CLOCK / PROFILER *************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;
    // ProfilerStop();
    /******************************************************************/
    return 0;
}
