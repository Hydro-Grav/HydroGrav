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

// Fluid profile
void example_FluidProfile(const std::string& filename) {
    const auto vw = 0.9; // detonation
    // const auto vw = 0.4; // deflagration
    // const auto vw = 0.6; // hybrid
    // const auto vw = 0.0378243;
    // const auto alN = 0.183018;
    
    // Will's benchmark point:
    // std::string veff_file = "flynn_eos.csv";
    // const auto Ts = 46.0096;
    // const auto gs = 106.75;
    // const auto alN = 0.120242;
    // const auto Hs = 3.2193e-15; // GeV
    // const auto betaHs = 588.135; // beta/Hs
    // const auto beta = betaHs * Hs;
    // const auto Rs = std::pow(8 * M_PI, 1. / 3.) * vw / beta;
    // const auto dtau = 10.0 * Rs;
    // const auto TN = Ts;
    // const auto nuc_type = "exp";
    
    // Xiao's benchmark point:
    const auto Ts = 53.370765185008004; // GeV
    const auto gs = 106.75;
    const auto alN = 0.11384915003991744;
    const auto beta = 5.794e+12 * (1.0 / 1.52e+24); // s^-1 * Gev/s^-1 = GeV;
    const auto Rs = std::pow(8 * M_PI, 1. / 3.) * vw / beta;
    const auto dtau = 10.0 * Rs;
    const auto TN = Ts; // GeV
    const auto cpsq = 1.0 / 3.0;
    const auto cmsq = cpsq;
    const auto nuc_type = "exp";
    std::string veff_file = "thermo.csv";

    const PhaseTransition::Universe un(Ts, gs);
    const PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, dtau, nuc_type, un, cpsq, cmsq);
    const PhaseTransition::PTParams_Veff params_veff(vw, alN, TN, beta, dtau, nuc_type, un, veff_file);

    un.print();
    params.print();

    // params_veff.plot_thermo();
    // params_veff.plot_csq();

    const Hydrodynamics::FluidProfile profile_bag(params); // bag model
    // profile_bag.plot("profile_bag.png");
    // profile_bag.write("prof_bag.csv");
    const auto xi_bag = profile_bag.xi_vals();
    const auto v_bag = profile_bag.v_vals();
    const auto w_bag = profile_bag.w_vals();
    const auto T_bag = profile_bag.T_vals();
    const auto la_bag = profile_bag.la_vals();

    const Hydrodynamics::FluidProfile profile_veff(params_veff); // veff
    // profile_veff.plot("profile_veff.png");
    // profile_veff.write("prof_veff.csv");
    const auto xi_veff = profile_veff.xi_vals();
    const auto v_veff = profile_veff.v_vals();
    const auto w_veff = profile_veff.w_vals();
    const auto T_veff = profile_veff.T_vals();
    const auto la_veff = profile_veff.la_vals();

    // plt::figure_size(2400, 800);

    // // v(xi)
    // plt::subplot2grid(2, 2, 0, 0);
    // plt::plot(xi_bag, v_bag, "r-");
    // plt::plot(xi_veff, v_veff, "b-");
    // plt::xlabel("xi");
    // plt::ylabel("v(xi)");
    // plt::xlim(0.0, 1.0);
    // plt::grid(true);

    // // w(xi)
    // plt::subplot2grid(2, 2, 0, 1);
    // plt::plot(xi_bag, w_bag, "r-");
    // plt::plot(xi_veff, w_veff, "b-");
    // plt::xlabel("xi");
    // plt::ylabel("w(xi)");
    // plt::xlim(0.0, 1.0);
    // plt::grid(true);

    // // T(xi)
    // plt::subplot2grid(2, 2, 1, 0);
    // plt::plot(xi_bag, T_bag, "r-");
    // plt::plot(xi_veff, T_veff, "b-");
    // plt::xlabel("xi");
    // plt::ylabel("T(xi)");
    // plt::xlim(0.0, 1.0);
    // plt::grid(true);

    // // la(xi)
    // plt::subplot2grid(2, 2, 1, 1);
    // plt::plot(xi_bag, la_bag, "r-");
    // plt::plot(xi_veff, la_veff, "b-");
    // plt::xlabel("xi");
    // plt::ylabel("la(xi)");
    // plt::xlim(0.0, 1.0);
    // plt::grid(true);

    // plt::suptitle("vw = " + to_string_with_precision(vw) + ", alpha = " + to_string_with_precision(alN));
    // plt::save(filename);

    std::cout << "Fluid profile saved to " << filename << "\n";
}

// // Gravitational wave power spectrum
void example_GW_Spec(const std::string& filename) {
    const auto vw = 0.9; // detonation
    // const auto vw = 0.4; // deflagration
    // const auto vw = 0.6; // hybrid

    
    // Will's benchmark point:
    // const auto Ts = 46.0096;
    // const auto gs = 106.75;

    // const auto alN = 0.120242;
    // const auto Hs = 3.2193e-15; // GeV
    // const auto betaHs = 588.135; // beta/Hs
    // const auto beta = betaHs * Hs;
    // const auto Rs = std::pow(8 * M_PI, 1. / 3.) * vw / beta;
    // const auto dtau = 10.0 * Rs;
    // const auto TN = Ts;
    // const auto nuc_type = "exp";
    
    // Xiao's benchmark point:
    const auto Ts = 53.370765185008004; // GeV
    const auto gs = 106.75;

    const auto alN = 0.11384915003991744;
    const auto beta = 5.794e+12 * (1.0 / 1.52e+24); // s^-1 * Gev/s^-1 = GeV;
    const auto Rs = std::pow(8 * M_PI, 1. / 3.) * vw / beta;
    const auto dtau = 10.0 * Rs;
    const auto TN = Ts; // GeV
    const auto cpsq = 1.0 / 3.0;
    const auto cmsq = cpsq;
    const auto nuc_type = "exp";
    std::string veff_file = "thermo.csv";

    const PhaseTransition::Universe un(Ts, gs);
    // const PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, dtau, nuc_type, un, cpsq, cmsq);
    const PhaseTransition::PTParams_Veff params(vw, alN, TN, beta, dtau, nuc_type, un, veff_file);

    un.print();
    params.print();

    // Define GW spectrum
    // const auto freq_vals = logspace(1e-7, 1.0, 100); // Hz
    // std::vector<double> kRs_vals;
    // for (const auto& f : freq_vals) {
    //     const auto kRs = f * (10.0 / 2.6e-5) * (un.Hs() * params.Rs()) * (100.0 / un.Ts()) * std::pow(100.0 / un.gs(), 1.0 / 6.0);
    //     kRs_vals.push_back(kRs);
    // }
    const auto kRs_vals = logspace(1e-3, 1e+3, 100);
    const auto OmegaGW = Spectrum::GWSpec2(kRs_vals, params);
    // const Hydrodynamics::FluidProfile profile(params);
    // const auto OmegaGW = Spectrum::GWSpec2(kRs_vals, profile);
    
    // Write/plot to disk
    OmegaGW.write(filename + ".csv");
    #ifdef ENABLE_MATPLOTLIB
    OmegaGW.plot(filename + ".png");
    // plt::figure_size(800, 600);
    // plt::loglog(freq_vals, OmegaGW.P(), "k-");
    // plt::suptitle("vw = " + to_string_with_precision(params.vw()) + ", alN = " + to_string_with_precision(params.alN()));
    // plt::xlabel("f [Hz]");
    // plt::ylabel("Omega_GW(f)");
    // plt::xlim(freq_vals.front(), freq_vals.back());
    // plt::grid(true);
    // plt::save(filename + ".png");
    #endif

    return;
}

// // Tests parameter space (vw, alN) for fluid profile calculation
void test_FluidProfile(const std::string& filename = "fluid_profile_test.csv") {
    std::cout << "Running fluid profiles tests for (vw, alN) parameter space...\n";

    const int n = 30000;
    // const int n = 100;

    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    
    std::uniform_real_distribution<double> vw_distr(0.001, 0.999);
    std::uniform_real_distribution<double> log_alN_distr(-6.0, 0.0);

    const std::array<std::string, 3> unphysical_exception = {"alN too small for shock!", "alpha_+ too small for shock", "alpha_+ too large for shock"};
    int pass_count = 0;
    int unphysical_count = 0;

    std::ofstream file(filename);
    file << "vw,alN,mode\n";

    // Suppress console output during testing
    std::streambuf* original_cout_buffer = std::cout.rdbuf();
    std::cout.rdbuf(nullptr);

    for (int i = 0; i < n; ++i) {
        const auto vw = vw_distr(gen);
        const auto alN = std::pow(10.0, log_alN_distr(gen));

        const PhaseTransition::PTParams_Bag params(vw, alN);
        // const PhaseTransition::PTParams_Veff params(vw, alN, "thermo.csv");
        
        file << vw << "," << alN << ",";

        try {
            const Hydrodynamics::FluidProfile profile(params);
            pass_count++;
            file << profile.mode_str() << "\n";
        } catch (const std::exception& e) {
            // flags unphysical parameter choices
            if (e.what() == unphysical_exception[0] || e.what() == unphysical_exception[1] || e.what() == unphysical_exception[2]) {
                file << "unphysical\n";
                unphysical_count++;
            } else {
                file << "fail\n";
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
    const PhaseTransition::Universe un(Ts, gs);
    const auto Hs = un.Hs();

    const auto cpsq = 1.0 / 3.0;
    const auto cmsq = cpsq;
    const auto dtau = PhaseTransition::dflt_PTParams::dtau;
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
        
        const PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, dtau, nuc_type, un, cpsq, cmsq);
        // const PhaseTransition::PTParams_Veff params(vw, alN, TN, beta, dtau, nuc_type, un, "thermo.csv");

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
            const auto OmegaGW = Spectrum::GWSpec2(kRs_vals, params);
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

void test_dSiCi_accuracy() {
    const std::vector<std::pair<double, double>> test_ranges = {
        {0.1, 1.0}, {1.0, 10.0}, {5.0, 15.0}, {0.01, 0.1}, {1e-6, 1e-3}
    };
    const double tolerance = 1e-10;

    for (const auto& [y, x] : test_ranges) {
        // Custom implementation
        std::vector<double> res_custom = dSiCi(x, y, 2000);
        double dSi_custom = res_custom[0];
        double dCi_custom = res_custom[1];

        // ALGLIB implementation
        double Si_x, Ci_x, Si_y, Ci_y;
        alglib::sinecosineintegrals(x, Si_x, Ci_x);
        alglib::sinecosineintegrals(y, Si_y, Ci_y);
        double dSi_lib = Si_x - Si_y;
        double dCi_lib = Ci_x - Ci_y;

        // Check results
        if (std::abs(dSi_custom - dSi_lib) > tolerance ||
            std::abs(dCi_custom - dCi_lib) > tolerance) {
            std::cerr << "Test failed for range (" << y << ", " << x << ")\n";
            std::cerr << std::setprecision(15) << "  dSi: custom = " << dSi_custom << ", lib = " << dSi_lib << "\n";
            std::cerr << std::setprecision(15) << "  dCi: custom = " << dCi_custom << ", lib = " << dCi_lib << "\n";
            assert(false && "Mismatch exceeds tolerance.");
        }
    }

    std::cout << "All dSiCi tests passed successfully.\n";
}

int main() {
    /************************ CLOCK / PROFILER *************************/
    // ProfilerStart("profile.out");
    const auto ti = std::chrono::high_resolution_clock::now();
    /******************************************************************/

    /*** targets for k, p, z, tau: ***/
    // const auto kRs_vals = logspace(1e-3, 1e+3, 100);
    // const auto pRs_vals = logspace(1e-2, 1e+3, 1000);
    // const auto Ttilde_vals = logspace(1e-2, 20, 1000);
    // const auto z_vals = linspace(-1.0, 1.0, 1000);
    // python code takes 2-4mins to run
    /****************************/

    // test_profile_params();
    // example_Kin_Spec("Ekin");
    // example_GW_Spec("GWSpec");
    example_FluidProfile("profile_combined.png");
    // test_FluidProfile("profile_test_bag.csv");
    // test_GWSpec("GWSpec_test_bag.csv");
    // test_rk4_solver();

    /************************ CLOCK / PROFILER *************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;
    // ProfilerStop();
    /******************************************************************/
    return 0;
}
