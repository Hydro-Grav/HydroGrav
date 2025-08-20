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

// tests program across a large parameter space
void test_FluidProfile_params() {
    // Fluid profile
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto dtau = PhaseTransition::dflt_PTParams::dtau;
    const auto TN = PhaseTransition::dflt_PTParams::TN;
    const auto wNeN_rat = PhaseTransition::dflt_PTParams::wNeN_rat;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    PhaseTransition::Universe un;
    
    const auto vw_vals = linspace(1e-5, 1.0 - 1e-5, 10);
    const auto alN_vals = linspace(1e-5, 1.0, 10);

    for (const auto vw : vw_vals) {
        for (const auto alN : alN_vals) {
            PhaseTransition::PTParams params(vw, alN, beta, dtau, TN, wNeN_rat, nuc_type, un);
            try {
                Hydrodynamics::FluidProfile prof(params);
            } catch (std::runtime_error& e) {
                std::cout << "Failed for vw=" << vw << ", alN=" << alN << ":\n";
                std::cout << e.what() << "\n";
            } catch (std::invalid_argument& e) {
                std::cout << e.what() << "\n";
            }
        }
    }

    std::cout << "Parameter test for FluidProfile passed!\n";    
}

// Fluid profile
void example_FluidProfie(const std::string& filename) {
    // const auto vw = 0.9; // detonation
    const auto vw = 0.4; // deflagration
    // const auto vw = 0.6; // hybrid

    const auto alN = 0.11384915003991744;
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto dtau = PhaseTransition::dflt_PTParams::dtau;
    const auto TN = 53.370765185008004;
    const auto wNeN_rat = PhaseTransition::dflt_PTParams::wNeN_rat;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    const PhaseTransition::Universe un;
    const PhaseTransition::PTParams params(vw, alN, beta, dtau, TN, wNeN_rat, nuc_type, un);
    params.print();

    // read veff data
    std::vector<double> veff_T_vals, veff_ps_vals, veff_pb_vals, veff_es_vals, veff_eb_vals;
    std::string veff_file = "thermo.csv";
    std::ifstream file(veff_file);
    if (!file) {
        throw std::runtime_error("Could not open file " + veff_file);
    }

    std::string line;
    std::getline(file, line); // Skip header

    state_type xi_vals, v_vals, w_vals, la_vals;

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::array<double,5> values;
        std::string token;

        for (auto& val : values) {
            if (!std::getline(ss, token, ',')) {
                throw std::runtime_error("Malformed line in " + veff_file + ": " + line);
            }
            val = std::stod(token);
        }

        veff_T_vals.push_back(values[0]);
        veff_pb_vals.push_back(values[1]);
        veff_ps_vals.push_back(values[2]);
        veff_eb_vals.push_back(values[3]);
        veff_es_vals.push_back(values[4]);
    }

    const Hydrodynamics::FluidProfile profile_bag(params); // bag model
    // profile_bag.plot("profile_bag.png");
    const auto xi_bag = profile_bag.xi_vals();
    const auto v_bag = profile_bag.v_vals();
    const auto w_bag = profile_bag.w_vals();
    const auto T_bag = profile_bag.T_vals();
    const auto la_bag = profile_bag.la_vals();

    const Hydrodynamics::FluidProfile profile_veff(params, veff_T_vals, veff_ps_vals, veff_pb_vals, veff_es_vals, veff_eb_vals); // veff
    // profile_veff.plot("profile_veff.png");
    const auto xi_veff = profile_veff.xi_vals();
    const auto v_veff = profile_veff.v_vals();
    const auto w_veff = profile_veff.w_vals();
    const auto T_veff = profile_veff.T_vals();
    const auto la_veff = profile_veff.la_vals();

    plt::figure_size(2400, 800);

    // v(xi)
    plt::subplot2grid(2, 2, 0, 0);
    plt::plot(xi_bag, v_bag, "r-");
    plt::plot(xi_veff, v_veff, "b-");
    plt::xlabel("xi");
    plt::ylabel("v(xi)");
    plt::xlim(0.0, 1.0);
    plt::grid(true);

    // w(xi)
    plt::subplot2grid(2, 2, 0, 1);
    plt::plot(xi_bag, w_bag, "r-");
    plt::plot(xi_veff, w_veff, "b-");
    plt::xlabel("xi");
    plt::ylabel("w(xi)");
    plt::xlim(0.0, 1.0);
    plt::grid(true);

    // T(xi)
    plt::subplot2grid(2, 2, 1, 0);
    plt::plot(xi_bag, T_bag, "r-");
    plt::plot(xi_veff, T_veff, "b-");
    plt::xlabel("xi");
    plt::ylabel("T(xi)");
    plt::xlim(0.0, 1.0);
    plt::grid(true);

    // la(xi)
    plt::subplot2grid(2, 2, 1, 1);
    plt::plot(xi_bag, la_bag, "r-");
    plt::plot(xi_veff, la_veff, "b-");
    plt::xlabel("xi");
    plt::ylabel("la(xi)");
    plt::xlim(0.0, 1.0);
    plt::grid(true);

    plt::suptitle("vw = " + to_string_with_precision(vw) + ", alpha = " + to_string_with_precision(alN));
    plt::save("profile_combined.png");
}
// Kinetic power spectrum
void example_Kin_Spec(const std::string& filename) {
    // Create default universe parameters (temperature, Hubble and DoF today and at PT)
    const PhaseTransition::Universe un;

    // define PT parameters
    const auto vw = PhaseTransition::dflt_PTParams::vw;
    const auto alN = PhaseTransition::dflt_PTParams::alN;
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto dtau = PhaseTransition::dflt_PTParams::dtau;
    const auto TN = PhaseTransition::dflt_PTParams::TN;
    const auto wNeN_rat = PhaseTransition::dflt_PTParams::wNeN_rat;

    const PhaseTransition::PTParams params1(vw, alN, beta, dtau, TN, wNeN_rat, "exp", un);
    const PhaseTransition::PTParams params2(vw, alN, beta, dtau, TN, wNeN_rat, "sim", un);

    // Momentum values
    const auto kRs_vals = logspace(1e-1, 1e+3, 500);

    // Kinetic spectrum (exponential bubble nucleation)
    const auto Ek1 = Spectrum::Ekin(kRs_vals, params1);
    const auto Eks1 = Spectrum::norm_spec(Ek1); // Normalised spectrum
    // Eks1.write(filename + ".csv");

    // Kinetic spectrum (simultaneous bubble nucleation)
    const auto Ek2 = Spectrum::Ekin(kRs_vals, params2);
    const auto Eks2 = Spectrum::norm_spec(Ek2); // Normalised spectrum
    // Eks2.write(filename + ".csv");

    // Plot spectrum (alternatively, use Ek.plot())
    // plt::figure_size(800, 600);
    // plt::loglog(Eks1.K(), Eks1.P(), "k-"); // exp
    // plt::loglog(Eks2.K(), Eks2.P(), "r-"); // sim
    // plt::suptitle("vw = " + to_string_with_precision(vw) + ", alN = " + to_string_with_precision(alN));
    // plt::xlabel("K=kRs");
    // plt::ylabel("Ekin(K)");
    // plt::xlim(kRs_vals.front(), kRs_vals.back());
    // plt::ylim(1e-5, 1e+0);
    // plt::grid(true);
    // plt::save(filename + ".png");

    return;
}

// Gravitational wave power spectrum
void example_GW_Spec(const std::string& filename) {
    // Create default universe parameters (temperature, Hubble and DoF today and at PT)
    const PhaseTransition::Universe un;

    // define PT parameters
    // const auto vw = 0.8;
    // const auto alN = 0.1;
    const auto vw = PhaseTransition::dflt_PTParams::vw;
    const auto alN = PhaseTransition::dflt_PTParams::alN;
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto dtau = PhaseTransition::dflt_PTParams::dtau;
    const auto TN = PhaseTransition::dflt_PTParams::TN;
    const auto wNeN_rat = PhaseTransition::dflt_PTParams::wNeN_rat;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    const PhaseTransition::PTParams params(vw, alN, beta, dtau, TN, wNeN_rat, nuc_type, un);

    // un.print();
    // params.print();

    // Define GW spectrum
    const auto kRs_vals = logspace(1e-3, 1e+3, 100);
    const auto OmegaGW = Spectrum::GWSpec2(kRs_vals, params);
    
    // Write/plot to disk
    // OmegaGW.write(filename + ".csv");
    #ifdef ENABLE_MATPLOTLIB
    OmegaGW.plot(filename + ".png");
    #endif

    return;
}

// compare with Caprini paper
// void GW_spec_comparison(const std::string& filename) {
//     // define universe parameters
//     const auto T0 = PhaseTransition::dflt_universe::T0;
//     const auto Ts = PhaseTransition::dflt_universe::Ts;
//     const auto H0 = PhaseTransition::dflt_universe::H0;
//     const auto Hs = PhaseTransition::dflt_universe::Hs;
//     const auto g0 = PhaseTransition::dflt_universe::g0;
//     const auto gs = PhaseTransition::dflt_universe::gs;

//     const PhaseTransition::Universe un(T0, Ts, H0, Hs, g0, gs);

//     // from Pol et al.
//     const auto vw = 0.5;
//     const auto alN = 0.1;
//     const auto RsHs = 1.0;
//     const auto Rs = RsHs / Hs;
//     const auto dtau = 10.0 * Rs;

//     std::cout << "vw= " << vw << ", alN = " << alN << ", RsHs = " << RsHs << ", Rs = " << Rs << ", dtau = " << dtau << "\n";

//     // conversion to PTParams input
//     const auto betaH = std::pow(8.0 * M_PI, 1./3.) * vw / RsHs;
//     const auto beta = betaH * Hs;
    
//     const auto wNeN_rat = PhaseTransition::dflt_PTParams::wNeN_rat;
//     const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

//     const PhaseTransition::PTParams params(vw, alN, beta, dtau, wNeN_rat, nuc_type, un);

//     // un.print();
//     // params.print();

//     // Define GW spectrum
//     const auto kRs_vals = logspace(1e-3, 1e+3, 100);
//     const auto OmegaGW = Spectrum::GWSpec(kRs_vals, params);
    
//     #ifdef ENABLE_MATPLOTLIB
//     OmegaGW.plot(filename + ".png");
//     #endif

//     return;
// }

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
    // example_FluidProfile("profile");
    // example_Kin_Spec("Ekin");
    // example_GW_Spec("GW_spec");
    example_FluidProfie("profile");
    // GW_spec_comparison("GW_spec_comparison");
    // test_dSiCi_accuracy();    

    /************************ CLOCK / PROFILER *************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;
    // ProfilerStop();
    /******************************************************************/
    return 0;
}
