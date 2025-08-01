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
    const auto wNeN_rat = PhaseTransition::dflt_PTParams::wNeN_rat;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    PhaseTransition::Universe un;
    
    const auto vw_vals = linspace(1e-5, 1.0 - 1e-5, 10);
    const auto alN_vals = linspace(1e-5, 1.0, 10);

    for (const auto vw : vw_vals) {
        for (const auto alN : alN_vals) {
            PhaseTransition::PTParams params(vw, alN, beta, dtau, wNeN_rat, nuc_type, un);
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

// Kinetic power spectrum
void example_Kin_Spec(const std::string& filename) {
    // Create default universe parameters (temperature, Hubble and DoF today and at PT)
    const PhaseTransition::Universe un;

    // define PT parameters
    const auto vw = PhaseTransition::dflt_PTParams::vw;
    const auto alN = PhaseTransition::dflt_PTParams::alN;
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto dtau = PhaseTransition::dflt_PTParams::dtau;
    const auto wNeN_rat = PhaseTransition::dflt_PTParams::wNeN_rat;

    const PhaseTransition::PTParams params1(vw, alN, beta, dtau, wNeN_rat, "exp", un);
    const PhaseTransition::PTParams params2(vw, alN, beta, dtau, wNeN_rat, "sim", un);

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
    const auto wNeN_rat = PhaseTransition::dflt_PTParams::wNeN_rat;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    const PhaseTransition::PTParams params(vw, alN, beta, dtau, wNeN_rat, nuc_type, un);

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
    example_GW_Spec("GW_spec");
    // test_dSiCi_accuracy();

    // const PhaseTransition::Universe un;

    // const auto vw = 0.8;
    // const auto alN = 0.1;
    // const auto beta = PhaseTransition::dflt_PTParams::beta;
    // const auto dtau = PhaseTransition::dflt_PTParams::dtau;
    // const auto wNeN_rat = PhaseTransition::dflt_PTParams::wNeN_rat;
    // const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    // const PhaseTransition::PTParams params(vw, alN, beta, dtau, wNeN_rat, nuc_type, un);
    // const Hydrodynamics::FluidProfile profile(params);

    // const auto nk = 100;
    // const auto kRs_vals = logspace(1e-3, 1e+3, nk);

    // const auto np = 1000;
    // const auto pRs_vals = logspace(1e-2, 1e+3, np);

    // const auto nz = 1000;
    // const auto z_vals = linspace(-1.0, 1.0, nz);

    // // const auto delta = Spectrum::dlt_SSM(kRs_vals, pRs_vals, z_vals, params);

    // double ptRs_min = 1.;
    // std::vector<double> kpz_min(3), kpz_max(3);
    // double ptRs_max = 2.;
    // for (size_t kk = 0; kk < nk; kk++ ) {
    //     const auto kRs = kRs_vals[kk];
    //     for (size_t pp = 0; pp < np; pp++) {
    //         const auto pRs = pRs_vals[pp];
    //         for (size_t zz = 0; zz < nz; zz++) {
    //             const auto z = z_vals[zz];
    //             const auto ptRs = Spectrum::ptilde(kRs, pRs, z);
    //             if (ptRs == 0.0)
    //                 continue;
                
    //             // ptRs_vals[kk * np * nz + pp * nz + zz] = ptRs;

    //             if (ptRs < ptRs_min) {
    //                 ptRs_min = ptRs;
    //                 kpz_min[0] = kRs;
    //                 kpz_min[1] = pRs;
    //                 kpz_min[2] = z;
    //             }

    //             if (ptRs > ptRs_max) {
    //                 ptRs_max = ptRs;
    //                 kpz_max[0] = kRs;
    //                 kpz_max[1] = pRs;
    //                 kpz_max[2] = z;
    //             }
    //         }
    //     }
    // }

    // std::cout << "ptRs_min = " << ptRs_min << ", ptRs_max = " << ptRs_max << std::endl;
    // std::cout << "kpz_min=" << kpz_min[0] << ", " << kpz_min[1] << ", " << kpz_min[2] << "\n";
    // std::cout << "kpz_max" << kpz_max[0] << ", " << kpz_max[1] << ", " << kpz_max[2] << "\n";


    // // const auto ptRs_vals_tmp = logspace(ptRs_min, ptRs_max, 2*np);
    // const auto ptRs_vals_tmp = logspace(1e-5, ptRs_max, 2*np);

    // // construct interpolating function for zetaKin(ptRs)
    // const auto zk_ptRs_spec = Spectrum::zetaKin(ptRs_vals_tmp, profile);
    // const auto zk_ptRs_K_vals = zk_ptRs_spec.K();
    // const auto zk_ptRs_P_vals = zk_ptRs_spec.P();

    // alglib::real_1d_array K_vals, P_vals;
    // K_vals.setcontent(zk_ptRs_K_vals.size(), zk_ptRs_K_vals.data());
    // P_vals.setcontent(zk_ptRs_P_vals.size(), zk_ptRs_P_vals.data());

    // alglib::spline1dinterpolant zk_ptRs_interp;;
    // alglib::spline1dbuildcubic(K_vals, P_vals, zk_ptRs_interp);

    // std::vector<double> zk_ptRs_vals(ptRs_vals_tmp.size());
    // for (size_t i = 0; i < ptRs_vals_tmp.size(); i++) {
    //     const auto ptRs = ptRs_vals_tmp[i];
    //     zk_ptRs_vals[i] = alglib::spline1dcalc(zk_ptRs_interp, ptRs);
    // }
    
    // plt::figure_size(800, 600);
    // plt::loglog(zk_ptRs_K_vals, zk_ptRs_P_vals, "k-");
    // plt::loglog(ptRs_vals_tmp, zk_ptRs_vals, "r--");
    // plt::suptitle("vw = " + to_string_with_precision(params.vw()) + ", alN = " + to_string_with_precision(params.alN()));
    // plt::xlabel("K=kRs");
    // // plt::xlim(zk_ptRs_K_vals.front(), zk_ptRs_K_vals.back());
    // plt::xlim(1e-7, 1e-4);
    // plt::grid(true);
    // plt::save("spectrum_interp.png");
    

    /************************ CLOCK / PROFILER *************************/
    const auto tf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = tf - ti;
    std::cout << "Timer: " << duration.count() << " s" << std::endl;
    // ProfilerStop();
    /******************************************************************/
    return 0;
}
