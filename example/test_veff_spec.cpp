#include "deepphase.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

int main() {

    const auto vw = 0.72787;
    const auto alN = 0.064721;
    const auto TN = 64.913715;
    const auto beta = 362.348213 * 6.333703e-15;
    const auto Rs =  0.006979/(8.784576e-15);
    const auto dtau = 10 * Rs;
    const auto cs_p = 0.573157;
    const auto cs_m = 0.560104;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    std::cout << "Rs * beta = " << Rs * beta << std::endl;

    // const auto vw = 0.7;
    // const auto alN = 0.1;
    // const auto TN = 100;
    // const auto beta = 100 * 9.868763e-15;
    // const auto Rs = 0.1/(9.868763e-15);
    // const auto dtau = 10 * Rs;
    // const auto cs_p = 0.573788;
    // const auto cs_m = 0.562375;
    // // const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;
    // const auto nuc_type = "exp";

    const PhaseTransition::Universe un;

    // PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, dtau, nuc_type, un, cs_p*cs_p, cs_m*cs_m);
    PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, dtau, nuc_type, un, 1./3., 1./3.);

    params.print();

    const Hydrodynamics::FluidProfile profile(params);

    std::cout << "xi_min = " << profile.xi_min() << ", xi_max = " << profile.xi_max() << std::endl;
    std::cout << "type = " << profile.mode_str() << std::endl;
    std::cout << "vw = " << params.vw() << std::endl;

    profile.write("fluid_profile.csv");

    #ifdef ENABLE_MATPLOTLIB
    profile.plot("fluid_profile.png");
    #endif

    const auto kRs_vals = logspace(1e-3, 1e+3, 100);
    Spectrum::PowerSpec OmegaGW = Spectrum::GWSpec(kRs_vals, params);

    OmegaGW.write("gw_spectrum.csv");

    #ifdef ENABLE_MATPLOTLIB
    OmegaGW.plot("gw_spectrum.png");
    #endif


    // fig 1
    // std::vector<double> vw_vals = {0.1, 0.3, 0.5, 0.6, 0.7, 0.9};
    // for ( double vw : vw_vals ) {
    //     const auto alN = 0.1;
    //     const auto TN = 100;
    //     const auto beta = 100 * 9.868763e-15;
    //     const auto Rs = PhaseTransition::Rs_approx(vw, beta);
    //     const auto dtau = 10 * Rs;
    //     const auto cs_p = 0.573788;
    //     const auto cs_m = 0.562375;

    //     const PhaseTransition::Universe un;

    //     PhaseTransition::PTParams_Bag params_exp(vw, alN, TN, beta, Rs, dtau, "exp", un, 1./3., 1./3.);
    //     PhaseTransition::PTParams_Bag params_sim(vw, alN, TN, beta, Rs, dtau, "sim", un, 1./3., 1./3.);

    //     const auto kRs_vals = logspace(1e-3, 1e+3, 100);
    //     Spectrum::PowerSpec OmegaGW_exp = Spectrum::zetaKin(kRs_vals, params_exp);
    //     Spectrum::PowerSpec OmegaGW_sim = Spectrum::zetaKin(kRs_vals, params_sim);

    //     OmegaGW_exp.write("data/" + std::to_string(vw) + "_exp.csv");
    //     OmegaGW_sim.write("data/" + std::to_string(vw) + "_sim.csv");
    // }

    // fig 6
    // std::vector<double> vw_vals = {0.1, 0.3, 0.5, 0.6, 0.7, 0.9};
    // std::vector<double> RsH_vals = {0.001, 0.01, 0.1, 1.0};

    // for ( double vw : vw_vals ) {
    //     for ( double RsH : RsH_vals ) {

    //         const auto alN = 0.1;
    //         const auto TN = PhaseTransition::dflt_PTParams::TN;
    //         const auto Hs = PhaseTransition::dflt_universe::Hs;
    //         const auto Rs = RsH / Hs;
    //         const auto beta = PhaseTransition::Rs_approx(vw, Rs);
    //         const auto dtau = 10 * Rs;

    //         std::cout << "vw = " << vw << ", RsH = " << RsH << ", Rs = " << Rs << ", beta = " << beta << std::endl;

    //         const PhaseTransition::Universe un;

    //         PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, dtau, "exp", un, 1./3., 1./3.);

    //         const auto kRs_vals = logspace(1e-3, 1e+3, 100);
    //         Spectrum::PowerSpec OmegaGW = Spectrum::GWSpec(kRs_vals, params);

    //         OmegaGW.write("data/fig_6/" + std::to_string(vw) + "_RsH_" + std::to_string(RsH) + ".csv");
    //     }
    // }

    return 0;
}