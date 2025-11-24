#include "deepphase.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

int main() {

    // const auto vw = 0.577499;
    // const auto alN = 0.010780;
    // const auto TN = 94.049935;
    // const auto beta = 1794.730595 * 9.868763e-15;
    // const auto Rs = 0.001409/(9.868763e-15);
    // const auto dtau = 10 * Rs;
    // const auto cs_p = 0.573788;
    // const auto cs_m = 0.562375;
    // const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    const auto vw = 0.7;
    const auto alN = 0.1;
    const auto TN = 100;
    const auto beta = 100 * 9.868763e-15;
    const auto Rs = 0.1/(9.868763e-15);
    const auto dtau = 10 * Rs;
    const auto cs_p = 0.573788;
    const auto cs_m = 0.562375;
    // const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;
    const auto nuc_type = "exp";

    const PhaseTransition::Universe un;

    PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, dtau, nuc_type, un, 1./3., 1./3.);

    const Hydrodynamics::FluidProfile profile(params);

    profile.write("fluid_profile.csv");

    #ifdef ENABLE_MATPLOTLIB
    profile.plot("fluid_profile.png");
    #endif

    const auto kRs_vals = logspace(1e-3, 1e+3, 100);
    for(double kRs : kRs_vals) { std::cout << "k = " << kRs << "\n";}
    Spectrum::PowerSpec OmegaGW = Spectrum::GWSpec(kRs_vals, params);

    OmegaGW.write("gw_spectrum.csv");

    #ifdef ENABLE_MATPLOTLIB
    OmegaGW.plot("gw_spectrum.png");
    #endif

    // std::vector<double> vw_vals = {0.1, 0.3, 0.5, 0.6, 0.7, 0.9};

    // for ( double vw : vw_vals ) {
    //     const auto alN = 0.1;
    //     const auto TN = 100;
    //     const auto beta = 100 * 9.868763e-15;
    //     const auto Rs = 0.1/(9.868763e-15);
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

    return 0;
}