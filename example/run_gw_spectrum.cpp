#include "deepphase.hpp"

int main(int argc, char* argv[]) {
    // Create default universe parameters (temperature, Hubble and DoF today and at PT)
    const PhaseTransition::Universe un;

    // define PT parameters
    double alN, beta, vw;
    if (argc == 4) {
        alN = std::stod(argv[1]);
        beta = std::stod(argv[2])*(1e-15);
        vw = std::stod(argv[3]);
    } else { // dlft values
        alN = PhaseTransition::dflt_PTParams::alN_bag;
        beta = PhaseTransition::dflt_PTParams::beta;
        vw = PhaseTransition::dflt_PTParams::vw;
    }

    auto Rs = PhaseTransition::dflt_PTParams::Rs;
    auto TN = PhaseTransition::dflt_PTParams::TN;
    auto cpsq = PhaseTransition::dflt_PTParams::cpsq;
    auto cmsq = PhaseTransition::dflt_PTParams::cmsq;
    auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    const PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, nuc_type, un, cpsq, cmsq);

    // Momentum values
    const auto kRs_vals = logspace(-3.0, 3.0, 100);

    // GW power spectrum
    Spectrum::PowerSpec OmegaGW = Spectrum::GWSpec(kRs_vals, params);

    OmegaGW.write("gw_spectrum.csv");

    #ifdef ENABLE_MATPLOTLIB
    OmegaGW.plot("gw_spectrum.png");
    #endif

    return 0;
}