#include "deepphase.hpp"

int main() {
    // Create default universe parameters (temperature, Hubble and DoF today and at PT)
    const PhaseTransition::Universe un;

    // Define PT parameters
    const auto vw = PhaseTransition::dflt_PTParams::vw;
    const auto alN = PhaseTransition::dflt_PTParams::alN;
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto dtau = PhaseTransition::dflt_PTParams::dtau;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    if (argc == 4) {;
        alN = std::stod(argv[1]);
        beta = std::stod(argv[2]);
        vw = std::stod(argv[3]);
        dtau = 1/beta;
    }

    const PhaseTransition::PTParams params(vw, alN, beta, dtau, nuc_type, un);
    if(argc == 4) {
        params.print();
    }

    // Construct GW power spectrum
    const auto kRs_vals = logspace(1e-3, 1e+3, 100);
    Spectrum::PowerSpec OmegaGW = Spectrum::GWSpec(kRs_vals, params);

    OmegaGW.write("gw_spectrum.csv");

    #ifdef ENABLE_MATPLOTLIB
    OmegaGW.plot("gw_spectrum.png");
    #endif

    return 0;
}