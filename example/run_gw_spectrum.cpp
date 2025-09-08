#include "deepphase.hpp"

int main(int argc, char* argv[]) {
    // Defaults
    auto vw = PhaseTransition::dflt_PTParams::vw;
    auto alN = PhaseTransition::dflt_PTParams::alN;
    auto betaH = PhaseTransition::dflt_PTParams::betaH;
    auto dtauH = PhaseTransition::dflt_PTParams::dtauH;
    auto wNeN_rat = PhaseTransition::dflt_PTParams::wNeN_rat;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    if (argc == 6) {
        alN = std::stod(argv[1]);
        betaH = std::stod(argv[2]);
        dtauH = std::stod(argv[3]);
        vw = std::stod(argv[4]);
        wNeN_rat = std::stod(argv[5]);
    }

    const PhaseTransition::Universe un;
    const PhaseTransition::PTParams params(vw, alN, betaH, dtauH, wNeN_rat, nuc_type, un);
    params.print();

    const auto kRs_vals = logspace(1e-3, 1e+2, 100);
    Spectrum::PowerSpec OmegaGW = Spectrum::GWSpec2(kRs_vals, params);

    OmegaGW.write("gw_spectrum.csv");

    #ifdef ENABLE_MATPLOTLIB
    OmegaGW.plot("gw_spectrum.png");
    #endif

    return 0;
}