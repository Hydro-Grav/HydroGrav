#include "deepphase.hpp"

int main(int argc, char* argv[]) {
    // Defaults
    auto vw = PhaseTransition::dflt_PTParams::vw;
    auto alN = PhaseTransition::dflt_PTParams::alN;
    auto beta = PhaseTransition::dflt_PTParams::beta;
    auto dtau = PhaseTransition::dflt_PTParams::dtau;
    auto TN = PhaseTransition::dflt_PTParams::TN;
    auto cpsq = PhaseTransition::dflt_PTParams::cpsq;
    auto cmsq = PhaseTransition::dflt_PTParams::cmsq;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    if (argc == 5) {
        alN = std::stod(argv[1]);
        beta = std::stod(argv[2]);
        vw = std::stod(argv[3]);
    }

    const PhaseTransition::Universe un;
    const PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, dtau, nuc_type, un, cpsq, cmsq);
    if(argc == 5) {
        params.print();
    }
    const Hydrodynamics::FluidProfile profile(params);

    profile.write("fluid_profile.csv");

    #ifdef ENABLE_MATPLOTLIB
    profile.plot("fluid_profile.png");
    #endif

    return 0;
}