#include "deepphase.hpp"

int main(int argc, char* argv[]) {
    // Defaults
    auto vw = PhaseTransition::dflt_PTParams::vw;
    auto alN = PhaseTransition::dflt_PTParams::alN;
    auto beta = PhaseTransition::dflt_PTParams::beta;
    auto dt = PhaseTransition::dflt_PTParams::dt;
    auto wNeN_rat = PhaseTransition::dflt_PTParams::wNeN_rat;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    if (argc == 5) {
        alN = std::stod(argv[1]);
        beta = std::stod(argv[2]);
        vw = std::stod(argv[3]);
        wNeN_rat = std::stod(argv[4]);
        dt = 1/beta;
    }

    const PhaseTransition::Universe un;
    const PhaseTransition::PTParams params(vw, alN, beta, dt, wNeN_rat, nuc_type, un);
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