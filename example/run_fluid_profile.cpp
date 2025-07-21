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

    const PhaseTransition::PTParams params(vw, alN, beta, dtau, nuc_type, un);

    // Create fluid profile
    const Hydrodynamics::FluidProfile profile(params);

    profile.write("fluid_profile.csv");

    #ifdef ENABLE_MATPLOTLIB
    profile.plot("fluid_profile.png");
    #endif

    return 0;
}