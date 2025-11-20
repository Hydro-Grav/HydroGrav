#include "deepphase.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

int main() {

    const auto vw = 0.577499;
    const auto alN = 0.010780;
    const auto TN = 94.049935;
    const auto beta = 1794.730595 * 1.198829e-14;
    const auto Rs = 0.001409/(9.974751e-15);
    const auto dtau = 10 * Rs;
    const auto cs_p = 0.573788;
    const auto cs_m = 0.562375;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    const PhaseTransition::Universe un;

    PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, dtau, nuc_type, un, 1./3., 1./3.);

    const Hydrodynamics::FluidProfile profile(params);

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

    return 0;
}