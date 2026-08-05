#include "hydrograv.hpp"

int main(int argc, char* argv[]) {
    // Create default universe parameters (temperature, Hubble and DoF today and at PT)
    const PhaseTransition::Universe un;

    // define PT parameters
    // double alN, beta, vw;
    // if (argc == 4) {
    //     alN = std::stod(argv[1]);
    //     beta = std::stod(argv[2])*(1e-15);
    //     vw = std::stod(argv[3]);
    // } else { // dlft values
    //     alN = PhaseTransition::dflt_PTParams::alN_bag;
    //     beta = PhaseTransition::dflt_PTParams::beta;
    //     vw = PhaseTransition::dflt_PTParams::vw;
    // }

    // BPs from Fig. 8 arXiv:2406.01596
    // deflagration
    // double alN = -0.12015;
    // double vw = 0.75;

    // detonation
    double alN = -0.164931; // approx
    double vw = 0.25;

    double cpsq = 1.0 / 3.0;
    double cmsq = cpsq;
    
    auto beta = PhaseTransition::dflt_PTParams::beta;
    auto Rs = PhaseTransition::dflt_PTParams::Rs;
    auto TN = PhaseTransition::dflt_PTParams::TN;
    auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    const PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, nuc_type, un, cpsq, cmsq);
    params.print();

    const Hydrodynamics::FluidProfile profile(params, 5000, true);
    profile.write();

    // profile.write("fluid_profile.csv");

    #ifdef ENABLE_MATPLOTLIB
    profile.plot("fluid_profile.png", 0.0, 0.8);
    #endif
    

    return 0;
}