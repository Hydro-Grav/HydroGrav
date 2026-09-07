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

    double cpsq = 1.0 / 3.0;
    double cmsq = cpsq;

    /********** BPs from Fig. 8 arXiv:2406.01596 **********/
    // Note: they use slightly different value for ap/am in detonation/hybrid profiles

    // detonation (vp=0.5, vm=0.25, alp=-0.277778)
    // double alN = -0.164931; // approx - they only give alp
    // double vw = 0.25;

    // hybrid (vp=0.55, vm=0.9, alp=-0.100173)
    // double alN = 
    // double vw = 

    // deflagration (vp=0.75, vm=0.89, alp=-0.12015)
    // double alN = -0.12015;
    // double vw = 0.75;
    /******************************************************/
    /**************** BPs from Xiao's code ****************/
    double alN = -0.1;
    // double vw = 0.3; // detonation
    double vw = 0.5; // hybrid
    // double vw = 0.8; // deflagration
    /******************************************************/

    
    
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
    profile.plot("fp.png", 0.0, 1.0);
    #endif
    

    return 0;
}