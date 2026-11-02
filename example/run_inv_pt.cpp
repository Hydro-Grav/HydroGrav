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

    /********** BPs from Fig. 8/9 arXiv:2406.01596 **********/
    // Note: they use slightly different value for ap/am in detonation/hybrid profiles
    // vp, vm, alp for hybrid are approx. values from Fig. 9)

    // detonation (vp=0.5, vm=0.25, alp=-0.277778)
    // double alN = -0.164931; // approx - they only give alp
    // double vw = 0.25;

    // hybrid (vp=0.577, vm=0.95, alp=-0.128)
    // double alN = -0.1;
    // double vw = 0.5;

    // deflagration (vp=0.75, vm=0.89, alp=-0.12015)
    // double alN = -0.12015;
    // double vw = 0.75;
    /******************************************************/
    /**************** BPs from Xiao's code ****************/
    // double alN = -0.1;
    // double vw = 0.3; // detonation
    // double vw = 0.5; // hybrid
    // double vw = 0.8; // deflagration
    /******************************************************/    
    
    auto beta = PhaseTransition::dflt_PTParams::beta;
    auto Rs = PhaseTransition::dflt_PTParams::Rs;
    auto TN = PhaseTransition::dflt_PTParams::TN;
    auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    const std::vector<double> vw_list = {0.25, 0.5, 0.75};
    const std::vector<double> alN_list = {-0.164931, -0.1, -0.12015};
    const std::vector<std::string> filenames = {"fp_inv_det", "fp_inv_hyb", "fp_inv_def"};

    for (int i = 0; i < vw_list.size(); i++) {
        const auto vw = vw_list[i];
        const auto alN = alN_list[i];
        const auto filename = filenames[i];
        
        const PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, nuc_type, un, cpsq, cmsq);
        const Hydrodynamics::FluidProfile profile(params, 5000, true);
        profile.write(filename + ".csv");
        profile.plot(filename + ".png");
    }

    // const PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, nuc_type, un, cpsq, cmsq);
    // params.print();

    // const Hydrodynamics::FluidProfile profile(params, 5000, true);
    // profile.write("fp_inv.csv");

    // #ifdef ENABLE_MATPLOTLIB
    // profile.plot("fp_inv.png");
    // #endif
    

    return 0;
}