#include "hydrograv.hpp"
#include <random>

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
    auto Rbar = PhaseTransition::dflt_PTParams::Rbar; // new dflt value for Rbar!
    auto TN = PhaseTransition::dflt_PTParams::TN;
    auto cpsq = PhaseTransition::dflt_PTParams::cpsq;
    auto cmsq = PhaseTransition::dflt_PTParams::cmsq;
    auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    // new PTParams ctor takes in Rbar!
    const PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, Rbar, nuc_type, un, cpsq, cmsq);

    // Momentum values
    const auto kRs_vals = logspace(-3.0, 3.0, 100);

    // construct LogNormal(0,1) distro to sample from
    std::random_device rd;
    std::mt19937 gen(rd());
    std::lognormal_distribution<> d(0.0, 1.0);

    // GW power spectrum
    Spectrum::PowerSpec OmegaGW = Spectrum::sample_GWSpec(kRs_vals, params, d(gen));
    std::cout << "dtau/Rs=" << OmegaGW.dtau() / OmegaGW.params()->Rs() << "\n";

    auto snr_list = get_SNR(OmegaGW.freq(), OmegaGW.P());
    for (const auto& snr_result : snr_list) {
        std::cout << snr_result.detector_name << " SNR: " << snr_result.snr << std::endl;
    }

    OmegaGW.write("gw_spectrum.csv");

    #ifdef ENABLE_MATPLOTLIB
    OmegaGW.profile().plot("fp.png");
    OmegaGW.plot("gw_spectrum.png");
    #endif

    return 0;
}