#include "deepphase.hpp"

int main(int argc, char* argv[]) {
    // Defaults
    auto vw = PhaseTransition::dflt_PTParams::vw;
    auto alN = PhaseTransition::dflt_PTParams::alN_bag;
    auto beta = PhaseTransition::dflt_PTParams::beta;
    auto Rs = PhaseTransition::dflt_PTParams::Rs;
    auto TN = PhaseTransition::dflt_PTParams::TN;
    auto cpsq = PhaseTransition::dflt_PTParams::cpsq;
    auto cmsq = PhaseTransition::dflt_PTParams::cmsq;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    if (argc == 4) {
        alN = std::stod(argv[1]);
        beta = std::stod(argv[2])*(1e-15);
        vw = std::stod(argv[3]);
    }

    // to populate the distribution, we assume T_tilde = logspace(-3, 3, 300)
    // and nu = exp(-t). This should produce a spectrum identical to the one 
    // from run_gw_spectrum.
    const auto Ttilde_values = logspace(-3.0, 3.0, 300);
    std::vector<double> nu_values;
    for (const auto& Ttilde : Ttilde_values) {
        const auto nu = exp(- Ttilde );
        nu_values.push_back(nu);
    }

    PhaseTransition::LifetimeDistribution ld_dist(Ttilde_values, nu_values);

    PhaseTransition::Universe un;
    PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, nuc_type, un, cpsq, cmsq);
    params.set_lifetime_distribution(ld_dist);
    if(argc == 5) {
        params.print();
    }

    const auto kRs_vals = logspace(-3.0, 3.0, 100);
    Spectrum::PowerSpec OmegaGW = Spectrum::GWSpec(kRs_vals, params);

    OmegaGW.write("gw_spectrum.csv");

    #ifdef ENABLE_MATPLOTLIB
    OmegaGW.plot("gw_spectrum.png");
    #endif

    return 0;
}