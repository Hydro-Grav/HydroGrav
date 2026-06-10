#include "hydrograv.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

int main() {

    // fig 1
    std::vector<double> vw_vals = {0.1, 0.3, 0.5, 0.6, 0.7, 0.9};
    for ( double vw : vw_vals ) {
        const auto alN = 0.1;
        const auto TN = 100;
        const auto beta = 100 * 9.868763e-15;
        const auto Rs = PhaseTransition::Rs_approx(vw, beta);
        // const auto cs_p = 0.573788;
        // const auto cs_m = 0.562375;

        const PhaseTransition::Universe un;

        PhaseTransition::PTParams_Bag params_exp(vw, alN, TN, beta, Rs, "exp", un, 1./3., 1./3.);
        PhaseTransition::PTParams_Bag params_sim(vw, alN, TN, beta, Rs, "sim", un, 1./3., 1./3.);

        const auto kRs_vals = logspace(-3.0, 3.0, 100);
        Spectrum::SSM_PowerSpec OmegaGW_exp = Spectrum::zetaKin(kRs_vals, params_exp);
        Spectrum::SSM_PowerSpec OmegaGW_sim = Spectrum::zetaKin(kRs_vals, params_sim);

        OmegaGW_exp.write("data/fig_1/" + std::to_string(vw) + "_exp.csv");
        OmegaGW_sim.write("data/fig_1/" + std::to_string(vw) + "_sim.csv");
    }

    //  fig 6
    std::vector<double> RsH_vals = {0.001, 0.01, 0.1, 1.0};
    for ( double vw : vw_vals ) {
        for ( double RsH : RsH_vals ) {

            const auto TN = PhaseTransition::dflt_PTParams::TN;
            const auto Hs = PhaseTransition::dflt_universe::Hs;
            const auto Rs = RsH / Hs;
            const auto beta = PhaseTransition::Rs_approx(vw, Rs);
            const auto dtau = 10.0 * Rs;

            std::cout << "vw = " << vw << ", RsH = " << RsH << ", Rs = " << Rs << ", beta = " << beta << std::endl;

            const PhaseTransition::Universe un;

            const auto kRs_vals = logspace(-3.0, 3.0, 100);

            try {
            PhaseTransition::PTParams_Bag params_01(vw, 0.1, TN, beta, Rs, "exp", un, 1./3., 1./3.);
            Spectrum::SSM_PowerSpec OmegaGW_01 = Spectrum::GWSpec(kRs_vals, params_01, dtau);
            OmegaGW_01.write("data/fig_6a/" + std::to_string(vw) + "_RsH_" + std::to_string(RsH) + ".csv");
            } catch(...) { std::cout << "Failed for fig_6a\n";}

            try {
            PhaseTransition::PTParams_Bag params_001(vw, 0.01, TN, beta, Rs, "exp", un, 1./3., 1./3.);
            Spectrum::SSM_PowerSpec OmegaGW_001 = Spectrum::GWSpec(kRs_vals, params_001, dtau);
            OmegaGW_001.write("data/fig_6b/" + std::to_string(vw) + "_RsH_" + std::to_string(RsH) + ".csv");
            } catch(...) { std::cout << "Failed for fig_6b\n";}

            try {
            PhaseTransition::PTParams_Bag params_0001(vw, 0.001, TN, beta, Rs, "exp", un, 1./3., 1./3.);
            Spectrum::SSM_PowerSpec OmegaGW_0001 = Spectrum::GWSpec(kRs_vals, params_0001, dtau);
            OmegaGW_0001.write("data/fig_6c/" + std::to_string(vw) + "_RsH_" + std::to_string(RsH) + ".csv");
            } catch(...) { std::cout << "Failed for fig_6c\n";}

            try {
                PhaseTransition::PTParams_Bag params_00001(vw, 0.0001, TN, beta, Rs, "exp", un, 1./3., 1./3.);
                Spectrum::SSM_PowerSpec OmegaGW_00001 = Spectrum::GWSpec(kRs_vals, params_00001, dtau);
                OmegaGW_00001.write("data/fig_6d/" + std::to_string(vw) + "_RsH_" + std::to_string(RsH) + ".csv");
            } catch(...) { std::cout << "Failed for fig_6d\n";}
        }
    }

    return 0;
}