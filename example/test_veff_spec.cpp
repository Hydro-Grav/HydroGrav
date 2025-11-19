#include "deepphase.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

int main() {

    const auto vw = 0.661453;
    const auto alN = 0.030803;
    const auto TN = 76.118057;
    const auto beta = 668.119813 * 8.485469e-15;
    const auto Rs = 0.003697/(8.485469e-15);
    const auto dtau = 10 * Rs;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    const PhaseTransition::Universe un;

    // std::vector<double> T, ps, pb, es, eb;

    // // Load data from CSV file
    // std::ifstream file("eos.csv");
    // if (!file.is_open()) {
    //     throw std::runtime_error("Failed to open CSV file: eos.csv");
    // }

    // std::string line;
    // std::getline(file, line);

    // while (std::getline(file, line)) {
    //     std::stringstream ss(line);
    //     std::string token;
    //     std::vector<double> values;

    //     while (std::getline(ss, token, ',')) {
    //         values.push_back(std::stod(token));
    //     }

    //     if (values.size() >= 5) {
    //         T.push_back(values[0]);   // T
    //         ps.push_back(values[1]);  // P_plus
    //         pb.push_back(values[2]);  // P_minus
    //         es.push_back(values[3]);  // e_plus
    //         eb.push_back(values[4]);  // e_minus
    //     }
    // }
    // file.close();

    // PhaseTransition::EquationOfState eos_data(T, ps, pb, es, eb);

    // PhaseTransition::PTParams_Veff params(vw, alN, TN, beta, Rs, dtau, nuc_type, un, eos_data);
    PhaseTransition::PTParams_Bag params(vw, alN, TN, beta, Rs, dtau, nuc_type, un, 1.0/3.0, 1.0/3.0);

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