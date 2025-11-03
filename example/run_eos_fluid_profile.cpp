#include "deepphase.hpp"

int main() {

    /*
        This example script generalises the run_fluid_profile to work
        for arbitrary equation of state. The generalised equation of state
        is stored in the EquationOfState object.
    */

    const auto vw = PhaseTransition::dflt_PTParams::vw;
    const auto alN = PhaseTransition::dflt_PTParams::alN;
    const auto TN = PhaseTransition::dflt_PTParams::TN;
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto dtau = PhaseTransition::dflt_PTParams::dtau;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    const PhaseTransition::Universe un;
    const std::string eos_path = "thermo.csv";

    /*
        Firstly, we create the EquationOfState object from the eos data
    */
    std::vector<double> T  = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0};
    std::vector<double> ps = {1.0,  4.0,  9.0,  16.0, 25.0, 36.0};
    std::vector<double> pb = {0.8,  3.2,  7.2,  12.8, 20.0, 28.8};
    std::vector<double> es = {3.0,  12.0, 27.0, 48.0, 75.0, 108.0};
    std::vector<double> eb = {2.4,  9.6,  21.6, 38.4, 60.0, 86.4};

    PhaseTransition::EquationOfState eos_data(T, ps, pb, es, eb);

    /* 
        Secondly, we construct the PTParams object
    */
    PhaseTransition::PTParams_Veff params(vw, alN, TN, beta, dtau, nuc_type, un, eos_data);

    /*
        Lastly, we use this to compute the fluid profile
        Note as we use mock data above, this results in an error.
    */
    try 
    { 
        const Hydrodynamics::FluidProfile profile(params);
    } catch (...) {
        std::cout << "!!! Expected error caught and handled. Carry on !!!\n";
    }

    /*
        Alternatively, the EoS can be loaded from a precomputed file. The format should be:
        T,pb,ps,eb,es
        10.0,0.8,1.0,2.4,3.0
        20.0,3.2,4.0,9.6,12.0
        30.0,7.2,9.0,21.6,27.0
        ...
        We do this using the from_file function in EquationOfState.
    */

    PhaseTransition::EquationOfState eos_data_from_path = PhaseTransition::EquationOfState::from_file(eos_path);

    /*
        This can then be passed to PTParams_Veff. Additionally, PTParams_Veff itself can be 
        initialised from the file path (for backwards compatibility, might not make it to release)
    */
    PhaseTransition::PTParams_Veff params_from_path(vw, alN, TN, beta, dtau, nuc_type, un, eos_path);

    const Hydrodynamics::FluidProfile profile_from_path(params_from_path);

    profile_from_path.write("fluid_profile_eos.csv");

    #ifdef ENABLE_MATPLOTLIB
    profile_from_path.plot("fluid_profile_eos.png");
    #endif

    return 0;
}