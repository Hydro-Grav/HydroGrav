#include "deepphase.hpp"

int main() {

    /*
        This example script generalises the run_fluid_profile to work
        for arbitrary equation of state. The generalised equation of state
        is stored in the EquationOfState object.
    */

    const auto vw = PhaseTransition::dflt_PTParams::vw;
    const auto alN = PhaseTransition::dflt_PTParams::alN_bag;
    const auto TN = PhaseTransition::dflt_PTParams::TN;
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto Rs = PhaseTransition::dflt_PTParams::Rs;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    const PhaseTransition::Universe un;

    /*
        Mock equation of state data based on Bag model 
        fit to above thermal params. In practice, this data
        would come from a dedicated calculation of the
        finite temperature effective potential.
    */
    std::vector<double> T, ps, pb, es, eb;
    const double Tmin = 20;
    const double Tmax = 100;
    const int n_points = 300;
    const double a_s = M_PI*M_PI/30.0 * un.gs();
    const double a_b = 0.9 * a_s;
    const double epsilon = alN/(1-alN) * a_s * std::pow(TN,4.);
    for (double tt = Tmin; tt <= Tmax; tt += (Tmax - Tmin)/(n_points-1)) {
        double tt4 = tt*tt*tt*tt;
        T.push_back(tt);
        ps.push_back(a_s/3.0 * tt4 - epsilon);
        pb.push_back(a_b/3.0 * tt4);
        es.push_back(a_s * tt4 + epsilon);
        eb.push_back(a_b * tt4);
    }

    /*
        Firstly, we create the EquationOfState object from the eos data
    */
    PhaseTransition::EquationOfState eos_data(T, ps, pb, es, eb);

    /* 
        Secondly, we construct the PTParams object
    */
    PhaseTransition::PTParams_Veff params(vw, alN, TN, beta, Rs, nuc_type, un, eos_data);

    /*
        Lastly, we use this to compute the fluid profile
        Note as we use mock data above, this results in an error.
    */
    const Hydrodynamics::FluidProfile profile_from_eos(params);

    profile_from_eos.write("fluid_profile_eos_from_vectors.csv");

    #ifdef ENABLE_MATPLOTLIB
    profile_from_eos.plot("fluid_profile_eos_from_vectors.png");
    #endif


    /*
        Alternatively, the EoS can be loaded from a precomputed file. The format should be:
        T,pb,ps,eb,es
        10.0,0.8,1.0,2.4,3.0
        20.0,3.2,4.0,9.6,12.0
        30.0,7.2,9.0,21.6,27.0
        ...
        We do this using the from_file function in EquationOfState.
    */

    const std::string eos_path = "thermo.csv";

    PhaseTransition::EquationOfState eos_data_from_path = PhaseTransition::EquationOfState::from_file(eos_path);

    /*
        This can then be passed to PTParams_Veff. Additionally, PTParams_Veff itself can be 
        initialised from the file path (for backwards compatibility, might not make it to release)
    */
    PhaseTransition::PTParams_Veff params_from_path(vw, alN, TN, beta, Rs, nuc_type, un, eos_path);

    const Hydrodynamics::FluidProfile profile_from_path(params_from_path);

    profile_from_path.write("fluid_profile_eos_from_file.csv");

    #ifdef ENABLE_MATPLOTLIB
    profile_from_path.plot("fluid_profile_eos_from_file.png");
    #endif

    return 0;
}