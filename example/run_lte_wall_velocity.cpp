#include <cstdlib>
#include <iostream>

#include "hydrograv.hpp"

/*
    This example computes the bubble-wall velocity in local thermal
    equilibrium (LTE) from a tabulated equation of state, then uses it to
    build the fluid profile and gravitational wave spectrum.

    In LTE, entropy is conserved across the wall, which adds the junction
    condition gamma+ T+ = gamma- T- to the usual energy-momentum conditions.
    That extra condition fixes the wall velocity
    (Ai, Laurent & van de Vis, arXiv:2303.10171).

    Usage:  ./bin/run_lte_wall_velocity [alN] [psiN]
      alN   template-model strength parameter at T_N     (default 0.05)
      psiN  enthalpy ratio w_b(T_N)/w_s(T_N)              (default 0.90)
*/

int main(int argc, char* argv[]) {
    const double alN  = (argc > 1) ? std::atof(argv[1]) : 0.05;
    const double psiN = (argc > 2) ? std::atof(argv[2]) : 0.90;

    const PhaseTransition::Universe un;
    const auto TN = PhaseTransition::dflt_PTParams::TN;
    const auto beta = PhaseTransition::dflt_PTParams::beta;
    const auto nuc_type = PhaseTransition::dflt_PTParams::nuc_type;

    /*
        Mock equation of state: the template (bag) model with c_s^2 = c_b^2 = 1/3,
            p_s = a_s T^4/3 - eps,  e_s = a_s T^4 + eps,
            p_b = a_b T^4/3,        e_b = a_b T^4,
        with a_b = psiN * a_s and eps = alN * a_s * TN^4, so that
        alN = (theta_s - theta_b) / (3 w_s(TN)) with theta = e - 3p.
        In practice this table would come from the finite-temperature
        effective potential of a particle physics model.
    */
    const double a_s = M_PI * M_PI / 30.0 * un.gs();
    const double a_b = psiN * a_s;
    const double eps = alN * a_s * std::pow(TN, 4);

    std::vector<double> T, ps, pb, es, eb;
    const int n_points = 400;
    const double Tmin = 0.5 * TN, Tmax = 2.0 * TN;
    for (int i = 0; i < n_points; ++i) {
        const double tt = Tmin + (Tmax - Tmin) * i / (n_points - 1);
        const double tt4 = std::pow(tt, 4);
        T.push_back(tt);
        ps.push_back(a_s / 3.0 * tt4 - eps);
        pb.push_back(a_b / 3.0 * tt4);
        es.push_back(a_s * tt4 + eps);
        eb.push_back(a_b * tt4);
    }
    PhaseTransition::EquationOfState eos(T, ps, pb, es, eb);

    /*
        This table is exact, so fit it without smoothing. The default smoothing
        (meant for noisy tables from an effective potential) biases p(T) and
        breaks w = T dp/dT, which the LTE condition relies on.
    */
    eos.smoothing = 0.0;

    /*
        The wall velocity passed here is a placeholder: the solver ignores it.
        Rs and Rbar are rescaled to the solved wall velocity.
    */
    const double vw_placeholder = 0.5;
    const auto Rs = PhaseTransition::Rs_approx(vw_placeholder, beta);
    const auto Rbar = PhaseTransition::Rbar_approx(vw_placeholder, beta);
    const PhaseTransition::PTParams_Veff base(vw_placeholder, alN, TN, beta, Rs, Rbar, nuc_type, un, eos);

    /*
        Solve for the LTE wall velocity. The scan of the entropy residual
        Delta(vw) = gamma+ T+ / (gamma- T-) - 1 is written to disk for inspection.
    */
    const Hydrodynamics::LTEWallSolver lte(base);
    const auto res = lte.solve();
    Hydrodynamics::LTEWallSolver::write_scan(res, "lte_scan.csv");

    std::cout << "\n  alN = " << alN << ", psiN = " << psiN << "\n"
              << "  status = " << Hydrodynamics::to_string(res.status) << "\n"
              << "  vJ     = " << res.vJ << "\n";

    if (res.status != Hydrodynamics::LTEStatus::Converged) {
        return 0;
    }

    std::cout << "  vw     = " << res.vw << " (mode " << res.mode << ")\n"
              << "  T+/TN = " << res.wall.TpTN << ", T-/TN = " << res.wall.TmTN
              << ", v+ = " << res.wall.vp << ", v- = " << res.wall.vm << "\n\n";

    /*
        res.params carries the solved wall velocity and must outlive anything
        built from it (FluidProfile and PowerSpec keep a pointer to it).
    */
    const Hydrodynamics::FluidProfile profile(*res.params);
    profile.write("fluid_profile_lte.csv");

    const auto kRs_vals = logspace(-3.0, 3.0, 100);
    const auto dtau = 10.0 * res.params->Rs();
    const auto OmegaGW = Spectrum::GWSpec(kRs_vals, *res.params, dtau);
    OmegaGW.write("gw_spectrum_lte.csv");

    #ifdef ENABLE_MATPLOTLIB
    profile.plot("fluid_profile_lte.png");
    OmegaGW.plot("gw_spectrum_lte.png");
    #endif

    return 0;
}
