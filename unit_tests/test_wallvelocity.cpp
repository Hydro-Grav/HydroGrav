#include "catch/catch.hpp"

#include <cmath>
#include <memory>
#include <vector>

#include "logger.hpp"
#include "phasetransition.hpp"
#include "profile.hpp"
#include "wallvelocity.hpp"

namespace {

/*
 * Template (bag) EoS with c_s^2 = c_b^2 = 1/3:
 *   p_s = a_s T^4/3 - eps,  e_s = a_s T^4 + eps,  p_b = a_b T^4/3,  e_b = a_b T^4,
 * with a_b = psiN a_s and eps = alN a_s TN^4, so that the template-model
 * strength alN = (theta_s - theta_b)/(3 w_s(TN)) with theta = e - 3p, and
 * psiN = w_b(TN)/w_s(TN).
 */
PhaseTransition::EquationOfState template_eos(double alN, double psiN, double TN, double smoothing = 0.0) {
    const double a_s = M_PI * M_PI / 30.0 * 106.75;
    const double a_b = psiN * a_s;
    const double eps = alN * a_s * std::pow(TN, 4);

    std::vector<double> T, ps, pb, es, eb;
    const int n = 400;
    for (int i = 0; i < n; ++i) {
        const double t = TN * (0.5 + 1.5 * i / (n - 1.0));
        const double t4 = std::pow(t, 4);
        T.push_back(t);
        ps.push_back(a_s / 3.0 * t4 - eps);
        pb.push_back(a_b / 3.0 * t4);
        es.push_back(a_s * t4 + eps);
        eb.push_back(a_b * t4);
    }

    PhaseTransition::EquationOfState eos(T, ps, pb, es, eb);
    eos.smoothing = smoothing;
    return eos;
}

std::unique_ptr<PhaseTransition::PTParams_Veff> template_params(double alN, double psiN, double smoothing = 0.0) {
    const double TN = PhaseTransition::dflt_PTParams::TN;
    const double beta = PhaseTransition::dflt_PTParams::beta;
    const double vw = 0.5; // placeholder, ignored by the LTE solver
    return std::make_unique<PhaseTransition::PTParams_Veff>(
        vw, alN, TN, beta,
        PhaseTransition::Rs_approx(vw, beta), PhaseTransition::Rbar_approx(vw, beta),
        "exp", PhaseTransition::default_universe(), template_eos(alN, psiN, TN, smoothing));
}

// Analytic Jouguet velocity of the bag model (c_b^2 = 1/3)
double vJ_bag(double alN) {
    return (std::sqrt(2.0 * alN / 3.0 + alN * alN) + std::sqrt(1.0 / 3.0)) / (1.0 + alN);
}

Hydrodynamics::LTEOptions fast_options() {
    return Hydrodynamics::LTEOptions{}; // defaults from config.hpp
}

} // namespace

TEST_CASE("PTParams_Veff template-model diagnostics", "[lte][ptparams]") {
    logging::ScopedLevel quiet(logging::Level::Error);

    const auto p = template_params(0.05, 0.9);

    SECTION("alN_munu and psiN reproduce the EoS inputs") {
        REQUIRE(p->alN_munu() == Approx(0.05).epsilon(1e-5));
        REQUIRE(p->psiN() == Approx(0.9).epsilon(1e-6));
    }

    SECTION("exact EoS is thermodynamically consistent") {
        for (double TTN : {0.7, 1.0, 1.3, 1.7}) {
            REQUIRE(std::abs(p->entropy_consistency(TTN, false)) < 1e-3);
            REQUIRE(std::abs(p->entropy_consistency(TTN, true)) < 1e-3);
        }
    }

    SECTION("default spline smoothing is flagged as inconsistent") {
        const auto smoothed = template_params(0.05, 0.9, 1e-2);
        const Hydrodynamics::LTEWallSolver lte(*smoothed);
        REQUIRE(lte.max_entropy_inconsistency() > config::lte_entropy_consistency_tol);
    }

    SECTION("with_vw sets vw and rescales Rs and Rbar") {
        const auto q = p->with_vw(0.3);
        const double beta = p->beta();
        REQUIRE(q->vw() == Approx(0.3));
        REQUIRE(q->Rs() == Approx(PhaseTransition::Rs_approx(0.3, beta)));
        REQUIRE(q->Rbar() == Approx(PhaseTransition::Rbar_approx(0.3, beta)));
        REQUIRE(q->ps_val(1.1) == Approx(p->ps_val(1.1)));
        REQUIRE(p->vw() == Approx(0.5)); // original untouched

        const auto r = p->with_vw(0.3, false);
        REQUIRE(r->Rs() == Approx(p->Rs()));
        REQUIRE(r->Rbar() == Approx(p->Rbar()));

        REQUIRE_THROWS_AS(p->with_vw(1.2), std::invalid_argument);
    }
}

TEST_CASE("jouguet_velocity_veff matches the bag model", "[lte]") {
    logging::ScopedLevel quiet(logging::Level::Error);
    for (double alN : {0.02, 0.05, 0.1}) {
        const auto p = template_params(alN, 0.9);
        REQUIRE(Hydrodynamics::jouguet_velocity_veff(*p) == Approx(vJ_bag(alN)).epsilon(1e-4));
    }
}

TEST_CASE("LTE wall velocity matches Ai, Laurent & van de Vis", "[lte]") {
    logging::ScopedLevel quiet(logging::Level::Error);

    /*
     * Reference values from the Python code in Appendix C of arXiv:2303.10171
     * (find_vw with cb2 = cs2 = 1/3), which solves the template-model matching
     * conditions analytically and shoots on the shock condition: an
     * independent formulation from the tabulated-EoS + entropy residual used here.
     */
    struct Ref { double alN, psiN, vw; int mode; };
    const std::vector<Ref> refs = {
        {0.04, 0.90, 0.51356792, 0}, // deflagration (vw < c_b)
        {0.05, 0.90, 0.59741829, 1}, // hybrid
        {0.10, 0.85, 0.65312869, 1}, // hybrid
    };

    for (const auto& ref : refs) {
        INFO("alN=" << ref.alN << ", psiN=" << ref.psiN);
        const auto p = template_params(ref.alN, ref.psiN);
        const Hydrodynamics::LTEWallSolver lte(*p, fast_options());
        const auto res = lte.solve();

        REQUIRE(res.status == Hydrodynamics::LTEStatus::Converged);
        REQUIRE(res.vw == Approx(ref.vw).margin(1e-3));
        REQUIRE(res.mode == ref.mode);
        REQUIRE(res.vw < res.vJ);

        // The LTE condition holds at the solution, to the precision of the inner
        // T-/TN shooting in FluidProfile: Delta(vw) is piecewise smooth with
        // jumps of O(1e-5), so the root can sit on one of them.
        REQUIRE(std::abs(lte.entropy_residual(res.vw)) < 1e-4);

        // returned parameters carry the solved wall velocity
        REQUIRE(res.params);
        REQUIRE(res.params->vw() == Approx(res.vw));
        REQUIRE(res.params->Rs() == Approx(PhaseTransition::Rs_approx(res.vw, p->beta())));
    }
}

TEST_CASE("LTE wall velocity classifies runaway and non-expanding walls", "[lte]") {
    logging::ScopedLevel quiet(logging::Level::Error);

    SECTION("alN below alpha_min = (1 - psiN)/3: no expansion") {
        const auto p = template_params(0.02, 0.9);
        const auto res = Hydrodynamics::LTEWallSolver(*p, fast_options()).solve();
        REQUIRE(res.status == Hydrodynamics::LTEStatus::NoExpansion);
        REQUIRE_FALSE(res.params);
    }

    SECTION("alN above alpha_max^hyb: runaway") {
        // reference code: alN_max^hyb(psiN=0.9) = 0.123
        const auto p = template_params(0.3, 0.9);
        const auto res = Hydrodynamics::LTEWallSolver(*p, fast_options()).solve();
        REQUIRE(res.status == Hydrodynamics::LTEStatus::Runaway);
        REQUIRE_FALSE(res.params);
    }
}
