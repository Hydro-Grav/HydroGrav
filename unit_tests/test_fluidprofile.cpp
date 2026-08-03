#include "catch/catch.hpp"
#include <cmath>
#include "profile.hpp"
#include "phasetransition.hpp"

// Shared fixture: default bag parameters from the paper benchmark point
static PhaseTransition::PTParams_Bag make_bag_params(double vw = 0.9, double alN = 0.1) {
    return PhaseTransition::PTParams_Bag(
        vw, alN,
        PhaseTransition::dflt_PTParams::TN,
        PhaseTransition::dflt_PTParams::beta,
        PhaseTransition::dflt_PTParams::Rs,
        PhaseTransition::dflt_PTParams::nuc_type,
        PhaseTransition::default_universe()
    );
}

// ============================================================
// mu() — Lorentz boost between frames
// ============================================================
TEST_CASE("mu() — Lorentz boost", "[mu]") {
    SECTION("mu(xi, xi) = 0 (identical velocities)") {
        REQUIRE(Hydrodynamics::mu(0.5, 0.5) == Approx(0.0).margin(1e-12));
    }
    SECTION("mu(xi, 0) = xi (wall at rest in UF)") {
        REQUIRE(Hydrodynamics::mu(0.7, 0.0) == Approx(0.7).epsilon(1e-12));
    }
    SECTION("result lies in (-1, 1) for physical velocities") {
        for (double xi = 0.1; xi < 1.0; xi += 0.1)
            for (double v = 0.05; v < 1.0; v += 0.1)
                if (xi != v) {
                    double m = Hydrodynamics::mu(xi, v);
                    REQUIRE(std::abs(m) < 1.0);
                }
    }
    SECTION("anti-symmetric: mu(xi,v) = -mu(v,xi)") {
        REQUIRE(Hydrodynamics::mu(0.6, 0.3) == Approx(-Hydrodynamics::mu(0.3, 0.6)).epsilon(1e-12));
    }
}

// ============================================================
// dxidv() — fluid characteristic slope
// ============================================================
TEST_CASE("dxidv()", "[dxidv]") {
    const double csq = 1.0 / 3.0; // bag EoS sound speed squared
    SECTION("finite for typical values") {
        double d = Hydrodynamics::dxidv(0.5, 0.2, csq);
        REQUIRE(std::isfinite(d));
    }
    SECTION("sign is positive for xi > v (outgoing rarefaction)") {
        // dxi/dv should be positive when the wall moves faster than the fluid
        double d = Hydrodynamics::dxidv(0.8, 0.1, csq);
        REQUIRE(d > 0.0);
    }
}

// ============================================================
// PTParams construction and accessors
// ============================================================
TEST_CASE("PTParams_Bag construction", "[ptparams]") {
    auto params = make_bag_params(0.9, 0.1);

    SECTION("vw stored correctly") {
        REQUIRE(params.vw() == Approx(0.9).epsilon(1e-12));
    }
    SECTION("alN stored correctly") {
        REQUIRE(params.alN() == Approx(0.1).epsilon(1e-12));
    }
    SECTION("default bag cpsq = 1/3") {
        REQUIRE(params.cpsq() == Approx(1.0/3.0).epsilon(1e-12));
    }
    SECTION("default bag cmsq = 1/3") {
        REQUIRE(params.cmsq() == Approx(1.0/3.0).epsilon(1e-12));
    }
    SECTION("EoS type is Bag") {
        REQUIRE(params.eos() == PhaseTransition::PTParams::ModelType::Bag);
    }
    SECTION("eos_to_string returns Bag") {
        REQUIRE(params.eos_to_string() == "Bag");
    }
    SECTION("Rs is positive") {
        REQUIRE(params.Rs() > 0.0);
    }
    SECTION("beta is positive") {
        REQUIRE(params.beta() > 0.0);
    }
    SECTION("tau_s is positive") {
        REQUIRE(params.tau_s() > 0.0);
    }
}

TEST_CASE("PTParams_Bag with explicit sound speeds", "[ptparams]") {
    PhaseTransition::PTParams_Bag params(
        0.7, 0.05,
        PhaseTransition::dflt_PTParams::TN,
        PhaseTransition::dflt_PTParams::beta,
        PhaseTransition::dflt_PTParams::Rs,
        PhaseTransition::dflt_PTParams::nuc_type,
        PhaseTransition::default_universe(),
        0.4, 0.3
    );
    REQUIRE(params.cpsq() == Approx(0.4).epsilon(1e-12));
    REQUIRE(params.cmsq() == Approx(0.3).epsilon(1e-12));
}

TEST_CASE("Rs_approx()", "[ptparams]") {
    double vw = 0.5, beta = 1e-13;
    double Rs = PhaseTransition::Rs_approx(vw, beta);
    REQUIRE(Rs > 0.0);
    // Rs ~ (8pi)^(1/3) * vw / beta
    double expected = std::pow(8.0 * M_PI, 1.0/3.0) * vw / beta;
    REQUIRE(Rs == Approx(expected).epsilon(1e-12));
}

// ============================================================
// Universe construction
// ============================================================
TEST_CASE("Universe construction and accessors", "[universe]") {
    const auto& un = PhaseTransition::default_universe();

    SECTION("temperatures are positive") {
        REQUIRE(un.T0() > 0.0);
        REQUIRE(un.Ts() > 0.0);
    }
    SECTION("Hubble constants are positive") {
        REQUIRE(un.H0() > 0.0);
        REQUIRE(un.Hs() > 0.0);
    }
    SECTION("dof are positive") {
        REQUIRE(un.g0() > 0.0);
        REQUIRE(un.gs() > 0.0);
    }
    SECTION("Ts >> T0 (universe was hotter at PT)") {
        REQUIRE(un.Ts() > un.T0());
    }
}

// ============================================================
// FluidProfile — detonation (vw=0.9, large alN)
// ============================================================
TEST_CASE("FluidProfile bag detonation (vw=0.9, alN=0.1)", "[fluidProfile]") {
    auto params = make_bag_params(0.9, 0.1);
    Hydrodynamics::FluidProfile profile(params);

    SECTION("mode is detonation (2)") {
        REQUIRE(profile.mode() == 2);
        REQUIRE(profile.mode_str() == "detonation");
    }

    SECTION("profile arrays are non-empty and equal in length") {
        REQUIRE(!profile.xi_vals().empty());
        REQUIRE(profile.v_vals().size()  == profile.xi_vals().size());
        REQUIRE(profile.w_vals().size()  == profile.xi_vals().size());
        REQUIRE(profile.la_vals().size() == profile.xi_vals().size());
        REQUIRE(profile.T_vals().size()  == profile.xi_vals().size());
    }

    SECTION("xi values lie in (0, 1)") {
        for (double xi : profile.xi_vals()) {
            REQUIRE(xi >= 0.0);
            REQUIRE(xi <= 1.0);
        }
    }

    SECTION("fluid velocities lie in [0, 1)") {
        for (double v : profile.v_vals()) {
            REQUIRE(v >= 0.0);
            REQUIRE(v < 1.0);
        }
    }

    SECTION("enthalpy ratios are positive") {
        for (double w : profile.w_vals())
            REQUIRE(w > 0.0);
    }

    SECTION("temperature ratios are positive") {
        for (double T : profile.T_vals())
            REQUIRE(T > 0.0);
    }

    SECTION("no NaN values in any profile") {
        for (double x : profile.xi_vals()) REQUIRE(!std::isnan(x));
        for (double v : profile.v_vals())  REQUIRE(!std::isnan(v));
        for (double w : profile.w_vals())  REQUIRE(!std::isnan(w));
        for (double T : profile.T_vals())  REQUIRE(!std::isnan(T));
        for (double l : profile.la_vals()) REQUIRE(!std::isnan(l));
    }

    SECTION("xi_min < vw (rarefaction starts inside wall)") {
        REQUIRE(profile.xi_min() <= params.vw());
    }

    SECTION("xi_max > vw (shock is ahead of wall)") {
        REQUIRE(profile.xi_max() >= params.vw());
    }

    SECTION("profile pointer matches params") {
        REQUIRE(profile.params() == &params);
    }
}

// ============================================================
// FluidProfile — deflagration (vw=0.4, alN=0.05)
// ============================================================
TEST_CASE("FluidProfile bag deflagration (vw=0.4, alN=0.05)", "[fluidProfile]") {
    auto params = make_bag_params(0.4, 0.05);
    Hydrodynamics::FluidProfile profile(params);

    SECTION("mode is deflagration (0)") {
        REQUIRE(profile.mode() == 0);
        REQUIRE(profile.mode_str() == "deflagration");
    }

    SECTION("profile arrays are non-empty") {
        REQUIRE(!profile.xi_vals().empty());
    }

    SECTION("no NaN values") {
        for (double x : profile.xi_vals()) REQUIRE(!std::isnan(x));
        for (double v : profile.v_vals())  REQUIRE(!std::isnan(v));
        for (double w : profile.w_vals())  REQUIRE(!std::isnan(w));
    }

    SECTION("xi values lie in (0, 1)") {
        for (double xi : profile.xi_vals()) {
            REQUIRE(xi >= 0.0);
            REQUIRE(xi <= 1.0);
        }
    }
}

// ============================================================
// FluidProfile — benchmark bag values from paper
// ============================================================
TEST_CASE("FluidProfile bag benchmark (dflt_PTParams)", "[fluidProfile][benchmark]") {
    PhaseTransition::PTParams_Bag params(
        PhaseTransition::dflt_PTParams::vw,
        PhaseTransition::dflt_PTParams::alN_bag,
        PhaseTransition::dflt_PTParams::TN,
        PhaseTransition::dflt_PTParams::beta,
        PhaseTransition::dflt_PTParams::Rs,
        PhaseTransition::dflt_PTParams::nuc_type,
        PhaseTransition::default_universe()
    );
    Hydrodynamics::FluidProfile profile(params);

    SECTION("profile is non-empty") {
        REQUIRE(!profile.xi_vals().empty());
    }
    SECTION("no NaNs") {
        for (double v : profile.v_vals()) REQUIRE(!std::isnan(v));
        for (double w : profile.w_vals()) REQUIRE(!std::isnan(w));
    }
    SECTION("enthalpy equals 1 far from the wall (unperturbed plasma)") {
        // The very last point (far from wall) should have w/wN ~ 1
        const auto& w = profile.w_vals();
        REQUIRE(w.back() == Approx(1.0).epsilon(0.05));
    }
}
