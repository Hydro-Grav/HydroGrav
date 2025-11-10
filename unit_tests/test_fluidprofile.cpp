#include "catch/catch.hpp"
#include <cmath>
#include "hydrodynamics.hpp"
#include "phasetransition.hpp"

TEST_CASE("Tests the FluidProfile class for correct behavior and integration.", "[fluidProfile]") {
    PhaseTransition::PTParams_Bag params(0.9, 0.1,
                                         PhaseTransition::dflt_PTParams::TN,
                                         PhaseTransition::dflt_PTParams::beta,
                                         PhaseTransition::dflt_PTParams::Rs,
                                         PhaseTransition::dflt_PTParams::dtau,
                                         PhaseTransition::dflt_PTParams::nuc_type,
                                         PhaseTransition::default_universe());
    Hydrodynamics::FluidProfile profile(params);

    auto xi_vals = profile.xi_vals();
    auto v_vals = profile.v_vals();
    auto w_vals = profile.w_vals();
    auto la_vals = profile.la_vals();

    for ( size_t i = 0; i < xi_vals.size(); i++) {
        const auto xi = xi_vals[i];
        const auto v = v_vals[i];
        const auto w = w_vals[i];
        const auto la = la_vals[i];

        SECTION("Check values are not NaN and within expected ranges") {
            REQUIRE(!std::isnan(xi));
            REQUIRE(!std::isnan(v));
            REQUIRE(!std::isnan(w));
            REQUIRE(!std::isnan(la));
        }

        SECTION("Check xi has expected range") {
            REQUIRE((xi >= 0.0 && xi <= 1.0));
        }

        // TODO: Add additional checks for v,w for being bubble wall and in front of shock
    }

}   
