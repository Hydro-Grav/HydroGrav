#include "catch/catch.hpp"
#include "spectrum.hpp"
#include "phasetransition.hpp"

TEST_CASE("Tests functionality of the PowerSpec class.", "[powerSpec]") {
    PhaseTransition::PTParams_Bag params(0.9, 0.1); // (vw, alN)
    Hydrodynamics::FluidProfile profile(params);

    std::vector<double> K_vals = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> P_vals = {10.0, 20.0, 15.0, 25.0, 5.0};

    Spectrum::PowerSpec spec(K_vals, P_vals, profile);

    SECTION("Test max") {
        REQUIRE(spec.max() == Approx(25.0).epsilon(1e-10));
    }

    // there's no test here...?
    // SECTION("Test interpolation") {
    //     auto interp = spec.interpolate();
    //     double test_k = 2.5;
    //     double interp_val = interp(test_k);
    // }
    
    SECTION("Test operations") {
        Spectrum::PowerSpec scaled = spec * 2.0;
        REQUIRE(scaled.P()[1] == Approx(40.0).epsilon(1e-10));

        spec *= 0.5;
        REQUIRE(spec.P()[1] == Approx(10.0).epsilon(1e-10));

        Spectrum::PowerSpec divided = scaled / 2.0;
        REQUIRE(divided.P()[1] == Approx(20.0).epsilon(1e-10));

        scaled /= 4.0;
        REQUIRE(scaled.P()[1] == Approx(10.0).epsilon(1e-10));
    }


}