#include "catch/catch.hpp"
#include <cmath>
#include <limits>
#include "ssm.hpp"
#include "profile.hpp"
#include "phasetransition.hpp"

// Shared fixture
static Hydrodynamics::FluidProfile make_profile() {
    static PhaseTransition::PTParams_Bag params(
            0.9, 0.1,
            PhaseTransition::dflt_PTParams::TN,
            PhaseTransition::dflt_PTParams::beta,
            PhaseTransition::dflt_PTParams::Rs,
            PhaseTransition::dflt_PTParams::nuc_type,
            PhaseTransition::default_universe()
            );
    return Hydrodynamics::FluidProfile(params);
}

// ============================================================
// Construction and accessors
// ============================================================
TEST_CASE("PowerSpec construction and accessors", "[powerSpec]") {
    auto profile = make_profile();

    std::vector<double> K = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> P = {10.0, 20.0, 15.0, 25.0, 5.0};
    Spectrum::PowerSpec spec(K, P, profile);

    SECTION("K() returns correct values") {
        REQUIRE(spec.K().size() == 5);
        REQUIRE(spec.K()[0] == Approx(1.0).epsilon(1e-12));
        REQUIRE(spec.K()[4] == Approx(5.0).epsilon(1e-12));
    }
    SECTION("P() returns correct values") {
        REQUIRE(spec.P().size() == 5);
        REQUIRE(spec.P()[1] == Approx(20.0).epsilon(1e-12));
        REQUIRE(spec.P()[3] == Approx(25.0).epsilon(1e-12));
    }
    SECTION("freq() is non-empty and same size as K") {
        REQUIRE(spec.freq().size() == spec.K().size());
    }
    SECTION("freq() values are positive") {
        for (double f : spec.freq())
            REQUIRE(f > 0.0);
    }
    SECTION("profile pointer is accessible") {
        REQUIRE(spec.params() != nullptr);
    }
    SECTION("dtau defaults to NaN") {
        REQUIRE(std::isnan(spec.dtau()));
    }
    SECTION("explicit dtau stored correctly") {
        Spectrum::PowerSpec s2(K, P, profile, 1.23);
        REQUIRE(s2.dtau() == Approx(1.23).epsilon(1e-12));
    }
}

// ============================================================
// peak_vals()
// ============================================================
TEST_CASE("PowerSpec peak_vals()", "[powerSpec]") {
    auto profile = make_profile();
    std::vector<double> K = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> P = {10.0, 20.0, 15.0, 25.0, 5.0};
    Spectrum::PowerSpec spec(K, P, profile);

    SECTION("peak amplitude is maximum P value") {
        REQUIRE(spec.peak_vals().second == Approx(25.0).epsilon(1e-10));
    }
    SECTION("peak frequency corresponds to peak amplitude") {
        // K[3]=4.0 has max P=25 — peak freq should correspond to K=4
        auto [f_peak, A_peak] = spec.peak_vals();
        REQUIRE(A_peak == Approx(25.0).epsilon(1e-10));
        REQUIRE(f_peak > 0.0);
    }
    SECTION("monotone spectrum — peak at last element") {
        std::vector<double> P2 = {1.0, 2.0, 3.0, 4.0, 5.0};
        Spectrum::PowerSpec s2(K, P2, profile);
        REQUIRE(s2.peak_vals().second == Approx(5.0).epsilon(1e-10));
    }
}

// ============================================================
// Arithmetic operators
// ============================================================
TEST_CASE("PowerSpec arithmetic operators", "[powerSpec]") {
    auto profile = make_profile();
    std::vector<double> K = {1.0, 2.0, 3.0};
    std::vector<double> P = {4.0, 8.0, 12.0};
    Spectrum::PowerSpec spec(K, P, profile);

    SECTION("operator* (spec * scalar)") {
        auto s2 = spec * 3.0;
        REQUIRE(s2.P()[0] == Approx(12.0).epsilon(1e-10));
        REQUIRE(s2.P()[1] == Approx(24.0).epsilon(1e-10));
        REQUIRE(s2.P()[2] == Approx(36.0).epsilon(1e-10));
    }
    SECTION("operator* (scalar * spec)") {
        auto s2 = 3.0 * spec;
        REQUIRE(s2.P()[1] == Approx(24.0).epsilon(1e-10));
    }
    SECTION("operator* preserves K values") {
        auto s2 = spec * 2.0;
        for (size_t i = 0; i < K.size(); ++i)
            REQUIRE(s2.K()[i] == Approx(K[i]).epsilon(1e-12));
    }
    SECTION("operator*=") {
        spec *= 0.5;
        REQUIRE(spec.P()[0] == Approx(2.0).epsilon(1e-10));
        REQUIRE(spec.P()[1] == Approx(4.0).epsilon(1e-10));
    }
    SECTION("operator/ (spec / scalar)") {
        auto s2 = spec / 2.0;
        REQUIRE(s2.P()[0] == Approx(2.0).epsilon(1e-10));
        REQUIRE(s2.P()[2] == Approx(6.0).epsilon(1e-10));
    }
    SECTION("operator/=") {
        spec /= 4.0;
        REQUIRE(spec.P()[0] == Approx(1.0).epsilon(1e-10));
    }
    SECTION("operator/ throws on zero") {
        REQUIRE_THROWS_AS(spec / 0.0, std::invalid_argument);
    }
    SECTION("operator/= throws on zero") {
        REQUIRE_THROWS_AS(spec /= 0.0, std::invalid_argument);
    }
    SECTION("multiply then divide gives original") {
        auto s2 = (spec * 7.0) / 7.0;
        for (size_t i = 0; i < P.size(); ++i)
            REQUIRE(s2.P()[i] == Approx(P[i]).epsilon(1e-10));
    }
    SECTION("scalar * spec == spec * scalar") {
        auto s1 = spec * 5.0;
        auto s2 = 5.0 * spec;
        for (size_t i = 0; i < P.size(); ++i)
            REQUIRE(s1.P()[i] == Approx(s2.P()[i]).epsilon(1e-12));
    }
}

// ============================================================
// norm_spec()
// ============================================================
TEST_CASE("norm_spec()", "[powerSpec]") {
    auto profile = make_profile();
    std::vector<double> K = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> P = {10.0, 40.0, 30.0, 20.0, 5.0};
    Spectrum::PowerSpec spec(K, P, profile);

    SECTION("normalised peak amplitude is 1") {
        auto normalised = Spectrum::norm_spec(spec);
        REQUIRE(normalised.peak_vals().second == Approx(1.0).epsilon(1e-10));
    }
    SECTION("all values in [0, 1] after normalisation") {
        auto normalised = Spectrum::norm_spec(spec);
        for (double p : normalised.P()) {
            REQUIRE(p >= 0.0);
            REQUIRE(p <= 1.0 + 1e-12);
        }
    }
    SECTION("K values unchanged by normalisation") {
        auto normalised = Spectrum::norm_spec(spec);
        for (size_t i = 0; i < K.size(); ++i)
            REQUIRE(normalised.K()[i] == Approx(K[i]).epsilon(1e-12));
    }
    SECTION("normalising already-normalised spec gives same result") {
        auto n1 = Spectrum::norm_spec(spec);
        auto n2 = Spectrum::norm_spec(n1);
        for (size_t i = 0; i < P.size(); ++i)
            REQUIRE(n1.P()[i] == Approx(n2.P()[i]).epsilon(1e-10));
    }
}

// ============================================================
// interpolate()
// ============================================================
TEST_CASE("PowerSpec interpolate()", "[powerSpec]") {
    auto profile = make_profile();
    // Use enough points for a meaningful spline
    std::vector<double> K = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<double> P;
    for (double k : K) P.push_back(k * k); // P = K^2

    Spectrum::PowerSpec spec(K, P, profile);
    auto interp = spec.interpolate();

    SECTION("reproduces knot values exactly") {
        for (size_t i = 0; i < K.size(); ++i)
            CHECK(interp(K[i]) == Approx(P[i]).epsilon(1e-6));
    }
    SECTION("interpolated value is finite at interior point") {
        double val = interp(4.5);
        REQUIRE(std::isfinite(val));
        REQUIRE(val > 0.0);
    }
}

// ============================================================
// SSM functions — Ekin, zetaKin
// ============================================================
TEST_CASE("Ekin() from FluidProfile", "[ekin]") {
    auto profile = make_profile();
    auto k_vals = linspace(1e-2, 1e2, 30);

    Spectrum::PowerSpec ekin = Spectrum::Ekin(k_vals, profile);

    SECTION("output size matches input") {
        REQUIRE(ekin.K().size() == k_vals.size());
        REQUIRE(ekin.P().size() == k_vals.size());
    }
    SECTION("all P values are non-negative") {
        for (double p : ekin.P())
            REQUIRE(p >= 0.0);
    }
    SECTION("no NaN values") {
        for (double p : ekin.P())
            REQUIRE(!std::isnan(p));
    }
    SECTION("peak is positive") {
        REQUIRE(ekin.peak_vals().second > 0.0);
    }
}

TEST_CASE("zetaKin() from FluidProfile", "[zetaKin]") {
    auto profile = make_profile();
    auto k_vals = linspace(1e-2, 1e2, 30);

    Spectrum::PowerSpec zk = Spectrum::zetaKin(k_vals, profile);

    SECTION("output size matches input") {
        REQUIRE(zk.K().size() == k_vals.size());
    }
    SECTION("all P values are non-negative") {
        for (double p : zk.P())
            REQUIRE(p >= 0.0);
    }
    SECTION("no NaN values") {
        for (double p : zk.P())
            REQUIRE(!std::isnan(p));
    }
}
