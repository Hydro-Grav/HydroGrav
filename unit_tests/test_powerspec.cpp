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

// ============================================================
// SoundShellKernel / dlt_SSM
// ============================================================

// Direct Si/Ci evaluation of the sound-shell time-correlation function. This is the
// definition dlt_SSM() must reproduce; the fast path in SoundShellKernel replaces it with a
// series in eps = dtau/(tau_s + tau_fin).
static double dlt_SSM_reference(double k, double p, double pt, double cs,
                                double tau_s, double tau_fin) {
    const double ptcs = pt * cs;
    const double pcs = p * cs;

    double dlt = 0.0;
    for (int m = -1; m < 2; m += 2) {
        const double pmn_1 = pcs + m * ptcs;
        for (int n = -1; n < 2; n += 2) {
            const double pmn = pmn_1 + n * k;

            double Si_1, Ci_1, Si_2, Ci_2;
            alglib::sinecosineintegrals(pmn * tau_fin, Si_1, Ci_1);
            alglib::sinecosineintegrals(pmn * tau_s, Si_2, Ci_2);

            const double dSi = Si_1 - Si_2;
            const double dCi = Ci_1 - Ci_2;
            dlt += 0.25 * (dCi * dCi + dSi * dSi);
        }
    }
    return dlt;
}

TEST_CASE("SoundShellKernel reproduces the Si/Ci definition", "[ssm_kernel]") {
    // tau_s = 1/Hs is ~900 Rs for typical parameters, so eps is small and the series path
    // is the one that matters; the largest ratio here sits at the edge of its validity.
    const double tau_s = 911.7;
    const double cs = 0.56705;

    SECTION("dlt_SSM matches the reference over a wide (k, p, z) range") {
        for (double dtau_over_tau_s : {0.010975, 0.05, 0.3}) {
            const double tau_fin = tau_s * (1.0 + dtau_over_tau_s);
            const Spectrum::SoundShellKernel kernel(tau_s, tau_fin);
            REQUIRE_FALSE(kernel.uses_exact());

            const double wk_scale = kernel.half_dtau();

            double worst = 0.0, peak = 0.0;
            for (int i = 0; i < 40; ++i) {
                const double k = std::pow(10.0, -3.0 + 6.0 * i / 39.0);
                const double sin_wk = std::sin(k * wk_scale);
                const double cos_wk = std::cos(k * wk_scale);
                for (int j = 0; j < 40; ++j) {
                    const double p = std::pow(10.0, -3.0 + 6.0 * j / 39.0);
                    for (int l = 0; l <= 20; ++l) {
                        const double z = -1.0 + 2.0 * l / 20.0;
                        const double pt = Spectrum::ptilde(k, p, z);

                        const double ref = dlt_SSM_reference(k, p, pt, cs, tau_s, tau_fin);
                        const double fast =
                            Spectrum::dlt_SSM(k, p, pt, cs, kernel, sin_wk, cos_wk);

                        peak = std::max(peak, ref);
                        worst = std::max(worst, std::abs(ref - fast));
                    }
                }
            }
            // Error is O(eps^5) relative to the peak of Delta.
            REQUIRE(peak > 0.0);
            REQUIRE(worst / peak < 1e-5);
        }
    }

    SECTION("the tau_s/tau_fin overload agrees with the kernel overload") {
        const double tau_fin = tau_s * 1.010975;
        const Spectrum::SoundShellKernel kernel(tau_s, tau_fin);
        const double wk_scale = kernel.half_dtau();

        for (double k : {1e-2, 1.0, 7.5, 250.0}) {
            const double p = 0.6 * k + 0.3;
            const double pt = Spectrum::ptilde(k, p, 0.25);
            const double a = Spectrum::dlt_SSM(k, p, pt, cs, tau_s, tau_fin);
            const double b = Spectrum::dlt_SSM(k, p, pt, cs, kernel,
                                               std::sin(k * wk_scale), std::cos(k * wk_scale));
            REQUIRE(a == Approx(b).epsilon(1e-12));
        }
    }

    SECTION("falls back to the exact path for long-lived sound waves") {
        // dtau comparable to tau_s pushes eps past the series' range of validity.
        const Spectrum::SoundShellKernel kernel(tau_s, 3.0 * tau_s);
        REQUIRE(kernel.uses_exact());

        for (double k : {1e-2, 1.0, 30.0}) {
            const double p = 0.6 * k + 0.3;
            const double pt = Spectrum::ptilde(k, p, -0.4);
            const double ref = dlt_SSM_reference(k, p, pt, cs, tau_s, 3.0 * tau_s);
            const double got = Spectrum::dlt_SSM(k, p, pt, cs, kernel, 0.0, 1.0);
            // Not bit-identical: the two group the product cs*(p + m*pt) differently, and
            // the Ci/Si differences amplify that last-bit difference.
            REQUIRE(got == Approx(ref).epsilon(1e-8));
        }
    }

    SECTION("G is even in its argument") {
        const Spectrum::SoundShellKernel kernel(tau_s, tau_s * 1.010975);
        for (double x : {1e-6, 1e-3, 0.5, 3.0, 40.0, 5000.0}) {
            REQUIRE(kernel(x) == Approx(kernel(-x)).epsilon(1e-12));
        }
    }
}

// ============================================================
// GWSpec regression
// ============================================================
TEST_CASE("GWSpec regression", "[gwspec]") {
    // Values pinned against a deliberately over-resolved evaluation of the same integral
    // (a dense uniform grid in ptRs, ~8 points per resonance half-width, 6000 log-p points).
    // The quadrature in GWSpec targets ~1e-3 relative accuracy, so compare at that level.
    //
    // NOTE: the momentum grid is part of the fixture, not a convenience. gw_prefac()
    // normalises by the peak of Ekin() *over the grid it is given*, so a coarser grid shifts
    // the overall amplitude.
    using D = PhaseTransition::dflt_PTParams;
    const PhaseTransition::Universe un;
    const PhaseTransition::PTParams_Bag params(D::vw, D::alN_bag, D::TN, D::beta, D::Rs,
                                               D::nuc_type, un, D::cpsq, D::cmsq);

    const auto kRs_vals = logspace(-3.0, 3.0, 100);
    const auto gw = Spectrum::GWSpec(kRs_vals, params, 10.0 * D::Rs);

    REQUIRE(gw.P().size() == kRs_vals.size());
    REQUIRE(gw.dtau() == Approx(10.0 * D::Rs).epsilon(1e-12));

    // {index into the grid, over-resolved value}
    const std::vector<std::pair<size_t, double>> expected = {
        {33, 3.517690674e-18},   // kRs = 0.1
        {50, 2.989474464e-16},   // kRs = 1.07227
        {66, 2.051501273e-15},   // kRs = 10
        {83, 2.670272866e-16},   // kRs = 107.227
        {99, 1.395283186e-19},   // kRs = 1000
    };

    for (const auto& [i, want] : expected) {
        INFO("kRs = " << kRs_vals[i]);
        CHECK(gw.P()[i] == Approx(want).epsilon(2e-3));
    }

    SECTION("spectrum is positive and finite everywhere") {
        for (double p : gw.P()) {
            REQUIRE(p > 0.0);
            REQUIRE(std::isfinite(p));
        }
    }
}
