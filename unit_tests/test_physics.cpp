#include "catch/catch.hpp"
#include <cmath>
#include <vector>
#include "physics.hpp"
#include "maths_ops.hpp"

// ============================================================
// gammaSq()
// ============================================================
TEST_CASE("gammaSq()", "[gammaSq]") {
    SECTION("v=0 gives 1") {
        REQUIRE(gammaSq(0.0) == Approx(1.0).epsilon(1e-12));
    }
    SECTION("v=0.5 gives 4/3") {
        REQUIRE(gammaSq(0.5) == Approx(4.0 / 3.0).epsilon(1e-12));
    }
    SECTION("v=1/sqrt(2) gives 2") {
        REQUIRE(gammaSq(1.0 / std::sqrt(2.0)) == Approx(2.0).epsilon(1e-10));
    }
    SECTION("increases monotonically with v") {
        double prev = gammaSq(0.0);
        for (double v = 0.1; v < 0.99; v += 0.1) {
            double curr = gammaSq(v);
            REQUIRE(curr > prev);
            prev = curr;
        }
    }
    SECTION("always >= 1 for 0 <= v < 1") {
        for (double v = 0.0; v < 1.0; v += 0.05)
            REQUIRE(gammaSq(v) >= 1.0);
    }
}

// ============================================================
// P_oms() and P_acc()
// ============================================================
TEST_CASE("P_oms()", "[lisa_noise]") {
    SECTION("positive for positive frequencies") {
        for (double f : {1e-4, 1e-3, 1e-2, 1e-1, 1.0})
            REQUIRE(P_oms(f) > 0.0);
    }
    SECTION("increases at low frequency due to 1/f^4 term") {
        // At very low f the (2e-3/f)^4 term dominates — P_oms should be larger at lower f
        REQUIRE(P_oms(1e-5) > P_oms(1e-2));
    }
    SECTION("baseline value at high frequency") {
        // At large f the correction term vanishes: P_oms -> (1.5e-11)^2
        double expected = (1.5e-11) * (1.5e-11);
        REQUIRE(P_oms(1.0) == Approx(expected).epsilon(1e-3));
    }
}

TEST_CASE("P_acc()", "[lisa_noise]") {
    SECTION("positive for positive frequencies") {
        for (double f : {1e-4, 1e-3, 1e-2, 1e-1})
            REQUIRE(P_acc(f) > 0.0);
    }
    SECTION("increases at low frequency (1/f^4 envelope)") {
        REQUIRE(P_acc(1e-4) > P_acc(1e-3));
        REQUIRE(P_acc(1e-3) > P_acc(1e-2));
    }
}

// ============================================================
// get_LISA_omegahsq()
// ============================================================
TEST_CASE("get_LISA_omegahsq()", "[lisa_sensitivity]") {
    SECTION("positive at LISA band frequencies") {
        for (double f : {1e-4, 1e-3, 1e-2, 1e-1})
            REQUIRE(get_LISA_omegahsq(f) > 0.0);
    }
    SECTION("has a minimum in the LISA band (around 3-10 mHz)") {
        // The LISA sensitivity curve has a bucket shape — minimum is near 3 mHz
        double val_low  = get_LISA_omegahsq(1e-4);
        double val_mid  = get_LISA_omegahsq(3e-3);
        double val_high = get_LISA_omegahsq(1e-1);
        REQUIRE(val_mid < val_low);
        REQUIRE(val_mid < val_high);
    }
    SECTION("proportional to f^3 * noise at high frequency") {
        // Increasing f should eventually cause Omega to increase
        double f1 = 0.05, f2 = 0.1;
        REQUIRE(get_LISA_omegahsq(f2) > get_LISA_omegahsq(f1));
    }
}

// ============================================================
// LISA_snr()
// ============================================================
TEST_CASE("LISA_snr()", "[lisa_snr]") {
    // Build a simple flat signal well above the noise
    const size_t N = 200;
    auto freqs = linspace(1e-4, 1e-1, N);

    SECTION("SNR is positive for a signal above noise") {
        // Signal 100x the noise floor at each frequency
        std::vector<double> amp(N);
        for (size_t i = 0; i < N; ++i)
            amp[i] = 100.0 * get_LISA_omegahsq(freqs[i]);
        double snr = LISA_snr(freqs, amp);
        REQUIRE(snr > 0.0);
    }

    SECTION("SNR scales with observation time") {
        std::vector<double> amp(N);
        for (size_t i = 0; i < N; ++i)
            amp[i] = 10.0 * get_LISA_omegahsq(freqs[i]);
        double snr_4yr  = LISA_snr(freqs, amp, 4.0);
        double snr_8yr  = LISA_snr(freqs, amp, 8.0);
        // SNR ~ sqrt(T), so doubling T should multiply SNR by sqrt(2)
        REQUIRE(snr_8yr == Approx(snr_4yr * std::sqrt(2.0)).epsilon(1e-10));
    }

    SECTION("larger signal gives larger SNR") {
        std::vector<double> amp_low(N), amp_high(N);
        for (size_t i = 0; i < N; ++i) {
            amp_low[i]  = 1.0  * get_LISA_omegahsq(freqs[i]);
            amp_high[i] = 10.0 * get_LISA_omegahsq(freqs[i]);
        }
        REQUIRE(LISA_snr(freqs, amp_high) > LISA_snr(freqs, amp_low));
    }

    SECTION("throws on size mismatch") {
        std::vector<double> amp(N - 1, 1e-10);
        REQUIRE_THROWS_AS(LISA_snr(freqs, amp), std::invalid_argument);
    }

    SECTION("throws on empty input") {
        REQUIRE_THROWS_AS(LISA_snr({}, {}), std::invalid_argument);
    }
}
