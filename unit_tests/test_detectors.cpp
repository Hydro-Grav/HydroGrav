#include "catch/catch.hpp"
#include <cmath>
#include <vector>
#include <stdexcept>

#include "detectors.hpp"
#include "snr.hpp"
#include "maths.hpp"

// ============================================================
// Detector::name() and ndet()
// ============================================================
TEST_CASE("Detector names", "[detectors]") {
    REQUIRE(LISA().name() == "LISA");
    REQUIRE(DECIGO().name() == "DECIGO");
}

TEST_CASE("Detector ndet() defaults and overrides", "[detectors]") {
    SECTION("default ndet is 1") {
        REQUIRE(LISA().ndet() == Approx(1.0).epsilon(1e-12));
    }
    SECTION("DECIGO overrides ndet to 2") {
        REQUIRE(DECIGO().ndet() == Approx(2.0).epsilon(1e-12));
    }
}

// ============================================================
// omega_noise() — positivity and basic shape for each detector
// ============================================================
TEST_CASE("LISA omega_noise()", "[detectors][lisa]") {
    LISA lisa;
    SECTION("positive at LISA band frequencies") {
        for (double f : {1e-4, 1e-3, 1e-2, 1e-1})
            REQUIRE(lisa.omega_noise(f) > 0.0);
    }
    SECTION("has a minimum in the LISA band (bucket shape)") {
        double val_low  = lisa.omega_noise(1e-4);
        double val_mid  = lisa.omega_noise(3e-3);
        double val_high = lisa.omega_noise(1e-1);
        REQUIRE(val_mid < val_low);
        REQUIRE(val_mid < val_high);
    }
    SECTION("increases with frequency above the bucket minimum") {
        REQUIRE(lisa.omega_noise(0.1) > lisa.omega_noise(0.05));
    }
}

TEST_CASE("DECIGO omega_noise()", "[detectors][decigo]") {
    DECIGO decigo;
    SECTION("positive across the DECIGO band") {
        for (double f : {1e-2, 1e-1, 1.0, 7.36, 10.0})
            REQUIRE(decigo.omega_noise(f) > 0.0);
    }
    SECTION("finite near the characteristic frequency f_p = 7.36 Hz") {
        REQUIRE(std::isfinite(decigo.omega_noise(7.36)));
    }
}

// ============================================================
// calculate_snr() — generic single-detector SNR
// ============================================================
TEST_CASE("calculate_snr() with LISA", "[snr]") {
    const size_t N = 200;
    auto freqs = linspace(1e-4, 1e-1, N);
    LISA lisa;

    SECTION("positive for a signal well above the noise floor") {
        std::vector<double> amp(N);
        for (size_t i = 0; i < N; ++i)
            amp[i] = 100.0 * lisa.omega_noise(freqs[i]);
        REQUIRE(calculate_snr(freqs, amp, lisa) > 0.0);
    }

    SECTION("SNR scales as sqrt(Tyear)") {
        std::vector<double> amp(N);
        for (size_t i = 0; i < N; ++i)
            amp[i] = 10.0 * lisa.omega_noise(freqs[i]);
        double snr_4yr = calculate_snr(freqs, amp, lisa, 4.0);
        double snr_8yr = calculate_snr(freqs, amp, lisa, 8.0);
        REQUIRE(snr_8yr == Approx(snr_4yr * std::sqrt(2.0)).epsilon(1e-10));
    }

    SECTION("larger signal gives larger SNR") {
        std::vector<double> amp_low(N), amp_high(N);
        for (size_t i = 0; i < N; ++i) {
            amp_low[i]  = 1.0  * lisa.omega_noise(freqs[i]);
            amp_high[i] = 10.0 * lisa.omega_noise(freqs[i]);
        }
        REQUIRE(calculate_snr(freqs, amp_high, lisa) > calculate_snr(freqs, amp_low, lisa));
    }

    SECTION("throws on size mismatch") {
        std::vector<double> amp(N - 1, 1e-10);
        REQUIRE_THROWS_AS(calculate_snr(freqs, amp, lisa), std::invalid_argument);
    }

    SECTION("throws on empty input") {
        REQUIRE_THROWS_AS(calculate_snr({}, {}, lisa), std::invalid_argument);
    }
}

// ============================================================
// calculate_all_snrs()
// ============================================================
TEST_CASE("calculate_all_snrs()", "[snr]") {
    const size_t N = 100;
    auto freqs = linspace(1e-3, 1e-1, N);

    static const LISA lisa;
    static const DECIGO decigo;
    
    std::vector<double> amp(N);
    for (size_t i = 0; i < N; ++i)
        amp[i] = 50.0 * lisa.omega_noise(freqs[i]);

    SECTION("returns one result per detector, in order") {
        std::vector<const Detector*> detectors = { &lisa, &decigo };
        auto results = calculate_all_snrs(freqs, amp, detectors);

        REQUIRE(results.size() == 2);
        REQUIRE(results[0].detector_name == "LISA");
        REQUIRE(results[1].detector_name == "DECIGO");
    }

    SECTION("each SNR matches a direct calculate_snr() call") {
        std::vector<const Detector*> detectors = { &lisa, &decigo };
        auto results = calculate_all_snrs(freqs, amp, detectors);

        REQUIRE(results[0].snr == Approx(calculate_snr(freqs, amp, lisa)).epsilon(1e-10));
        REQUIRE(results[1].snr == Approx(calculate_snr(freqs, amp, decigo)).epsilon(1e-10));
    }

    SECTION("empty detector list returns empty result") {
        std::vector<const Detector*> detectors;
        auto results = calculate_all_snrs(freqs, amp, detectors);
        REQUIRE(results.empty());
    }
}

// ============================================================
// get_SNR() — all four standard detectors
// ============================================================
TEST_CASE("get_SNR() returns LISA, DECIGO in order", "[snr]") {
    const size_t N = 100;
    auto freqs = linspace(1e-3, 1e2, N);

    LISA lisa;
    std::vector<double> amp(N);
    for (size_t i = 0; i < N; ++i)
        amp[i] = 50.0 * lisa.omega_noise(std::max(freqs[i], 1e-4));

    auto results = get_SNR(freqs, amp);

    REQUIRE(results.size() == 2);
    REQUIRE(results[0].detector_name == "LISA");
    REQUIRE(results[1].detector_name == "DECIGO");

    SECTION("all SNRs are finite and non-negative") {
        for (const auto& r : results) {
            REQUIRE(std::isfinite(r.snr));
            REQUIRE(r.snr >= 0.0);
        }
    }
}

// ============================================================
// Individual convenience functions: LISA_snr, DECIGO_snr
// ============================================================
TEST_CASE("Individual detector convenience functions match calculate_snr()", "[snr]") {
    const size_t N = 100;
    auto freqs = linspace(1e-2, 10.0, N);

    LISA lisa;
    DECIGO decigo;

    std::vector<double> amp_lisa(N), amp_decigo(N);
    for (size_t i = 0; i < N; ++i) {
        amp_lisa[i]   = 20.0 * lisa.omega_noise(std::max(freqs[i], 1e-4));
        amp_decigo[i] = 20.0 * decigo.omega_noise(freqs[i]);
    }

    SECTION("LISA_snr matches calculate_snr") {
        REQUIRE(LISA_snr(freqs, amp_lisa) == Approx(calculate_snr(freqs, amp_lisa, lisa)).epsilon(1e-10));
    }
    SECTION("DECIGO_snr matches calculate_snr") {
        REQUIRE(DECIGO_snr(freqs, amp_decigo) == Approx(calculate_snr(freqs, amp_decigo, decigo)).epsilon(1e-10));
    }

    SECTION("all throw on size mismatch") {
        std::vector<double> bad(N - 1, 1e-10);
        REQUIRE_THROWS_AS(LISA_snr(freqs, bad), std::invalid_argument);
        REQUIRE_THROWS_AS(DECIGO_snr(freqs, bad), std::invalid_argument);
    }
}
