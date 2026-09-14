#include "catch/catch.hpp"
#include <string>
#include <utility>
#include <vector>
#include "logger.hpp"

namespace {

/// Captures records so a test can assert on what was emitted.
struct Capture {
    std::vector<std::pair<logging::Level, std::string>> records;

    void install() {
        logging::set_sink([this](logging::Level level, const std::string& msg) {
            records.push_back({level, msg});
        });
    }
};

/// Restores the level and the default sink however a test exits.
struct LoggingFixture {
    logging::Level saved = logging::get_level();
    ~LoggingFixture() {
        logging::set_level(saved);
        logging::set_sink(nullptr);
    }
};

} // namespace

// ============================================================
// sink
// ============================================================
TEST_CASE("custom sink receives records", "[logging]") {
    LoggingFixture fixture;
    Capture capture;
    capture.install();
    logging::set_level(logging::Level::Trace);

    SECTION("level and text are forwarded") {
        HG_LOG(Warn) << "shock residual " << 42;

        REQUIRE(capture.records.size() == 1);
        CHECK(capture.records[0].first == logging::Level::Warn);
        CHECK(capture.records[0].second == "shock residual 42");
    }

    SECTION("each record is delivered separately") {
        HG_LOG(Info) << "first";
        HG_LOG(Info) << "second";

        REQUIRE(capture.records.size() == 2);
        CHECK(capture.records[0].second == "first");
        CHECK(capture.records[1].second == "second");
    }

    SECTION("passing nullptr restores the default sink") {
        logging::set_sink(nullptr);
        HG_LOG(Info) << "goes to cout, not the capture";

        CHECK(capture.records.empty());
    }
}

// ============================================================
// severity threshold
// ============================================================
TEST_CASE("level filters records", "[logging]") {
    LoggingFixture fixture;
    Capture capture;
    capture.install();

    SECTION("records below the threshold are dropped") {
        logging::set_level(logging::Level::Warn);

        HG_LOG(Trace) << "trace";
        HG_LOG(Debug) << "debug";
        HG_LOG(Info) << "info";
        HG_LOG(Warn) << "warn";
        HG_LOG(Error) << "error";

        REQUIRE(capture.records.size() == 2);
        CHECK(capture.records[0].first == logging::Level::Warn);
        CHECK(capture.records[1].first == logging::Level::Error);
    }

    SECTION("Off silences everything") {
        logging::set_level(logging::Level::Off);

        HG_LOG(Trace) << "trace";
        HG_LOG(Info) << "info";
        HG_LOG(Warn) << "warn";
        HG_LOG(Error) << "error";

        CHECK(capture.records.empty());
    }

    SECTION("Trace emits everything") {
        logging::set_level(logging::Level::Trace);

        HG_LOG(Trace) << "trace";
        HG_LOG(Error) << "error";

        CHECK(capture.records.size() == 2);
    }

    SECTION("enabled() agrees with what is emitted") {
        logging::set_level(logging::Level::Info);

        CHECK_FALSE(logging::enabled(logging::Level::Trace));
        CHECK_FALSE(logging::enabled(logging::Level::Debug));
        CHECK(logging::enabled(logging::Level::Info));
        CHECK(logging::enabled(logging::Level::Warn));
        CHECK(logging::enabled(logging::Level::Error));
    }
}

// ============================================================
// argument evaluation
// ============================================================
TEST_CASE("filtered records do not evaluate their arguments", "[logging]") {
    LoggingFixture fixture;
    Capture capture;
    capture.install();
    logging::set_level(logging::Level::Warn);

    int calls = 0;
    const auto expensive = [&calls]() { return ++calls; };

    SECTION("below the threshold the expression never runs") {
        HG_LOG(Debug) << "value " << expensive();

        CHECK(calls == 0);
        CHECK(capture.records.empty());
    }

    SECTION("at or above the threshold it runs exactly once") {
        HG_LOG(Error) << "value " << expensive();

        CHECK(calls == 1);
        REQUIRE(capture.records.size() == 1);
        CHECK(capture.records[0].second == "value 1");
    }
}

// ============================================================
// ScopedLevel
// ============================================================
TEST_CASE("ScopedLevel restores the previous level", "[logging]") {
    LoggingFixture fixture;
    logging::set_level(logging::Level::Warn);

    SECTION("restores on normal scope exit") {
        {
            logging::ScopedLevel raise(logging::Level::Debug);
            CHECK(logging::get_level() == logging::Level::Debug);
        }
        CHECK(logging::get_level() == logging::Level::Warn);
    }

    SECTION("restores when unwound by an exception") {
        try {
            logging::ScopedLevel raise(logging::Level::Trace);
            throw std::runtime_error("boom");
        } catch (const std::runtime_error&) {
        }
        CHECK(logging::get_level() == logging::Level::Warn);
    }

    SECTION("nests") {
        {
            logging::ScopedLevel outer(logging::Level::Info);
            {
                logging::ScopedLevel inner(logging::Level::Trace);
                CHECK(logging::get_level() == logging::Level::Trace);
            }
            CHECK(logging::get_level() == logging::Level::Info);
        }
        CHECK(logging::get_level() == logging::Level::Warn);
    }
}

// ============================================================
// macro hygiene
// ============================================================
TEST_CASE("HG_LOG is safe as an unbraced if body", "[logging]") {
    LoggingFixture fixture;
    Capture capture;
    capture.install();
    logging::set_level(logging::Level::Trace);

    // The else must bind to this if, not to anything inside the macro.
    const bool condition = false;
    if (condition)
        HG_LOG(Info) << "taken";
    else
        HG_LOG(Warn) << "not taken";

    REQUIRE(capture.records.size() == 1);
    CHECK(capture.records[0].first == logging::Level::Warn);
    CHECK(capture.records[0].second == "not taken");
}
