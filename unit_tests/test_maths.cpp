#include "catch/catch.hpp"
#include <random>
#include <chrono>
#include <cmath>
#include <interpolation.h> // ALGLIB
#include "maths_ops.hpp"

// ============================================================
// linspace / logspace
// ============================================================
TEST_CASE("linspace", "[linspace]") {
    SECTION("correct number of points") {
        auto v = linspace(0.0, 1.0, 11);
        REQUIRE(v.size() == 11);
    }
    SECTION("endpoints correct") {
        auto v = linspace(2.0, 5.0, 50);
        REQUIRE(v.front() == Approx(2.0).epsilon(1e-12));
        REQUIRE(v.back()  == Approx(5.0).epsilon(1e-12));
    }
    SECTION("uniform spacing") {
        auto v = linspace(0.0, 1.0, 5);
        REQUIRE(v[1] == Approx(0.25).epsilon(1e-12));
        REQUIRE(v[2] == Approx(0.50).epsilon(1e-12));
        REQUIRE(v[3] == Approx(0.75).epsilon(1e-12));
    }
    SECTION("single point") {
        auto v = linspace(3.0, 3.0, 1);
        REQUIRE(v.size() == 1);
        REQUIRE(v[0] == Approx(3.0).epsilon(1e-12));
    }
}

TEST_CASE("logspace", "[logspace]") {
    SECTION("correct number of points") {
        auto v = logspace(-3.0, 3.0, 7);
        REQUIRE(v.size() == 7);
    }
    SECTION("endpoints correct") {
        auto v = logspace(-2.0, 2.0, 100);
        REQUIRE(v.front() == Approx(1e-2).epsilon(1e-10));
        REQUIRE(v.back()  == Approx(1e+2).epsilon(1e-10));
    }
    SECTION("values are strictly positive") {
        auto v = logspace(-3.0, 3.0, 50);
        for (double x : v) REQUIRE(x > 0.0);
    }
    SECTION("values are strictly increasing") {
        auto v = logspace(-3.0, 3.0, 50);
        for (size_t i = 1; i < v.size(); ++i)
            REQUIRE(v[i] > v[i-1]);
    }
}

// ============================================================
// power()
// ============================================================
TEST_CASE("power()", "[power]") {
    REQUIRE(power(2.0, 0)  == Approx(1.0).epsilon(1e-12));
    REQUIRE(power(2.0, 3)  == Approx(8.0).epsilon(1e-12));
    REQUIRE(power(3.0, 4)  == Approx(81.0).epsilon(1e-12));
    REQUIRE(power(2.0, -2) == Approx(0.25).epsilon(1e-12));
    REQUIRE(power(0.0, 5)  == Approx(0.0).margin(1e-15));
}

// ============================================================
// simpson_integrate()
// ============================================================
TEST_CASE("simpson_integrate()", "[simpson]") {
    SECTION("integral of constant 1 over [0,1] = 1") {
        auto x = linspace(0.0, 1.0, 101);
        std::vector<double> y(x.size(), 1.0);
        REQUIRE(simpson_integrate(x, y) == Approx(1.0).epsilon(1e-6));
    }
    SECTION("integral of x over [0,1] = 0.5") {
        auto x = linspace(0.0, 1.0, 101);
        std::vector<double> y(x.size());
        for (size_t i = 0; i < x.size(); ++i) y[i] = x[i];
        REQUIRE(simpson_integrate(x, y) == Approx(0.5).epsilon(1e-6));
    }
    SECTION("integral of x^2 over [0,1] = 1/3") {
        auto x = linspace(0.0, 1.0, 201);
        std::vector<double> y(x.size());
        for (size_t i = 0; i < x.size(); ++i) y[i] = x[i] * x[i];
        REQUIRE(simpson_integrate(x, y) == Approx(1.0/3.0).epsilon(1e-5));
    }
    SECTION("integral of sin(x) over [0, pi] = 2") {
        auto x = linspace(0.0, M_PI, 501);
        std::vector<double> y(x.size());
        for (size_t i = 0; i < x.size(); ++i) y[i] = std::sin(x[i]);
        REQUIRE(simpson_integrate(x, y) == Approx(2.0).epsilon(1e-6));
    }
    SECTION("throws on size mismatch") {
        std::vector<double> x = {0.0, 1.0, 2.0};
        std::vector<double> y = {0.0, 1.0};
        REQUIRE_THROWS_AS(simpson_integrate(x, y), std::invalid_argument);
    }
    SECTION("throws on fewer than two points") {
        std::vector<double> x = {1.0};
        std::vector<double> y = {1.0};
        REQUIRE_THROWS_AS(simpson_integrate(x, y), std::invalid_argument);
    }
}

// ============================================================
// CubicSpline
// ============================================================
TEST_CASE("CubicSpline — cubic polynomial reproduced exactly", "[cubicSpline]") {
    auto f = [](double x) { return x*x*x - 2*x*x + x - 5; };

    auto x_vals = linspace(-2.0, 2.0, 50);
    std::vector<double> y_vals;
    for (double x : x_vals) y_vals.push_back(f(x));

    CubicSpline<double> spline(x_vals, y_vals);

    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(-2.0, 2.0);
    for (int i = 0; i < 50; ++i) {
        double xi = dis(gen);
        CHECK(spline(xi) == Approx(f(xi)).epsilon(1e-4));
    }
}

TEST_CASE("CubicSpline — reproduces knot values exactly", "[cubicSpline]") {
    auto x_vals = linspace(0.0, 5.0, 20);
    std::vector<double> y_vals;
    for (double x : x_vals) y_vals.push_back(std::sin(x));

    CubicSpline<double> spline(x_vals, y_vals);

    for (size_t i = 0; i < x_vals.size(); ++i)
        CHECK(spline(x_vals[i]) == Approx(y_vals[i]).epsilon(1e-10));
}

TEST_CASE("CubicSpline — decreasing x input", "[cubicSpline]") {
    // Build on decreasing grid — spline should handle reversal internally
    std::vector<double> x_vals, y_vals;
    for (int i = 50; i >= 0; --i) {
        double x = -2.0 + i * (4.0 / 50.0);
        x_vals.push_back(x);
        y_vals.push_back(x * x);
    }
    CubicSpline<double> spline(x_vals, y_vals);
    CHECK(spline(0.0)  == Approx(0.0).margin(1e-4));
    CHECK(spline(1.0)  == Approx(1.0).epsilon(1e-4));
    CHECK(spline(-1.0) == Approx(1.0).epsilon(1e-4));
}

TEST_CASE("CubicSpline — scalar arithmetic operators", "[cubicSpline]") {
    auto x_vals = linspace(0.0, 2.0, 20);
    std::vector<double> y_vals;
    for (double x : x_vals) y_vals.push_back(x * x);

    CubicSpline<double> s(x_vals, y_vals);

    SECTION("spline + scalar") {
        auto s2 = s + 3.0;
        CHECK(s2(1.0) == Approx(s(1.0) + 3.0).epsilon(1e-10));
    }
    SECTION("scalar + spline") {
        auto s2 = 3.0 + s;
        CHECK(s2(1.0) == Approx(s(1.0) + 3.0).epsilon(1e-10));
    }
    SECTION("spline - scalar") {
        auto s2 = s - 1.5;
        CHECK(s2(1.0) == Approx(s(1.0) - 1.5).epsilon(1e-10));
    }
    SECTION("scalar - spline") {
        auto s2 = 5.0 - s;
        CHECK(s2(1.0) == Approx(5.0 - s(1.0)).epsilon(1e-10));
    }
    SECTION("spline * scalar") {
        auto s2 = s * 2.0;
        CHECK(s2(1.0) == Approx(s(1.0) * 2.0).epsilon(1e-10));
    }
    SECTION("scalar * spline") {
        auto s2 = 2.0 * s;
        CHECK(s2(1.0) == Approx(s(1.0) * 2.0).epsilon(1e-10));
    }
    SECTION("spline / scalar") {
        auto s2 = s / 4.0;
        CHECK(s2(1.0) == Approx(s(1.0) / 4.0).epsilon(1e-10));
    }
}

// ============================================================
// golden_section_minimize()
// ============================================================
TEST_CASE("golden_section_minimize()", "[golden_section]") {
    SECTION("minimum of x^2 at x=0") {
        auto f = [](double x) { return x * x; };
        double xmin = golden_section_minimize(f, -2.0, 2.0, 1e-8);
        REQUIRE(xmin == Approx(0.0).margin(1e-6));
    }
    SECTION("minimum of (x-1.5)^2 at x=1.5") {
        auto f = [](double x) { return (x - 1.5) * (x - 1.5); };
        double xmin = golden_section_minimize(f, 0.0, 3.0, 1e-8);
        REQUIRE(xmin == Approx(1.5).epsilon(1e-6));
    }
    SECTION("minimum of cos(x) on [pi/2, 3pi/2] at x=pi") {
        auto f = [](double x) { return std::cos(x); };
        double xmin = golden_section_minimize(f, M_PI / 2.0, 3.0 * M_PI / 2.0, 1e-8);
        REQUIRE(xmin == Approx(M_PI).epsilon(1e-6));
    }
}

// ============================================================
// rk4_solver()
// ============================================================
TEST_CASE("rk4_solver — simple harmonic oscillator", "[rk4]") {
    // dy0/dx = y1, dy1/dx = -y0 => y0 = cos(x), y1 = -sin(x)
    auto dydx = [](double, const std::vector<double>& y) -> std::vector<double> {
        return { y[1], -y[0] };
    };

    auto [x_vals, y_vals] = rk4_solver(dydx, 0.0, M_PI / 2.0,
                                        std::vector<double>{1.0, 0.0}, size_t(5000));
    const double tol = 1e-3;
    const double eps = 1e-10;

    for (size_t i = 0; i < x_vals.size(); ++i) {
        double x = x_vals[i];
        double e0 = std::abs(y_vals[i][0] - std::cos(x));
        double e1 = std::abs(y_vals[i][1] + std::sin(x));
        double d0 = std::max(std::abs(std::cos(x)), eps);
        double d1 = std::max(std::abs(std::sin(x)), eps);
        REQUIRE(e0 / d0 <= tol);
        REQUIRE(e1 / d1 <= tol);
    }
}

TEST_CASE("rk4_solver — exponential growth", "[rk4]") {
    // dy/dx = y, y(0) = 1 => y = exp(x)
    auto dydx = [](double, const std::vector<double>& y) -> std::vector<double> {
        return { y[0] };
    };
    auto [x_vals, y_vals] = rk4_solver(dydx, 0.0, 2.0,
                                        std::vector<double>{1.0}, size_t(1000));
    for (size_t i = 0; i < x_vals.size(); ++i)
        REQUIRE(y_vals[i][0] == Approx(std::exp(x_vals[i])).epsilon(1e-4));
}

TEST_CASE("rk4_solver — correct number of steps", "[rk4]") {
    auto dydx = [](double, const std::vector<double>& y) -> std::vector<double> {
        return { y[0] };
    };
    const size_t n = 200;
    auto [x_vals, y_vals] = rk4_solver(dydx, 0.0, 1.0,
                                        std::vector<double>{1.0}, n);
    REQUIRE(x_vals.size() == n + 1);
    REQUIRE(y_vals.size() == n + 1);
}

// ============================================================
// L1_norm() and L2_norm()
// ============================================================
TEST_CASE("L1_norm()", "[norms]") {
    SECTION("identical distributions give zero") {
        auto x = linspace(0.0, 1.0, 101);
        std::vector<double> y(101, 1.0);
        REQUIRE(L1_norm(x, y, y) == Approx(0.0).margin(1e-10));
    }
    SECTION("known analytic value: |f1 - f2| = |x - (1-x)| = |2x-1| over [0,1]") {
        auto x = linspace(0.0, 1.0, 1001);
        std::vector<double> f1(x.size()), f2(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            f1[i] = x[i];
            f2[i] = 1.0 - x[i];
        }
        // integral of |2x-1| over [0,1] = 0.5
        REQUIRE(L1_norm(x, f1, f2) == Approx(0.5).epsilon(1e-3));
    }
    SECTION("throws on size mismatch") {
        std::vector<double> x = {0.0, 1.0};
        std::vector<double> y1 = {0.0, 1.0};
        std::vector<double> y2 = {0.0};
        REQUIRE_THROWS_AS(L1_norm(x, y1, y2), std::invalid_argument);
    }
}

TEST_CASE("L2_norm()", "[norms]") {
    SECTION("identical distributions give zero") {
        auto x = linspace(0.0, 1.0, 101);
        std::vector<double> y(101, 1.0);
        REQUIRE(L2_norm(x, y, y) == Approx(0.0).margin(1e-10));
    }
    SECTION("non-negative") {
        auto x = linspace(0.0, 1.0, 51);
        std::vector<double> f1(51, 2.0), f2(51, 1.0);
        REQUIRE(L2_norm(x, f1, f2) >= 0.0);
    }
    SECTION("throws on size mismatch") {
        std::vector<double> x = {0.0, 1.0};
        std::vector<double> y1 = {0.0, 1.0};
        std::vector<double> y2 = {0.0};
        REQUIRE_THROWS_AS(L2_norm(x, y1, y2), std::invalid_argument);
    }
}

// ============================================================
// Wasserstein distance
// ============================================================
TEST_CASE("wasserstein_distance_1d()", "[wasserstein]") {
    const size_t N = 401;
    auto log_freq = linspace(-6.0, 6.0, N);

    auto gaussian = [](double x, double mu, double sig) {
        return std::exp(-0.5 * std::pow((x - mu) / sig, 2));
    };

    const double mu1 = 0.0, mu2 = 0.5, sigma = 0.4;
    std::vector<double> s1(N), s2(N), s3(N);
    for (size_t i = 0; i < N; ++i) {
        s1[i] = gaussian(log_freq[i], mu1, sigma);
        s2[i] = gaussian(log_freq[i], mu2, sigma);
        s3[i] = s1[i];
    }

    SECTION("zero distance for identical spectra") {
        double W = wasserstein_distance_1d(log_freq, s1, log_freq, s3);
        REQUIRE(W == Approx(0.0).margin(1e-6));
    }
    SECTION("correct shift distance") {
        double W = wasserstein_distance_1d(log_freq, s1, log_freq, s2);
        REQUIRE(W == Approx(std::abs(mu2 - mu1)).epsilon(1e-3));
    }
    SECTION("symmetry") {
        double W12 = wasserstein_distance_1d(log_freq, s1, log_freq, s2);
        double W21 = wasserstein_distance_1d(log_freq, s2, log_freq, s1);
        REQUIRE(W12 == Approx(W21).epsilon(1e-10));
    }
    SECTION("non-negative") {
        double W = wasserstein_distance_1d(log_freq, s1, log_freq, s2);
        REQUIRE(W >= 0.0);
    }
    SECTION("scale invariance of weights") {
        std::vector<double> s1_scaled = s1;
        for (auto& v : s1_scaled) v *= 3.0;
        double W1 = wasserstein_distance_1d(log_freq, s1, log_freq, s2);
        double W2 = wasserstein_distance_1d(log_freq, s1_scaled, log_freq, s2);
        REQUIRE(W1 == Approx(W2).epsilon(1e-10));
    }
    SECTION("throws on size mismatch") {
        std::vector<double> short_v(10, 1.0);
        REQUIRE_THROWS(wasserstein_distance_1d(log_freq, short_v, log_freq, s2));
        REQUIRE_THROWS(wasserstein_distance_1d(log_freq, s1, log_freq, short_v));
    }
}

// ============================================================
// ALGLIB spline (regression test — checks library is linked)
// ============================================================
TEST_CASE("ALGLIB spline regression", "[alglib]") {
    auto f = [](double x) { return x*x*x - 2*x*x + x - 5; };
    const int n = 50;
    alglib::real_1d_array x_arr, y_arr;
    std::vector<double> xv(n), yv(n);
    for (int i = 0; i < n; ++i) {
        xv[i] = -2.0 + i * (4.0 / (n - 1));
        yv[i] = f(xv[i]);
    }
    x_arr.setcontent(n, xv.data());
    y_arr.setcontent(n, yv.data());

    alglib::spline1dinterpolant spline;
    alglib::spline1dbuildcubic(x_arr, y_arr, spline);

    std::mt19937 gen(0);
    std::uniform_real_distribution<> dis(-2.0, 2.0);
    for (int i = 0; i < 30; ++i) {
        double xi = dis(gen);
        CHECK(alglib::spline1dcalc(spline, xi) == Approx(f(xi)).epsilon(1e-4));
    }
}
