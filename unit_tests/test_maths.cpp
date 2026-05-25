#include "catch/catch.hpp"
#include <random>
#include <chrono>
#include <interpolation.h> //ALGLIB
#include "maths_ops.hpp"

TEST_CASE("Test rk4 Coupled ODEs", "[rk4]") {

    // dy0/dx = y1, dy1/dx = -y0 with y(0) = [1, 0] => y0(x) = cos(x), y1(x) = -sin(x)

    auto dydx = [](double x, const std::vector<double>& y) -> std::vector<double> {
        return { y[1], -y[0] };
    };

    const double x0 = 0.0;
    const double xf = M_PI / 2; // pi/2
    const size_t steps = 5000;
    const double tol = 1e-3;
    const double epsilon = 1e-10; // for safe relative error denom

    const std::vector<double> y0 = {1.0, 0.0};

    auto [x_vals, y_vals] = rk4_solver(dydx, x0, xf, y0, steps);

    for (size_t i = 0; i < x_vals.size(); ++i) {
        double x = x_vals[i];
        double y0_num = y_vals[i][0];
        double y1_num = y_vals[i][1];

        double y0_exact = std::cos(x);
        double y1_exact = -std::sin(x);

        double denom0 = std::abs(y0_exact) > epsilon ? std::abs(y0_exact) : 1.0;
        double denom1 = std::abs(y1_exact) > epsilon ? std::abs(y1_exact) : 1.0;

        double rel_err0 = std::abs(y0_num - y0_exact) / denom0;
        double rel_err1 = std::abs(y1_num - y1_exact) / denom1;

        REQUIRE(rel_err0 <= tol);
        REQUIRE(rel_err1 <= tol);
    }

}

TEST_CASE("Test rk4 Solver", "[rk4]") {

    // dy0/dx = y1, dy1/dx = -y0 with y(0) = [1, 0] => y0(x) = cos(x), y1(x) = -sin(x)

    auto dydx = [](double x, const std::vector<double>& y) -> std::vector<double> {
        return { y[1], -y[0] };
    };

    const double x0 = 0.0;
    const double xf = M_PI / 2; // pi/2
    const size_t steps = 5000;
    const double tol = 1e-3;
    const double epsilon = 1e-10; // for safe relative error denom

    const std::vector<double> y0 = {1.0, 0.0};

    auto [x_vals, y_vals] = rk4_solver(dydx, x0, xf, y0, steps);

    for (size_t i = 0; i < x_vals.size(); ++i) {
        double x = x_vals[i];
        double y0_num = y_vals[i][0];
        double y1_num = y_vals[i][1];

        double y0_exact = std::cos(x);
        double y1_exact = -std::sin(x);

        double denom0 = std::abs(y0_exact) > epsilon ? std::abs(y0_exact) : 1.0;
        double denom1 = std::abs(y1_exact) > epsilon ? std::abs(y1_exact) : 1.0;

        double rel_err0 = std::abs(y0_num - y0_exact) / denom0;
        double rel_err1 = std::abs(y1_num - y1_exact) / denom1;

        REQUIRE(rel_err0 <= tol);
        REQUIRE(rel_err1 <= tol);
    }
    
}

TEST_CASE("Test cubic spline interpolation", "[cubicSpline]") {

    auto f = [](double x) { return x*x*x - 2*x*x + x - 5; };

    std::vector<double> x_vals = linspace(-2.0, 2.0, 50);
    std::vector<double> y_vals;
    for (double x : x_vals) {
        y_vals.push_back(f(x));
    }

    auto start = std::chrono::high_resolution_clock::now();
    CubicSpline<double> spline(x_vals, y_vals);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "CubicSpline build time: " << duration.count() << " s\n";

    auto start_eval = std::chrono::high_resolution_clock::now();
    spline(1.0);
    auto end_eval = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_eval = end_eval - start_eval;
    std::cout << "CubicSpline run_time: " << duration_eval.count() << " s\n";

    std::vector<double> test_points;
    std::mt19937 gen(0);
    std::uniform_real_distribution<> dis(-2.0, 2.0);
    for (int i = 0; i < 50; ++i) {test_points.push_back(dis(gen));}

    for (double xi : test_points) {

        double yi_exact = f(xi);
        double yi_interp = spline(xi);
        
        CHECK(yi_interp == Approx(yi_exact).epsilon(1e-4));
    }
}

TEST_CASE("ALGLIB function", "[alglib]") {

    auto f = [](double x) { return x*x*x - 2*x*x + x - 5; };

    alglib::real_1d_array x_vals, y_vals;
    const int n = 50;
    std::vector<double> x_vec(n), y_vec(n);

    for (int i = 0; i < n; i++) {
        double x = -2.0 + i * 0.08; // -2 to 2 in 50 steps
        double y = f(x);
        x_vec[i] = x;
        y_vec[i] = y;
    }

    x_vals.setcontent(n, x_vec.data());
    y_vals.setcontent(n, y_vec.data());

    alglib::spline1dinterpolant spline;

    auto start = std::chrono::high_resolution_clock::now();
    alglib::spline1dbuildcubic(x_vals, y_vals, spline);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "AlgLib build time: " << duration.count() << " s\n";

    auto start_eval = std::chrono::high_resolution_clock::now();
    alglib::spline1dcalc(spline, 1.0);
    auto end_eval = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_eval = end_eval - start_eval;
    std::cout << "AlgLib run_time: " << duration_eval.count() << " s\n";

    std::vector<double> test_points;
    std::mt19937 gen(0);
    std::uniform_real_distribution<> dis(-2.0, 2.0);
    for (int i = 0; i < 50; ++i) {test_points.push_back(dis(gen));}

    double tol = 1e-4;
    for (double xi : test_points) {

        double yi_exact = f(xi);
        double yi_interp = alglib::spline1dcalc(spline, xi);
        double error = std::abs(yi_interp - yi_exact);
        
        CHECK(error < tol);
    }

}

TEST_CASE("Test 1D Wasserstein distance on spectra", "[wasserstein]") {

    const size_t N = 401;
    const double fmin = 1e-3;
    const double fmax = 1e3;
    const double tol  = 1e-3;

    // Log-spaced frequencies
    std::vector<double> freq(N);
    for (size_t i = 0; i < N; ++i) {
        double t = static_cast<double>(i) / (N - 1);
        freq[i] = fmin * std::pow(fmax / fmin, t);
    }

    // Transform to log-frequency (positions for Wasserstein distance)
    std::vector<double> log_freq(N);
    for (size_t i = 0; i < N; ++i) {
        log_freq[i] = std::log(freq[i]);
    }

    // Gaussian in log-frequency
    auto log_gaussian = [](double x, double mu, double sigma) {
        return std::exp(-0.5 * std::pow((x - mu) / sigma, 2));
    };

    std::vector<double> s1(N), s2(N), s3(N);

    const double mu1 = 0.0;
    const double mu2 = 0.5;     // known shift in log-space
    const double sigma = 0.4;

    for (size_t i = 0; i < N; ++i) {
        double lx = std::log(freq[i]);
        s1[i] = log_gaussian(lx, mu1, sigma);
        s2[i] = log_gaussian(lx, mu2, sigma);
        s3[i] = log_gaussian(lx, mu1, sigma); // identical to s1
    }

    SECTION("Zero distance for identical spectra") {
        // log_freq is positions, s1 and s3 are weights
        double W = wasserstein_distance_1d(log_freq, s1, log_freq, s3);
        REQUIRE(W == Approx(0.0).margin(1e-6));
    }

    SECTION("Correct distance for shifted distributions") {
        // For identical shapes shifted by Δ in log-frequency space:
        // W1 ≈ |Δ| when distributions are well-separated relative to their width
        double W = wasserstein_distance_1d(log_freq, s1, log_freq, s2);
        
        // The Wasserstein distance for Gaussians with identical variance
        // but shifted means is exactly |μ2 - μ1|
        REQUIRE(W == Approx(std::abs(mu2 - mu1)).epsilon(tol));
    }

    SECTION("Symmetry of Wasserstein distance") {
        double W12 = wasserstein_distance_1d(log_freq, s1, log_freq, s2);
        double W21 = wasserstein_distance_1d(log_freq, s2, log_freq, s1);

        REQUIRE(W12 == Approx(W21).epsilon(1e-10));
    }

    SECTION("Triangle inequality") {
        // For three distributions: W(s1,s3) <= W(s1,s2) + W(s2,s3)
        std::vector<double> s4(N);
        const double mu4 = 1.0;
        for (size_t i = 0; i < N; ++i) {
            double lx = std::log(freq[i]);
            s4[i] = log_gaussian(lx, mu4, sigma);
        }
        
        double W13 = wasserstein_distance_1d(log_freq, s1, log_freq, s3);
        double W12 = wasserstein_distance_1d(log_freq, s1, log_freq, s2);
        double W23 = wasserstein_distance_1d(log_freq, s2, log_freq, s3);
        
        REQUIRE(W13 <= W12 + W23 + 1e-10);
    }

    SECTION("Handles non-normalized weights correctly") {
        // Wasserstein distance normalizes weights internally
        std::vector<double> s1_scaled = s1;
        for (auto& val : s1_scaled) val *= 2.0;
        
        double W_normal = wasserstein_distance_1d(log_freq, s1, log_freq, s2);
        double W_scaled = wasserstein_distance_1d(log_freq, s1_scaled, log_freq, s2);
        
        REQUIRE(W_normal == Approx(W_scaled).epsilon(1e-10));
    }

    SECTION("Throws on size mismatch") {
        std::vector<double> short_freq(10);
        std::vector<double> short_spec(10);
        
        // Mismatch between u_values and u_weights
        REQUIRE_THROWS(wasserstein_distance_1d(log_freq, short_spec, log_freq, s2));
        
        // Mismatch between v_values and v_weights
        REQUIRE_THROWS(wasserstein_distance_1d(log_freq, s1, short_freq, s2));
    }

    SECTION("Throws on insufficient points") {
        std::vector<double> f = {1.0, 2.0};
        std::vector<double> s = {1.0, 1.0};
        REQUIRE_THROWS(wasserstein_distance_1d(f, s, f, s));
    }

    SECTION("Throws on negative or zero weights") {
        auto bad_spec = s1;
        bad_spec[10] = -1.0;
        REQUIRE_THROWS(wasserstein_distance_1d(log_freq, bad_spec, log_freq, s2));
        
        bad_spec[10] = 0.0;
        // Zero weights should be allowed (just ignored)
        REQUIRE_NOTHROW(wasserstein_distance_1d(log_freq, bad_spec, log_freq, s2));
    }

    SECTION("Works with different frequency grids") {
        // Create a different frequency grid for s2
        std::vector<double> freq2(N);
        std::vector<double> log_freq2(N);
        std::vector<double> s2_different_grid(N);
        
        for (size_t i = 0; i < N; ++i) {
            double t = static_cast<double>(i) / (N - 1);
            freq2[i] = fmin * std::pow(fmax / fmin, t) * 1.1; // slightly shifted grid
            log_freq2[i] = std::log(freq2[i]);
            s2_different_grid[i] = log_gaussian(log_freq2[i], mu2, sigma);
        }
        
        // Should still compute distance correctly
        double W = wasserstein_distance_1d(log_freq, s1, log_freq2, s2_different_grid);
        REQUIRE(W > 0.0);
        REQUIRE(W == Approx(std::abs(mu2 - mu1)).epsilon(0.1)); // larger tolerance due to grid mismatch
    }
}