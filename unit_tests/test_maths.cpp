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

TEST_CASE("Test si(x) and ci(x)", "[trigInt]") {

    struct TestContainer {
        double x;
        double expected_si;
        double expected_ci;
    };

    TestContainer tests[] = {
        {-5.0, -1.5499312449, -0.19003},
        {-4.0, -1.7582031389, -0.140982},
        {-3.0, -1.8486525274, 0.11963},
        {-2.0, -1.6054129768, 0.422981},
        {-1.0, -0.9460830704, 0.337404 },
        {-0.5, -0.4931074180, -0.177784},
        { 0.5,  0.4931074180, -0.177784},
        { 1.0,  0.9460830704, 0.337404},
        { 2.0,  1.6054129768,  0.422981},
        { 3.0,  1.8486525274,  0.11963},
        { 4.0,  1.7582031389,  -0.140982},
        { 5.0,  1.5499312449,  -0.19003}
    };

    for (const auto& t : tests) {

        const auto [si, ci] = SiCi(t.x);

        CHECK(si == Approx(t.expected_si).epsilon(1e-5));
        CHECK(ci == Approx(t.expected_ci).epsilon(1e-5));
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