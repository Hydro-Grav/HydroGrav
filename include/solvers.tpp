#ifndef INCLUDE_SOLVERS_TPP_H
#define INCLUDE_SOLVERS_TPP_H

#include <stdexcept>
#include <concepts>
#include <vector>
#include <array>
#include <utility>
#include <functional>
#include <cassert>
#include <cmath>

template <typename T, typename Func>
T newton_solve_1d(
    const Func& F,  // F(x): returns scalar
    T x0,           // initial guess
    T tol,          // tolerance
    int max_iter,   // maximum iterations
    T h             // finite difference step for derivative
) {
    static_assert(std::is_floating_point<T>::value,
                  "newton_solve_1d requires floating-point T");

    for (int iter = 0; iter < max_iter; ++iter) {
        T fx = F(x0);
        if (std::fabs(fx) < tol)
            return x0;

        // Numerical derivative via finite difference
        T dfdx = (F(x0 + h) - fx) / h;
        if (std::fabs(dfdx) < T(1e-14))
            throw std::runtime_error("Derivative is near zero (Jacobian singular).");

        // Newton update
        T dx = -fx / dfdx;
        x0 += dx;

        if (std::fabs(dx) < tol)
            return x0;
    }

    throw std::runtime_error("Newton's method 1D did not converge");
}

template <typename T, typename Func>
std::array<T,2> newton_solve_2d(
    const Func& F,
    std::array<T,2> x0,
    T tol,
    int max_iter,
    T h
) {
    // Check T is floating point
    static_assert(std::is_floating_point<T>::value,
                  "rk4_solver requires floating-point T");

    // Check that calling dydx(x0, y0) returns a State
    using result_t = decltype(F(x0));
    static_assert(std::is_same<result_t, std::array<T, 2>>::value,
                  "rk4_solver requires dydx(T, std::array<T, 2>&) -> std::array<T, 2>");

    // Inline solver for 2×2 linear systems
    auto solve_linear_2x2 = [](const std::array<std::array<T,2>,2>& A,
                               const std::array<T,2>& b) {
        T det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
        if (std::fabs(det) < T(1e-14))
            throw std::runtime_error("Jacobian is singular!");
        std::array<T,2> x;
        x[0] = ( b[0]*A[1][1] - b[1]*A[0][1]) / det;
        x[1] = ( A[0][0]*b[1] - A[1][0]*b[0]) / det;
        return x;
    };

    for (int iter = 0; iter < max_iter; ++iter) {
        // Evaluate F(x0) once
        auto fx = F(x0);

        T norm = std::sqrt(fx[0]*fx[0] + fx[1]*fx[1]);
        if (norm < tol) return x0;

        // Numerical Jacobian with reuse of F(x0)
        std::array<std::array<T,2>,2> J{};
        for (int j = 0; j < 2; ++j) {
            auto xh = x0;
            xh[j] += h;
            auto fh = F(xh);
            for (int i = 0; i < 2; ++i) {
                J[i][j] = (fh[i] - fx[i]) / h;
            }
        }

        std::array<T,2> minus_fx = { -fx[0], -fx[1] };
        auto dx = solve_linear_2x2(J, minus_fx);

        x0[0] += dx[0];
        x0[1] += dx[1];

        if (std::sqrt(dx[0]*dx[0] + dx[1]*dx[1]) < tol)
            return x0;
    }

    throw std::runtime_error("Newton's method 2D did not converge");
}

// RK4 solver: works with std::vector<T> or std::array<T, N>
template <typename T, typename State, typename Func>
std::pair<std::vector<T>, std::vector<State>>
rk4_solver(
    const Func& dydx,
    T x0,
    T xf,
    const State& y0,
    size_t n
) {
    // Check T is floating point
    static_assert(std::is_floating_point<T>::value,
                  "rk4_solver requires floating-point T");

    // Check that calling dydx(x0, y0) returns a State
    using result_t = decltype(dydx(x0, y0));
    static_assert(std::is_same<result_t, State>::value,
                  "rk4_solver requires dydx(T, const State&) -> State");

    assert(n >= 2 && "Number of steps must be at least 2.");

    const T h = (xf - x0) / static_cast<T>(n - 1);
    std::vector<T> x_vals(n);
    std::vector<State> y_vals(n);

    x_vals[0] = x0;
    y_vals[0] = y0;

    T x = x0;
    State y = y0;

    for (size_t i = 1; i < n; ++i) {
        State k1 = dydx(x, y);

        State y_tmp = y;
        for (size_t j = 0; j < y.size(); ++j) y_tmp[j] = y[j] + 0.5*h*k1[j];
        State k2 = dydx(x + 0.5*h, y_tmp);

        for (size_t j = 0; j < y.size(); ++j) y_tmp[j] = y[j] + 0.5*h*k2[j];
        State k3 = dydx(x + 0.5*h, y_tmp);

        for (size_t j = 0; j < y.size(); ++j) y_tmp[j] = y[j] + h*k3[j];
        State k4 = dydx(x + h, y_tmp);

        for (size_t j = 0; j < y.size(); ++j) {
            y[j] += (h/6.0) * (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);
        }

        x += h;
        x_vals[i] = x;
        y_vals[i] = y;
    }

    return {x_vals, y_vals};
}

#endif