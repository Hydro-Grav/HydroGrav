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
T newton_solve_1d_bounded(
    const Func& f,
    T x0,
    T x_min,
    T x_max,
    T tol,
    int max_iter,
    T h,
    T bound_margin
) {
    static_assert(std::is_floating_point<T>::value,
                  "newton_solve_1d_bounded requires floating-point T");

    // Clamp initial guess to valid range
    x0 = std::clamp(x0, x_min + bound_margin, x_max - bound_margin);

    for (int iter = 0; iter < max_iter; ++iter) {
        // Evaluate function
        T fx = f(x0);
        
        // Check convergence
        if (std::abs(fx) < tol) return x0;

        // Compute derivative with adaptive one-sided differences
        T df;
        bool use_forward = (x0 + h <= x_max - bound_margin);
        bool use_backward = (x0 - h >= x_min + bound_margin);
        
        if (use_forward) {
            // Forward difference
            T fh = f(x0 + h);
            df = (fh - fx) / h;
        } else if (use_backward) {
            // Backward difference
            T fh = f(x0 - h);
            df = (fx - fh) / h;
        } else {
            // At boundary - use small central difference
            T h_small = std::min(x0 - x_min, x_max - x0) * 0.1;
            T fh_plus = f(x0 + h_small);
            T fh_minus = f(x0 - h_small);
            df = (fh_plus - fh_minus) / (2.0 * h_small);
        }

        // Check for zero derivative
        if (std::abs(df) < 1e-14) {
            // Derivative too small - try bisection step instead
            if (fx > 0) {
                x0 -= 0.01 * (x_max - x_min);
            } else {
                x0 += 0.01 * (x_max - x_min);
            }
            x0 = std::clamp(x0, x_min + bound_margin, x_max - bound_margin);
            continue;
        }

        // Compute Newton step
        T dx = -fx / df;

        // Step limiting
        T alpha = 1.0;
        T x_new = x0 + dx;
        
        if (x_new > x_max - bound_margin) {
            // Would exceed upper bound
            alpha = (x_max - bound_margin - x0) / dx;
        } else if (x_new < x_min + bound_margin) {
            // Would exceed lower bound
            alpha = (x_min + bound_margin - x0) / dx;
        }

        // Apply limited step
        x0 += alpha * dx;

        // Safety clamp
        x0 = std::clamp(x0, x_min + bound_margin, x_max - bound_margin);

        // Check step convergence
        if (std::abs(alpha * dx) < tol) {
            return x0;
        }
    }

    throw std::runtime_error("newton_solve_1d_bounded did not converge");
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

template <typename T, typename Func>
std::array<T,2> newton_solve_2d_bounded(
    const Func& F,
    std::array<T,2> x0,
    std::array<T,2> x_min,  // Lower bounds
    std::array<T,2> x_max,  // Upper bounds
    T tol,
    int max_iter,
    T h,
    T bound_margin,
    bool dev
) {
    if (dev) {
        std::cout << "Running newton_solve_2d_bounded() in dev mode!\n";
    }

    static_assert(std::is_floating_point<T>::value,
                  "newton_solve_2d_bounded requires floating-point T");

    using result_t = decltype(F(x0));
    static_assert(std::is_same<result_t, std::array<T, 2>>::value,
                  "F must return std::array<T, 2>");

    if (isnan(x0[0]) || isnan(x0[1])) throw std::runtime_error("x0 is nan in newton_solve_2d_bounded");
    if (isnan(x_min[0]) || isnan(x_min[1])) throw std::runtime_error("x_min is nan in newton_solve_2d_bounded");
    if (isnan(x_max[0]) || isnan(x_max[1])) throw std::runtime_error("x_max is nan in newton_solve_2d_bounded");

    // Clamp initial guess to valid range
    if (x0[0] < x_min[0] || x0[0] > x_max[0] || x0[1] < x_min[1] || x0[1] > x_max[1]) {
        std::cout << "WARNING: Initial guess x0 outside bounds of [x_min, x_max] in newton_solve_2d()! Clamping x0 to bounded interval.\n";
    }
    for (int j = 0; j < 2; ++j) {
        x0[j] = std::clamp(x0[j], x_min[j] + bound_margin, x_max[j] - bound_margin);
    }

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
        // Evaluate F(x0)
        auto fx = F(x0);

        if (isnan(fx[0]) || isnan(fx[1])) throw std::runtime_error("F(x) is nan in newton_solve_2d_bounded");

        T norm = std::sqrt(fx[0]*fx[0] + fx[1]*fx[1]);
        if (norm < tol) {
            if (dev) {
                std::cout << "norm < tol achieved for x0=(" << x0[0] << "," << x0[1] << "), terminating newton_solve_2d_bounded()!\n";
            }
            return x0;
        }

        // Compute Jacobian with adaptive one-sided derivatives
        std::array<std::array<T,2>,2> J{};
        for (int j = 0; j < 2; ++j) {
            auto xh = x0;
            
            // Choose forward or backward difference based on proximity to bounds
            bool use_forward = (x0[j] + h <= x_max[j] - bound_margin);
            bool use_backward = (x0[j] - h >= x_min[j] + bound_margin);
            
            if (use_forward) {
                // Forward difference
                xh[j] += h;
                auto fh = F(xh);
                for (int i = 0; i < 2; ++i) {
                    J[i][j] = (fh[i] - fx[i]) / h;
                }
            } else if (use_backward) {
                // Backward difference
                xh[j] -= h;
                auto fh = F(xh);
                for (int i = 0; i < 2; ++i) {
                    J[i][j] = (fx[i] - fh[i]) / h;
                }
            } else {
                // Stuck at boundary - use small central difference
                T h_small = std::min(x0[j] - x_min[j], x_max[j] - x0[j]) * 0.1;
                auto xh_plus = x0;
                auto xh_minus = x0;
                xh_plus[j] += h_small;
                xh_minus[j] -= h_small;
                auto fh_plus = F(xh_plus);
                auto fh_minus = F(xh_minus);
                for (int i = 0; i < 2; ++i) {
                    J[i][j] = (fh_plus[i] - fh_minus[i]) / (2.0 * h_small);
                }
            }
        }

        // Compute Newton step
        std::array<T,2> minus_fx = { -fx[0], -fx[1] };
        std::array<T,2> dx;
        try {
            dx = solve_linear_2x2(J, minus_fx);
        } catch (const std::runtime_error&) {
            // Singular Jacobian - try gradient descent step instead
            if (dev) {
                std::cout << "Warning: Singular Jacobian for x=(" << x0[0] << ", " << x0[1] << ")\n";
            }
            T grad_norm = std::sqrt(fx[0]*fx[0] + fx[1]*fx[1]);
            if (grad_norm < tol) return x0;
            
            dx[0] = -0.01 * fx[0] / grad_norm;
            dx[1] = -0.01 * fx[1] / grad_norm;
        }

        // Step limiting: ensure we don't go out of bounds
        T alpha = 1.0;  // Step size multiplier
        
        for (int j = 0; j < 2; ++j) {            
            // Calculate maximum allowed step to stay within bounds
            if (dx[j] > 0) {
                // Moving toward upper bound
                T max_step = (x_max[j] - bound_margin - x0[j]);
                if (dx[j] > max_step) {
                    T alpha_j = max_step / dx[j];
                    alpha = std::min(alpha, alpha_j);
                }
            } else if (dx[j] < 0) {
                // Moving toward lower bound
                T max_step = (x_min[j] + bound_margin - x0[j]);
                if (dx[j] < max_step) {
                    T alpha_j = max_step / dx[j];
                    alpha = std::min(alpha, alpha_j);
                }
            }
        }

        if (dev) {
            auto det_J = J[0][0] * J[1][1] - J[0][1] * J[1][0];
            std::cout << "iteration=" << iter << "\n"
                      << "x0=(" << x0[0] << "," << x0[1] << ")\n"
                      << "|J|=" << det_J << "\n"
                      << "dx=(" << dx[0] << "," << dx[1] << ")\n"
                      << "f(x0)=(" << fx[0] << "," << fx[1] << ")\n"
                      << "norm=" << norm << "\n"
                      << "alpha=" << alpha << "\n\n";
        }

        if (alpha < 1e-12) {
            throw std::runtime_error("newton_solve_2d_bounded failed (stuck at boundary)!");
        }

        // Apply limited step
        x0[0] += alpha * dx[0];
        x0[1] += alpha * dx[1];

        // Safety clamp (shouldn't be needed but just in case)
        for (int j = 0; j < 2; ++j) {
            x0[j] = std::clamp(x0[j], x_min[j] + bound_margin, x_max[j] - bound_margin);
        }

        // Check convergence
        T step_norm = std::sqrt(dx[0]*dx[0] + dx[1]*dx[1]) * alpha;
        if (step_norm < tol) {
            if (dev) {
                std::cout <<"dx[0]=" << dx[0] << ", dx[1]=" << dx[1] << ", alpha=" << alpha << "\n";
                std::cout << "step_norm=" << step_norm << " < tol=" << tol << " achieved, terminating newton_solve_2d_bounded()!\n";
            }
            return x0;
        }
    }

    throw std::runtime_error("newton_solve_2d_bounded did not converge");
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