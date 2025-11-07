// maths_ops.cpp
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <functional>
#include <cassert>
// #include <matplotlibcpp.h>

#include "ap.h"
#include "solvers.h"

#include "constants.hpp"
#include "maths_ops.hpp"

/*
TO DO:
- make vector class for vector arithmetic and include print_vector() function there so i can use vec.print()
- clean up use of std::vector<double>: any ICs should be std::vector<double>=std::array<double,3> and vectors of xi_vals should be std::vector<double>
*/

// void plot(std::vector<double> x_vals, std::vector<double> y_vals, const std::string& filename) {
//     namespace plt = matplotlibcpp;

//     plt::figure_size(800, 600);
//     plt::plot(x_vals, y_vals);
//     plt::grid(true);
//     plt::save("../" + filename);

//     return;
// }

std::vector<double> linspace(double start, double end, std::size_t num) {
    std::vector<double> result;

    if (num == 0)
        return result;
    if (num == 1) {
        result.push_back(start);
        return result;
    }

    double step = (end - start) / (num - 1);
    result.reserve(num);

    for (std::size_t i = 0; i < num; ++i) {
        result.push_back(start + i * step);
    }

    return result;
}

std::vector<double> logspace(double start, double stop, std::size_t num) {
    std::vector<double> result;
    double log_start = std::log10(start);  // Log of the start value
    double log_stop = std::log10(stop);    // Log of the stop value

    double step = (log_stop - log_start) / (num - 1);  // Step size in log space

    for (size_t i = 0; i < num; ++i) {
        double log_val = log_start + i * step;     // Calculate log value at step i
        result.push_back(std::pow(10, log_val));   // Convert back to linear space
    }

    return result;
}

// faster than std::pow
double power(double x, int exp) {
    if (exp == 0) return 1.0;  // x^0 = 1
    double result = x;
    for (int i = 1; i < std::abs(exp); ++i) {
        result *= x;
    }
    return (exp > 0) ? result : 1.0 / result;  // Handle negative exponents
}

std::string to_string_with_precision(double value, int precision) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

// simpson with non-uniform spacing - might be unstable?
// add error estimation somehow... (or maybe check against boost.math integration)
// could write check_stability() function that performs the integration a number of times with different number of pts to check when it converges
double simpson_integrate(const std::vector<double>& x, const std::vector<double>& y) {
    const size_t n = x.size();
    if (n != y.size()) {
        throw std::invalid_argument("x and y must be the same size");
    }
    if (n < 2) {
        throw std::invalid_argument("Need at least two points for integration");
    }

    double integral = 0.0;
    size_t limit = (n % 2 == 0) ? n - 1 : n;  // If even, stop at n-1 to apply Simpson's

    for (size_t i = 0; i + 2 < limit; i += 2) {
        double h = x[i + 2] - x[i];
        double h1 = x[i + 1] - x[i];
        double h2 = x[i + 2] - x[i + 1];

        // If spacing is non-uniform, adjust accordingly
        if (std::abs(h1 - h2) > 1e-8) {
            // fallback to composite trapezoid for irregular spacing
            integral += 0.5 * (x[i + 1] - x[i]) * (y[i] + y[i + 1]);
            integral += 0.5 * (x[i + 2] - x[i + 1]) * (y[i + 1] + y[i + 2]);
        } else {
            integral += (h / 6.0) * (y[i] + 4 * y[i + 1] + y[i + 2]);
        }
    }

    // Trapezoid on the last interval if n is even
    if (n % 2 == 0) {
        double h = x[n - 1] - x[n - 2];
        integral += 0.5 * h * (y[n - 2] + y[n - 1]);
    }

    return integral;
}

double simpson_2d_integrate(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f) {
    const size_t nx = x.size();
    const size_t ny = y.size();

    if (f.size() != ny || f[0].size() != nx) {
        throw std::invalid_argument("Size of f must be [y.size()][x.size()]");
    }

    double total = 0.0;

    size_t nx_lim = (nx % 2 == 0) ? nx - 1 : nx;
    size_t ny_lim = (ny % 2 == 0) ? ny - 1 : ny;

    for (size_t j = 0; j + 2 < ny_lim; j += 2) {
        double hy1 = y[j + 1] - y[j];
        double hy2 = y[j + 2] - y[j + 1];
        double hy = y[j + 2] - y[j];

        for (size_t i = 0; i + 2 < nx_lim; i += 2) {
            double hx1 = x[i + 1] - x[i];
            double hx2 = x[i + 2] - x[i + 1];
            double hx = x[i + 2] - x[i];

            if (std::abs(hx1 - hx2) > 1e-8 || std::abs(hy1 - hy2) > 1e-8) {
                // Use composite trapezoidal if spacing is non-uniform
                double area =
                0.25 * (x[i + 1] - x[i]) * (y[j + 1] - y[j]) * (f[j][i] + f[j + 1][i] + f[j][i + 1] + f[j + 1][i + 1]) +
                0.25 * (x[i + 2] - x[i + 1]) * (y[j + 1] - y[j]) * (f[j][i + 1] + f[j + 1][i + 1] + f[j][i + 2] + f[j + 1][i + 2]) +
                0.25 * (x[i + 1] - x[i]) * (y[j + 2] - y[j + 1]) * (f[j + 1][i] + f[j + 2][i] + f[j + 1][i + 1] + f[j + 2][i + 1]) +
                0.25 * (x[i + 2] - x[i + 1]) * (y[j + 2] - y[j + 1]) * (f[j + 1][i + 1] + f[j + 2][i + 1] + f[j + 1][i + 2] + f[j + 2][i + 2]);
                total += area;
            } else {
                // Use 2D Simpson's rule
                total += (hx * hy / 36.0) * (
                f[j][i] + 4 * f[j][i + 1] + f[j][i + 2] +
                4 * (f[j + 1][i] + 4 * f[j + 1][i + 1] + f[j + 1][i + 2]) +
                f[j + 2][i] + 4 * f[j + 2][i + 1] + f[j + 2][i + 2]
                );
            }
        }
    }

    // Handle remaining strip in x if nx is even
    if (nx % 2 == 0) {
        for (size_t j = 0; j + 1 < ny; ++j) {
            double hy = y[j + 1] - y[j];
            double hx = x[nx - 1] - x[nx - 2];
            total += 0.25 * hx * hy * (
            f[j][nx - 2] + f[j + 1][nx - 2] +
            f[j][nx - 1] + f[j + 1][nx - 1]);
        }
    }

    // Handle remaining strip in y if ny is even
    if (ny % 2 == 0) {
        for (size_t i = 0; i + 1 < nx; ++i) {
            double hx = x[i + 1] - x[i];
            double hy = y[ny - 1] - y[ny - 2];
            total += 0.25 * hx * hy * (
            f[ny - 2][i] + f[ny - 2][i + 1] +
            f[ny - 1][i] + f[ny - 1][i + 1]);
        }
    }

    return total;
}

double simpson_2d_integrate_flat(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& f_flat) {
    const size_t nx = x.size();
    const size_t ny = y.size();

    if (f_flat.size() != nx * ny) {
        throw std::invalid_argument("Size of f_flat must be x.size() * y.size()");
    }

    auto at = [&](size_t j, size_t i) -> double {
        return f_flat[j * nx + i];
    };

    double total = 0.0;

    size_t nx_lim = (nx % 2 == 0) ? nx - 1 : nx;
    size_t ny_lim = (ny % 2 == 0) ? ny - 1 : ny;

    for (size_t j = 0; j + 2 < ny_lim; j += 2) {
        double hy1 = y[j + 1] - y[j];
        double hy2 = y[j + 2] - y[j + 1];
        double hy = y[j + 2] - y[j];

        for (size_t i = 0; i + 2 < nx_lim; i += 2) {
            double hx1 = x[i + 1] - x[i];
            double hx2 = x[i + 2] - x[i + 1];
            double hx = x[i + 2] - x[i];

            if (std::abs(hx1 - hx2) > 1e-8 || std::abs(hy1 - hy2) > 1e-8) {
                // Composite trapezoidal rule
                total += 0.25 * hx1 * hy1 * (at(j, i) + at(j + 1, i) + at(j, i + 1) + at(j + 1, i + 1));
                total += 0.25 * hx2 * hy1 * (at(j, i + 1) + at(j + 1, i + 1) + at(j, i + 2) + at(j + 1, i + 2));
                total += 0.25 * hx1 * hy2 * (at(j + 1, i) + at(j + 2, i) + at(j + 1, i + 1) + at(j + 2, i + 1));
                total += 0.25 * hx2 * hy2 * (at(j + 1, i + 1) + at(j + 2, i + 1) + at(j + 1, i + 2) + at(j + 2, i + 2));
            } else {
                // 2D Simpson's rule
                total += (hx * hy / 36.0) * (
                    at(j, i) + 4 * at(j, i + 1) + at(j, i + 2) +
                    4 * (at(j + 1, i) + 4 * at(j + 1, i + 1) + at(j + 1, i + 2)) +
                    at(j + 2, i) + 4 * at(j + 2, i + 1) + at(j + 2, i + 2)
                );
            }
        }
    }

    // Handle right edge strip in x
    if (nx % 2 == 0) {
        for (size_t j = 0; j + 1 < ny; ++j) {
            double hx = x[nx - 1] - x[nx - 2];
            double hy = y[j + 1] - y[j];
            total += 0.25 * hx * hy * (
                at(j, nx - 2) + at(j + 1, nx - 2) +
                at(j, nx - 1) + at(j + 1, nx - 1)
            );
        }
    }

    // Handle top edge strip in y
    if (ny % 2 == 0) {
        for (size_t i = 0; i + 1 < nx; ++i) {
            double hx = x[i + 1] - x[i];
            double hy = y[ny - 1] - y[ny - 2];
            total += 0.25 * hx * hy * (
                at(ny - 2, i) + at(ny - 2, i + 1) +
                at(ny - 1, i) + at(ny - 1, i + 1)
            );
        }
    }

    return total;
}

void precompute_1d_weights(
    const std::vector<double>& coords,
    std::vector<std::vector<double>>& weights,
    std::vector<double>& intervals)
{
    const size_t n = coords.size();
    if (n < 3) {
        throw std::invalid_argument("Grid must have at least 3 points for Simpson's rule");
    }

    weights.resize(n - 2);
    intervals.resize(n - 2);

    for (size_t i = 0; i + 2 < n; ++i) {
        const double x0 = coords[i];
        const double x1 = coords[i+1];
        const double x2 = coords[i+2];

        const double hx0 = x1 - x0;
        const double hx1 = x2 - x1;
        const double denom = hx0 * hx1 * (hx0 + hx1);

        // Compute weights for points i, i+1, i+2
        double A0 = -hx1 * hx1 / denom;
        double A1 = (hx1 * hx1 - hx0 * hx0) / denom;
        double A2 = hx0 * hx0 / denom;

        weights[i] = {A0, A1, A2};
        intervals[i] = x2 - x0;
    }
}

// define this without function call to 1d version for speed
SimpsonWeights2D precompute_simpson_weights_2d(
    const std::vector<double>& x,
    const std::vector<double>& y)
{
    SimpsonWeights2D w;

    precompute_1d_weights(x, w.Ax_weights, w.dx);
    precompute_1d_weights(y, w.Ay_weights, w.dy);

    return w;
}

double simpson_2d_nonuniform_flat_weighted(
    const std::vector<double>& x,                        // full x vector, size nx
    const std::vector<double>& y,                        // full y vector, size ny
    const std::vector<double>& f_flat,                   // flattened f array, size nx * ny
    const std::vector<std::vector<double>>& Ax_weights,  // size (nx-2)/2 x 3
    const std::vector<std::vector<double>>& Ay_weights,  // size (ny-2)/2 x 3
    const std::vector<double>& dx,                       // size (nx-2)/2
    const std::vector<double>& dy                        // size (ny-2)/2
) {
    const size_t nx = x.size(); // pass this into integrator for slight speedup
    const size_t ny = y.size();

    if (f_flat.size() != nx * ny)
        throw std::invalid_argument("f_flat size must equal nx * ny");

    auto f = [&](size_t j, size_t i) -> double {
        return f_flat[j * nx + i];
    };

    double total = 0.0;

    // Bulk integration with precomputed Simpson weights
    for (size_t i = 0; i + 2 < nx; i += 2) {
        size_t wx_idx = i / 2;
        for (size_t j = 0; j + 2 < ny; j += 2) {
            size_t wy_idx = j / 2;

            double local = 0.0;
            for (int jj = 0; jj < 3; ++jj) {
                for (int ii = 0; ii < 3; ++ii) {
                    local += Ax_weights[wx_idx][ii] * Ay_weights[wy_idx][jj] * f(j + jj, i + ii);
                }
            }
            total += dx[wx_idx] * dy[wy_idx] * local / 4.0;
        }
    }

    // Handle last row (if ny even) with trapezoidal rule along y
    if (ny % 2 == 0) {
        size_t last_y = ny - 2;
        double hy = y[ny - 1] - y[last_y];

        for (size_t i = 0; i + 1 < nx; ++i) {
            double hx = x[i + 1] - x[i];

            total += 0.25 * hx * hy * (
                f(last_y, i) + f(last_y + 1, i) +
                f(last_y, i + 1) + f(last_y + 1, i + 1)
            );
        }
    }

    // Handle last column (if nx even) with trapezoidal rule along x
    if (nx % 2 == 0) {
        size_t last_x = nx - 2;
        double hx = x[nx - 1] - x[last_x];

        for (size_t j = 0; j + 1 < ny; ++j) {
            double hy = y[j + 1] - y[j];

            total += 0.25 * hx * hy * (
                f(j, last_x) + f(j + 1, last_x) +
                f(j, last_x + 1) + f(j + 1, last_x + 1)
            );
        }
    }

    return total;
}

// Im(Si(x))=0 for real x
// Im(Ci(x))=pi for x<0 (0 otherwise)
double Si(double x) {
    if (x == 0.0) return 0.0;

    auto sin_integrand = [](double t) {
        return (t == 0.0) ? 1.0 : std::sin(t) / t;
    };

    const int nt = 200; // integration steps
    std::vector<double> t_vals = linspace(0.0, x, nt);
    std::vector<double> sin_integrand_vals(nt);
    
    for (int j = 0; j < nt; j++) {
        const auto t = t_vals[j];
        sin_integrand_vals[j] = sin_integrand(t);
    }

    return simpson_integrate(t_vals, sin_integrand_vals);
}

double Ci(double x) {
    if (x == 0.0) return -INFINITY;

    auto cos_integrand = [](double t) {
        if (std::abs(t) < 1e-4) return -0.5 * t; // Taylor exp near zero (otherwise simpson integration breaks)
        return (std::cos(t) - 1.0) / t;
    };

    const int nt = 200; // integration steps
    std::vector<double> t_vals = linspace(0.0, x, nt);
    std::vector<double> cos_integrand_vals(nt);
    
    for (int j = 0; j < nt; j++) {
        const auto t = t_vals[j];
        cos_integrand_vals[j] = cos_integrand(t);
    }

    return gamma_euler + std::log(abs(x)) + simpson_integrate(t_vals, cos_integrand_vals);
}

std::pair<double, double> SiCi(double x, const size_t n) {
    if (x == 0.0) return {0.0, -INFINITY};

    auto sin_integrand = [](double t) {
        return (t == 0.0) ? 1.0 : std::sin(t) / t;
    };

    auto cos_integrand = [](double t) {
        if (std::abs(t) < 1e-4) return -0.5 * t; // Taylor exp near zero (otherwise simpson integration breaks)
        return (std::cos(t) - 1.0) / t;
    };

    std::vector<double> t_vals = linspace(0.0, x, n);
    std::vector<double> sin_integrand_vals(n);
    std::vector<double> cos_integrand_vals(n);
    
    for (size_t j = 0; j < n; j++) {
        const auto t = t_vals[j];
        sin_integrand_vals[j] = sin_integrand(t);
        cos_integrand_vals[j] = cos_integrand(t);
    }

    const auto sin_int = simpson_integrate(t_vals, sin_integrand_vals);
    const auto cos_int = gamma_euler + std::log(abs(x)) + simpson_integrate(t_vals, cos_integrand_vals);

    return {sin_int, cos_int};
}

std::vector<double> dSiCi(double x, double y, const size_t n) {
    // calculates dSi = Si(x) - Si(y) and dCi = Ci(x) - Ci(y)
    auto sin_integrand = [] (double t) {
        return std::sin(t) / t;
    };
    auto cos_integrand = [] (double t) {
        return (std::cos(t) - 1.0) / t;
    };

    // fill integrand vals
    const std::vector<double> t_vals = linspace(y, x, n);
    std::vector<double> sin_integrand_vals(n), cos_integrand_vals(n);
    for (size_t i = 0; i < n; i++) {
        const auto t = t_vals[i];
        sin_integrand_vals[i] = sin_integrand(t);
        cos_integrand_vals[i] = cos_integrand(t);
    }

    // integrate
    const auto dSi = simpson_integrate(t_vals, sin_integrand_vals);
    const auto dCi = std::log(x/y) + simpson_integrate(t_vals, cos_integrand_vals);

    return {dSi, dCi};
}

// write my own/make this better
// make it suitable to pass in just dvdxi (rather than only {dvdxi, dwdxi})
// std::pair<prof_type, std::vector<state_type>> rk4_solver(const deriv_func& dydx, double x0, double xf, const state_type& y0, size_t n) {
//     assert(n >= 2 && "Number of steps must be at least 2.");

//     const double h = (xf - x0) / static_cast<double>(n - 1);
//     prof_type x_vals(n);
//     std::vector<state_type> y_vals(n);

//     x_vals[0] = x0;
//     y_vals[0] = y0;

//     double x = x0;
//     state_type y = y0;

//     for (size_t i = 1; i < n; ++i) {
//         const auto k1 = dydx(x, y);

//         state_type y_tmp;
//         for (int j = 0; j < 3; ++j) y_tmp[j] = y[j] + 0.5*h*k1[j];
//         const auto k2 = dydx(x + 0.5*h, y_tmp);

//         for (int j = 0; j < 3; ++j) y_tmp[j] = y[j] + 0.5*h*k2[j];
//         const auto k3 = dydx(x + 0.5*h, y_tmp);

//         for (int j = 0; j < 3; ++j) y_tmp[j] = y[j] + h*k3[j];
//         const auto k4 = dydx(x + h, y_tmp);

//         for (int j = 0; j < 3; ++j) {
//             y[j] += (h/6.0) * (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);
//         }

//         x += h;
//         x_vals[i] = x;
//         y_vals[i] = y;
//     }

//     return {x_vals, y_vals};
// }

// modified bisection method only to be used for functions with exactly one root!
double root_finder(std::function<double(double)> f, double a, double b, double tol, int max_iter) {
    double fa = f(a);
    double fb = f(b);

    // std::cout << "f(a=" << a <<")=" << fa << ", f(b=" << b << ")=" << fb << std::endl;

    if (fa * fb > 0.0) {
        throw std::runtime_error("Bisection method interval not bracketed!");
    }

    if (!std::isfinite(fa) || !std::isfinite(fb)) {
        throw std::runtime_error("f(a) or f(b) is not finite.");
    }

    for (int i = 0; i < max_iter; ++i) {
        double c = 0.5 * (a + b);
        // std::cout << "c=" << c << "\n";
        double fc = f(c);

        if (!std::isfinite(fc)) {
            // throw std::runtime_error("f(c) is not finite during bisection.");
        }

        // Check convergence
        // since interval is not necessarily bracketed, need to remove second condition
        // otherwise it thinks it successfully found the root once it has gone through 
        // the entire interval
        if (std::abs(fc) < tol) {
            return c;
        }

        if (fa * fc < 0.0) {
        // if (fc > 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    // std::cout << "Bisection method failed, using smallest root algorithm.\n";
    // return find_smallest_root(f, a, b);
    throw std::runtime_error("Bisection method did not converge.");
}

// std::array<double, 2> newton_solve_2d(const std::function<std::array<double, 2>(std::array<double, 2>)>& F, std::array<double, 2> x0, double tol, int max_iter, double h) {
//     if (x0.size() != 2)
//         throw std::invalid_argument("newton_solve requires exactly 2 variables.");

//     // Inline solver for 2×2 linear systems
//     auto solve_linear_2x2 = [](const std::array<std::array<double,2>,2>& A,
//                                const std::array<double,2>& b) {
//         double det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
//         if (std::fabs(det) < 1e-14)
//             throw std::runtime_error("Jacobian is singular!");
//         std::array<double,2> x;
//         x[0] = ( b[0]*A[1][1] - b[1]*A[0][1]) / det;
//         x[1] = ( A[0][0]*b[1] - A[1][0]*b[0]) / det;
//         return x;
//     };

//     for (int iter = 0; iter < max_iter; ++iter) {
//         // Evaluate F(x) once
//         auto fx = F(x0);
//         if (fx.size() != 2)
//             throw std::invalid_argument("Function F must return exactly 2 values.");

//         double norm = std::sqrt(fx[0]*fx[0] + fx[1]*fx[1]);
//         if (norm < tol) return x0;

//         // Numerical Jacobian with reuse of F(x0)
//         std::array<std::array<double,2>,2> J{};
//         for (int j = 0; j < 2; ++j) {
//             auto xh = x0;
//             xh[j] += h;
//             auto fh = F(xh);
//             for (int i = 0; i < 2; ++i) {
//                 J[i][j] = (fh[i] - fx[i]) / h;
//             }
//         }

//         std::array<double,2> minus_fx = { -fx[0], -fx[1] };
//         auto dx = solve_linear_2x2(J, minus_fx);

//         x0[0] += dx[0];
//         x0[1] += dx[1];

//         if (std::sqrt(dx[0]*dx[0] + dx[1]*dx[1]) < tol)
//             return x0;
//     }

//     throw std::runtime_error("Newton's method solver did not converge");
// }

double golden_section_minimize(std::function<double(double)> f, double a, double b, double tol, int max_iter) {
    const double gr = (std::sqrt(5.0) - 1.0) / 2.0; // 0.618...
    double c = b - gr * (b - a);
    double d = a + gr * (b - a);
    double fc = f(c);
    double fd = f(d);

    int iter = 0;
    while ((b - a) > tol && iter++ < max_iter) {
        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = f(d);
        }
    }
    double xmid = 0.5*(a + b);
    return xmid;
}

// unused
double brent_minimize(std::function<double(double)> f, double a, double b, double tol, int max_iter) {
    const double phi = (3.0 - std::sqrt(5.0)) / 2.0; // ~0.3819660
    double x = a + phi * (b - a);
    double w = x;
    double v = x;
    double fx = f(x);
    double fw = fx;
    double fv = fx;

    double d = 0.0; // step size
    double e = 0.0; // distance moved on step before last

    for (int iter = 0; iter < max_iter; ++iter) {
        double m = 0.5 * (a + b);
        double tol1 = tol * std::fabs(x) + 1e-12;
        double tol2 = 2.0 * tol1;

        // Check convergence
        if (std::fabs(x - m) <= tol2 - 0.5 * (b - a)) {
            return x;
        }

        bool accept_parabolic = false;
        double p = 0.0, q = 0.0, r = 0.0;

        if (std::fabs(e) > tol1) {
            // Fit parabola
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p;
            q = std::fabs(q);

            if (std::fabs(p) < std::fabs(0.5 * q * e) &&
                p > q * (a - x) &&
                p < q * (b - x)) {
                // Parabolic step
                d = p / q;
                accept_parabolic = true;
            }
        }

        if (!accept_parabolic) {
            // Golden-section step
            if (x < m) e = b - x;
            else e = a - x;
            d = phi * e;
        }

        double u = x + (std::fabs(d) >= tol1 ? d : (d > 0 ? tol1 : -tol1));
        double fu = f(u);

        // Update points
        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w; fv = fw;
            w = x; fw = fx;
            x = u; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }

    // If max_iter reached
    return x;
}

alglib::real_1d_array vector_to_real_1d_array(const std::vector<double>& vec) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << vec[i];
        if (i != vec.size() - 1) oss << ",";
    }
    oss << "]";
    
    alglib::real_1d_array result;
    result = oss.str().c_str();  // real_1d_array supports string assignment
    return result;
}

std::array<double, 2> find_bracket(const std::function<double(double)>& residual_func, double a, double b) {
    constexpr int NSAMPLES = 200;      // coarse scan resolution
    constexpr double PENALTY = 1e+5;      // anything larger than realistic residuals
    constexpr double PAD_MULT = 2.0;      // how many sample spacings to pad the bracket

    // --- 1. Sample the interval ------------------------------------------
    std::vector<double> xs(NSAMPLES + 1), ys(NSAMPLES + 1);
    for (int i = 0; i <= NSAMPLES; ++i) {
        double x = a + (b - a) * (double(i) / NSAMPLES);
        xs[i] = x;
        double y = residual_func(x);
        ys[i] = (std::isfinite(y) ? y : PENALTY * 10.0); // ensure numeric
    }

    // --- 2. Identify contiguous finite segments --------------------------
    struct Segment { int i0, i1; };
    std::vector<Segment> segs;
    int i = 0;
    while (i <= NSAMPLES) {
        while (i <= NSAMPLES && ys[i] >= PENALTY) ++i;
        if (i > NSAMPLES) break;
        int start = i;
        while (i <= NSAMPLES && ys[i] < PENALTY) ++i;
        int end = i - 1;
        segs.push_back({start, end});
    }

    if (segs.empty()) {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
    }

    // --- 3. Find the segment with the lowest sampled residual ------------
    int idx_min = std::min_element(ys.begin(), ys.end()) - ys.begin();
    int chosen_seg = -1;
    for (int k = 0; k < (int)segs.size(); ++k) {
        if (idx_min >= segs[k].i0 && idx_min <= segs[k].i1) {
            chosen_seg = k;
            break;
        }
    }
    if (chosen_seg == -1) {
        // fallback: pick segment containing smallest finite residual
        double best_val = std::numeric_limits<double>::infinity();
        for (int k = 0; k < (int)segs.size(); ++k) {
            for (int j = segs[k].i0; j <= segs[k].i1; ++j) {
                if (ys[j] < best_val) { best_val = ys[j]; chosen_seg = k; }
            }
        }
    }

    const auto& S = segs[chosen_seg];

    // --- 4. Construct bracket from that segment --------------------------
    double dx  = (b - a) / NSAMPLES;
    double left  = std::max(a, xs[S.i0] - PAD_MULT * dx);
    double right = std::min(b, xs[S.i1] + PAD_MULT * dx);

    return {left, right};
}