// maths_ops.hpp
#ifndef INCLUDE_MATHS_OPS_HPP_H
#define INCLUDE_MATHS_OPS_HPP_H

/**
 * @file maths_ops.hpp
 * @brief Mathematical utilities including a custom vector class and cubic spline interpolation.
 *
 * Provides a template `vec<T>` class with arithmetic operations, dot product, and norm.
 * Includes a `CubicSpline<T>` class for 1D spline interpolation and utility functions like `linspace`.
 */


#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <initializer_list>
#include <functional>
#include <streambuf>
#include <limits>
#include <type_traits>
#include <array>
#include <complex>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include <gsl/gsl_integration.h>

#include "interpolation.h"

#include "ap.h"

/**
 * @brief Cubic spline interpolation class.
 * 
 * @tparam T Type for input/output values (typically double).
 */
template <typename T>
class CubicSpline {
  static_assert(std::is_floating_point<T>::value, "CubicSpline only supports floating-point types.");

  public:
    CubicSpline();  // Default constructor
    CubicSpline(const std::vector<T>& x, const std::vector<T>& y);  // Construct and compute spline

    void build(const std::vector<T>& x, const std::vector<T>& y); // Build splines
    bool is_initialised() const {return initialised_; } // check if spline is initialised
    void check_convergence() const; // check convergence

    T operator()(T xi) const;  // Evaluate spline at point xi

    // Scalar operations from the right
    CubicSpline<T> operator+(T scalar) const;
    CubicSpline<T> operator-(T scalar) const;
    CubicSpline<T> operator*(T scalar) const;
    CubicSpline<T> operator/(T scalar) const;

    // Scalar operations from the left
    template <typename U>
    friend CubicSpline<U> operator+(U scalar, const CubicSpline<U>& spline);

    template <typename U>
    friend CubicSpline<U> operator-(U scalar, const CubicSpline<U>& spline);

    template <typename U>
    friend CubicSpline<U> operator*(U scalar, const CubicSpline<U>& spline);

  private:
    std::vector<T> x_, y_;          // Input data
    std::vector<T> h_;              // Interval widths
    std::vector<T> a_, b_, c_, d_;  // Spline coefficients

    const T tol_ = static_cast<T>(1e-4); // convergence tolerance
    bool initialised_ = false;

    bool is_strictly_monotonic(const std::vector<T>& x) const; // Check monotonicity
};
#include "CubicSpline.tpp"

/**
 * @brief Creates a linearly spaced vector between `start` and `end` (C++ equivalent of python's linspace function).
 *
 * @param start Initial value.
 * @param end Final value.
 * @param num Number of points.
 * 
 * @return std::vector<double> Linearly spaced values.
 */
std::vector<double> linspace(double start, double end, std::size_t num=100);

std::vector<double> logspace(double start, double stop, std::size_t num=100);

/**
 * @brief Computes x raised to an integer exponent.
 * 
 * @param x Base value.
 * @param exp Integer exponent.
 * 
 * @return double Result of x^exp.
 * 
 */
double power(double x, int exp);

std::string to_string_with_precision(double value, int precision = 2);

double simpson_integrate(const std::vector<double>& x, const std::vector<double>& y);
double simpson_2d_integrate(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& f);
double simpson_2d_integrate_flat(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& f_flat);

struct SimpsonWeights2D {
    std::vector<std::vector<double>> Ax_weights; // size: (nx-2) x 3
    std::vector<std::vector<double>> Ay_weights; // size: (ny-2) x 3
    std::vector<double> dx;  // size: nx-2
    std::vector<double> dy;  // size: ny-2
};

void precompute_1d_weights(
    const std::vector<double>& coords,
    std::vector<std::vector<double>>& weights,
    std::vector<double>& intervals);

SimpsonWeights2D precompute_simpson_weights_2d(
    const std::vector<double>& x,
    const std::vector<double>& y);

double simpson_2d_nonuniform_flat_weighted(
    const std::vector<double>& x,                        // full x vector, size nx
    const std::vector<double>& y,                        // full y vector, size ny
    const std::vector<double>& f_flat,                   // flattened f array, size nx * ny
    const std::vector<std::vector<double>>& Ax_weights,  // size (nx-2)/2 x 3
    const std::vector<std::vector<double>>& Ay_weights,  // size (ny-2)/2 x 3
    const std::vector<double>& dx,                       // size (nx-2)/2
    const std::vector<double>& dy                        // size (ny-2)/2
);

double Si(double x);
double Ci(double x);
std::pair<double, double> SiCi(double x, const size_t n=1000);
std::vector<double> dSiCi(double x, double y, const size_t n);

double root_finder(std::function<double(double)> f, double a, double b, double tol = 1e-8, int max_iter = 100);
double root_finder_new(std::function<double(double)> f, double a, double b, double tol = 1e-8, int max_iter = 100);

double golden_section_minimize(std::function<double(double)> f, double a, double b, double tol = 1e-8, int max_iter = 100);
double brent_minimize(std::function<double(double)> f, double a, double b, double tol = 1e-8, int max_iter = 100);

alglib::real_1d_array vector_to_real_1d_array(const std::vector<double>& vec);

std::array<double, 2> find_bracket(const std::function<double(double)>& residual_func, double a, double b);

// solvers
template <typename T, typename State, typename Func>
std::pair<std::vector<T>, std::vector<State>> rk4_solver(const Func& dydx, T x0, T xf, const State& y0, size_t n=1000);

template <typename T, typename Func>
T newton_solve_1d(const Func& F, T x0, T tol=1e-8, int max_iter=100, T h=1e-8);

template <typename T, typename Func>
std::array<T,2> newton_solve_2d(const Func& F, std::array<T,2> x0, T tol=1e-8, int max_iter=100, T h=1e-8);

namespace Levin {

// Build Chebyshev differentiation matrix (Trefethen-style) on N nodes (N>=2).
// Nodes returned are physical nodes on [a,b], ordered k=0..N-1 where
// x(0) = b, x(N-1) = a (descending with cos).
inline void chebyshev_nodes_and_D(int N, double a, double b,
                                  Eigen::VectorXd &x_out,
                                  Eigen::MatrixXd &D_out)
{
    if (N < 2) throw std::invalid_argument("N must be >= 2");
    x_out.resize(N);
    D_out.resize(N, N);

    // Chebyshev nodes in t in [-1,1]: t_k = cos(pi * k / (N-1)), k=0..N-1
    // Map to x in [a,b]: x = 0.5*(a+b) + 0.5*(b-a)*t
    const double mid = 0.5 * (a + b);
    const double half = 0.5 * (b - a);

    Eigen::VectorXd t(N);
    for (int k = 0; k < N; ++k) {
        t[k] = std::cos(M_PI * k / double(N - 1));
        x_out[k] = mid + half * t[k];
    }

    // Build differentiation matrix D_t for derivative with respect to t (in [-1,1])
    // Then convert to derivative with respect to x: d/dx = (2/(b-a)) d/dt
    Eigen::VectorXd c(N);
    for (int k = 0; k < N; ++k) {
        if (k == 0 || k == N - 1) c[k] = 2.0;
        else c[k] = 1.0;
        if ((k & 1) != 0) c[k] = -c[k];
    }

    Eigen::MatrixXd D_t(N, N);
    D_t.setZero();
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            double denom = t[i] - t[j];
            D_t(i, j) = (c[i] / c[j]) * (std::pow(-1.0, i + j) / denom);
        }
    }

    // diagonal entries
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            sum += D_t(i, j);
        }
        D_t(i, i) = -sum;
    }

    // scale to physical x
    const double scale = 2.0 / (b - a);
    D_out = scale * D_t;
}

/*
    DISCLAIMER :
    ChatGPT 5.0 created the code below. 
    I have no idea if it is correct or efficient, but it seems to work.
    I have no idea how it works.

    I don't know if I should be proud or worried.
*/

// Levin integrator for sine kernel
// Returns integral I = ∫_a^b f(x) sin(chi * x) dx
// spline : alglib::spline1dinterpolant (must be built before calling)
// N_collocation : number of Chebyshev collocation nodes (suggest 24..64)
// tol_print : if true prints condition diagnostics
inline std::pair<double, double> levin_integrate(const alglib::spline1dinterpolant &spline,
                                 double chi, double a, double b,
                                 int N_collocation = 40,
                                 bool tol_print = false)
{
    if (a >= b) throw std::invalid_argument("a must be < b");
    if (N_collocation < 2) throw std::invalid_argument("need N_collocation >= 2");

    // 1) Build nodes and differentiation matrix
    Eigen::VectorXd x;
    Eigen::MatrixXd D;
    chebyshev_nodes_and_D(N_collocation, a, b, x, D);

    // 2) Evaluate f at nodes (real)
    Eigen::VectorXcd fvec(N_collocation);
    for (int k = 0; k < N_collocation; ++k) {
        double xv = x[k];
        double fv = alglib::spline1dcalc(spline, xv);
        fvec[k] = std::complex<double>(fv, 0.0);
    }

    // 3) Build complex linear system: (D + i*chi*I) u = f
    Eigen::MatrixXcd A = D.cast<std::complex<double> >();
    std::complex<double> imagchi(0.0, chi);
    for (int i = 0; i < N_collocation; ++i) A(i, i) += imagchi;

    // 4) Solve system (use PartialPivLU for robustness)
    Eigen::PartialPivLU<Eigen::MatrixXcd> solver(A);
    Eigen::VectorXcd u = solver.solve(fvec);

    // optional diagnostics: condition estimation via matrix inverse norm ≈ ||U||*||inv(U)||
    if (tol_print) {
        double anorm = A.cwiseAbs().sum(); // crude norm
        // attempt to compute a modest condition estimate via solving A x = e_j is expensive;
        std::cerr << "[Levin] N=" << N_collocation << " matrix rowsum-norm ~ " << anorm << "\n";
    }

    // 5) nodes ordering: recall x[0] = b, x[N-1] = a because of cos(pi*k/(N-1))
    std::complex<double> u_b = u[0];                // at x = b
    std::complex<double> u_a = u[N_collocation-1]; // at x = a

    // 6) compute integral: Im( u(b) e^{i chi b} - u(a) e^{i chi a} )
    std::complex<double> eb = std::exp(std::complex<double>(0.0, chi * x[0])); // e^{i chi b}
    std::complex<double> ea = std::exp(std::complex<double>(0.0, chi * x[N_collocation-1])); // e^{i chi a}
    std::complex<double> val = u_b * eb - u_a * ea;
    
    double integral_sin = std::imag(val);
    double integral_cos = std::real(val);

    return {integral_sin, integral_cos};
}
    
} // namespace Levin

#include "solvers.tpp"

#endif // INCLUDE_MATHS_OPS_HPP_H
