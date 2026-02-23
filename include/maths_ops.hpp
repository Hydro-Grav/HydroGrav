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
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>

// #include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include "interpolation.h"

#include "ap.h"

// for EoS with p(T), e(T) ~ O(1e+9), need residuals to converge to at least O(1e+10)
const static double minimiser_tol = 1e-10;

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
 * @brief Filon-type integrator for highly oscillatory integrals.
 */
class FilonQuadrature {
public:
    /**
     * @brief Construct a Filon-type integrator.
     * @param n_points Order of polynomial approximation (default: 8)
     */
    FilonQuadrature(int n_points = 8);
    
    /**
     * @brief Integrate f(x) * sin(omega * x) from a to b.
     * @param f Function to integrate (should accept double and return double)
     * @param omega Oscillatory frequency parameter
     * @param a Lower integration limit
     * @param b Upper integration limit
     * @return Integral value
     */
    double integrate_sin(const std::function<double(double)>& f, double omega, double a, double b);
    
    /**
     * @brief Integrate f(x) * cos(omega * x) from a to b.
     * @param f Function to integrate (should accept double and return double)
     * @param omega Oscillatory frequency parameter
     * @param a Lower integration limit
     * @param b Upper integration limit
     * @return Integral value
     */
    double integrate_cos(const std::function<double(double)>& f, double omega, double a, double b);
    
private:
    int n_points_;
    
    /**
     * @brief Integrate a single subinterval using Filon-type quadrature.
     * @param f_vals Function values at quadrature points
     * @param x_vals Quadrature point locations
     * @param omega Oscillatory frequency
     * @param a Lower limit
     * @param b Upper limit
     * @param use_sin True for sin integral, false for cos integral
     */
    double filon_integrate_interval(
        const std::vector<double>& f_vals,
        const std::vector<double>& x_vals,
        double omega,
        double a,
        double b,
        bool use_sin);
};

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
T newton_solve_1d_bounded(
    const Func& f,
    T x0,
    T x_min,
    T x_max,
    T tol = 1e-8,
    int max_iter = 100,
    T h = 1e-8,
    T bound_margin = 1e-6
);

template <typename T, typename Func>
std::array<T,2> newton_solve_2d(const Func& F, std::array<T,2> x0, T tol=1e-8, int max_iter=100, T h=1e-8);

template <typename T, typename Func>
std::array<T,2> newton_solve_2d_bounded(
    const Func& F,
    std::array<T,2> x0,
    std::array<T,2> x_min,  // Lower bounds
    std::array<T,2> x_max,  // Upper bounds
    T tol = 1e-8,
    int max_iter = 100,
    T h = 1e-8,
    T bound_margin = 1e-6,  // Stay this far from boundaries
    bool dev = false
);

std::array<double, 2> bisection_2d(
    std::function<std::array<double, 2>(std::array<double, 2>)> F,
    std::array<double, 2> lower,
    std::array<double, 2> upper,
    double tol=1e-8,
    int max_iter=50,
    bool dev=false);

#include "solvers.tpp"

double wasserstein_distance_1d(std::vector<double> u_values, 
                                std::vector<double> u_weights,
                                std::vector<double> v_values, 
                                std::vector<double> v_weights);

double L1_norm(const std::vector<double>& x_values,
               const std::vector<double>& dist1,
               const std::vector<double>& dist2);

double L2_norm(const std::vector<double>& x_values,
               const std::vector<double>& dist1,
               const std::vector<double>& dist2);

std::array<double, 2> grid_search_2d(std::function<std::array<double, 2>(std::array<double, 2>)> F, 
                                     const std::array<double, 2>& bounds_min, 
                                     const std::array<double, 2>& bounds_max, 
                                     const int n1=50, const int n2=50,
                                     const bool write_search=false,
                                     const std::string& filename="grid_search.csv");
                        
// nelder-mead 2d minimiser
struct NelderMead2DParams {
    const std::function<std::array<double, 2>(const std::array<double, 2>&)>& func;
};

static double nelder_mead_2d_objective(const gsl_vector* x, void* params);

std::array<double, 2> nelder_mead_minimise_2d(
    const std::function<std::array<double, 2>(const std::array<double, 2>&)>& func,
    double x0,
    double x1,
    double step0 = 0.1,
    double step1 = 0.01,
    double tol = 1e-8,
    int max_iter = 10000);

#endif // INCLUDE_MATHS_OPS_HPP_H