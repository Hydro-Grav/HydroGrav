/**
 * @file maths.hpp
 * @brief Mathematical utilities used in constructing fluid profiles and in the sound shell model.
 */
#ifndef INCLUDE_MATHS_HPP_H
#define INCLUDE_MATHS_HPP_H

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

#include "config.hpp"
#include <gsl/gsl_vector.h>
#include "interpolation.h"
#include "ap.h"

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

/**
 * @brief Generate a vector of logarithmically spaced values.
 *
 * Similar to Python's `numpy.logspace`, this function returns `num`
 * points spaced evenly on a log scale between `start` and `stop`.
 *
 * @param log_start Starting exponent value (base e).
 * @param log_end  Ending exponent value (base e).
 * @param num   Number of points to generate.
 * @return std::vector<double> Log spaced values.
 */
std::vector<double> logspace(double log_start, double log_end, std::size_t num=100);

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

/**
 * @brief Convert a floating-point value to string with fixed precision.
 *
 * @param value     Value to convert.
 * @param precision Number of decimal places.
 * @return Formatted string.
 */
std::string to_string_with_precision(double value, int precision = 2);

/**
 * @brief Integrate y(x) using Simpson's rule on a non-uniform grid.
 *
 * @param x Independent variable samples.
 * @param y Function values corresponding to `x`.
 * @return Approximate integral \f$\int y\,dx\f$.
 */
double simpson_integrate(const std::vector<double>& x, const std::vector<double>& y);

/**
 * @brief Minimise a scalar function using the golden-section search method.
 *
 * Locates a local minimum of @p f in the interval [a, b] without
 * requiring derivative information.
 *
 * @param f        Objective function.
 * @param a        Left bound of search interval.
 * @param b        Right bound of search interval.
 * @param tol      Convergence tolerance.
 * @param max_iter Maximum number of iterations.
 * @return Approximate minimiser.
 */
double golden_section_minimize(std::function<double(double)> f, double a, double b, double tol = config::golden_sec_tol, int max_iter = config::golden_sec_max_its);

/**
 * @brief Find a bracketing interval [a', b'] ⊆ [a, b] in which @p residual_func changes sign.
 *
 * Scans the interval adaptively until a sign change is found, then
 * returns the tightest bracket located.
 *
 * @param residual_func Function whose root is being bracketed.
 * @param a             Initial lower bound.
 * @param b             Initial upper bound.
 * @return Array {a', b'} such that residual_func(a') and residual_func(b') have opposite signs.
 */
std::array<double, 2> find_bracket(const std::function<double(double)>& residual_func, double a, double b);

/**
 * @brief Fourth-order Runge–Kutta ODE solver.
 *
 * Integrates the system \f$dy/dx = f(x, y)\f$ from @p x0 to @p xf
 * with @p n uniform steps, starting from initial condition @p y0.
 *
 * @tparam T     Floating-point type for the independent variable.
 * @tparam State Container type for the dependent variable (e.g. std::array).
 * @tparam Func  Callable with signature State(T, const State&).
 *
 * @param dydx  Right-hand-side function.
 * @param x0    Initial value of independent variable.
 * @param xf    Final value of independent variable.
 * @param y0    Initial state.
 * @param n     Number of steps (default: 1000).
 * @return Pair of {x_values, y_values} over the integration range.
 */
template <typename T, typename State, typename Func>
std::pair<std::vector<T>, std::vector<State>> rk4_solver(const Func& dydx, T x0, T xf, const State& y0, size_t n=config::rk4_steps);

/**
 * @brief Bounded two-dimensional Newton–Raphson root finder.
 *
 * Extends @c newton_solve_2d by clamping iterates within [x_min, x_max]
 * component-wise at each step.
 *
 * @tparam T    Floating-point type.
 * @tparam Func Callable with signature std::array<T,2>(std::array<T,2>).
 *
 * @param F            2D residual function.
 * @param x0           Initial guess.
 * @param x_min        Component-wise lower bounds.
 * @param x_max        Component-wise upper bounds.
 * @param tol          Convergence tolerance.
 * @param max_iter     Maximum iterations.
 * @param h            Finite-difference step.
 * @param bound_margin Minimum distance from domain boundaries.
 * @param dev          Enable verbose diagnostic output.
 * @return Approximate root.
 */
template <typename T, typename Func>
std::array<T,2> newton_solve_2d_bounded(
    const Func& F,
    std::array<T,2> x0,
    std::array<T,2> x_min,  // Lower bounds
    std::array<T,2> x_max,  // Upper bounds
    T tol = config::newton_solver_tol,
    int max_iter = config::newton_solver_max_its,
    T h = 1e-8,
    T bound_margin = 1e-6,  // Stay this far from boundaries
    bool dev = false
);

#include "solvers.tpp"

/**
 * @brief Compute the 1D Wasserstein (Earth Mover's) distance between two weighted distributions.
 *
 * @param u_values  Support points of the first distribution.
 * @param u_weights Weights of the first distribution (need not sum to 1).
 * @param v_values  Support points of the second distribution.
 * @param v_weights Weights of the second distribution.
 * @return Wasserstein-1 distance between the two distributions.
 */
double wasserstein_distance_1d(std::vector<double> u_values, 
                                std::vector<double> u_weights,
                                std::vector<double> v_values, 
                                std::vector<double> v_weights);

/**
 * @brief Compute the L1 norm between two distributions sampled on a common grid.
 *
 * @param x_values Common sample points.
 * @param dist1    Values of the first distribution at @p x_values.
 * @param dist2    Values of the second distribution at @p x_values.
 * @return \f$\int |f_1(x) - f_2(x)|\,dx\f$ approximated via trapezoidal rule.
 */
double L1_norm(const std::vector<double>& x_values,
               const std::vector<double>& dist1,
               const std::vector<double>& dist2);

/**
 * @brief Compute the L2 norm between two distributions sampled on a common grid.
 *
 * @param x_values Common sample points.
 * @param dist1    Values of the first distribution at @p x_values.
 * @param dist2    Values of the second distribution at @p x_values.
 * @return \f$\left(\int |f_1(x) - f_2(x)|^2\,dx\right)^{1/2}\f$ approximated via trapezoidal rule.
 */
double L2_norm(const std::vector<double>& x_values,
               const std::vector<double>& dist1,
               const std::vector<double>& dist2);

/**
 * @brief Find the approximate zero of a 2D residual function by exhaustive grid search.
 *
 * Evaluates @p F on a regular n1 × n2 grid over the specified bounds
 * and returns the grid point with the smallest combined residual magnitude.
 * Useful as an initial guess for Newton or Nelder–Mead refinement.
 *
 * @param F             2D residual function.
 * @param bounds_min    Lower corner of search domain.
 * @param bounds_max    Upper corner of search domain.
 * @param n1            Grid resolution along first axis (default: 50).
 * @param n2            Grid resolution along second axis (default: 50).
 * @param write_search  Write the full search grid to file if @c true.
 * @param filename      Output filename for the search grid (default: "grid_search.csv").
 * @return Best-fit point as a 2-element array.
 */
std::array<double, 2> grid_search_2d(std::function<std::array<double, 2>(std::array<double, 2>)> F, 
                                     const std::array<double, 2>& bounds_min, 
                                     const std::array<double, 2>& bounds_max, 
                                     const int n1=50, const int n2=50,
                                     const bool write_search=false,
                                     const std::string& filename="grid_search.csv");

/**
 * @struct NelderMead2DParams
 * @brief Parameter bundle passed to the GSL Nelder–Mead objective function.
 *
 * Wraps the user-supplied 2D residual function so it can be forwarded
 * through the GSL C interface via a @c void* pointer.
 */
struct NelderMead2DParams {
    const std::function<std::array<double, 2>(const std::array<double, 2>&)>& func;
};

/**
 * @brief GSL-compatible objective function for 2D Nelder–Mead minimisation.
 *
 * Evaluates the squared L2 norm of the residual vector returned by the
 * user function stored in @p params.  This function conforms to the
 * @c gsl_multimin_function signature and is used internally by
 * @c nelder_mead_minimise_2d.
 *
 * @param x      Current iterate as a GSL vector.
 * @param params Pointer to a @c NelderMead2DParams instance.
 * @return Scalar objective value \f$\sqrt{r_0^2 + r_1^2}\f$.
 */
double nelder_mead_2d_objective(const gsl_vector* x, void* params);

/**
 * @brief Minimise a 2D residual function using the Nelder–Mead simplex algorithm.
 *
 * Wraps the GSL @c nmsimplex2 minimiser.  The objective is the squared
 * L2 norm of the 2-component residual returned by @p func.
 *
 * @param func     2D residual function to minimise.
 * @param x0       Initial guess for the first component.
 * @param x1       Initial guess for the second component.
 * @param step0    Initial simplex step size for the first component (default: 0.1).
 * @param step1    Initial simplex step size for the second component (default: 0.01).
 * @param tol      Simplex size convergence tolerance (default: 1e-8).
 * @param max_iter Maximum number of iterations (default: 10000).
 * @return Minimiser as a 2-element array.
 * @throws std::runtime_error if the minimiser does not converge within @p max_iter.
 */
std::array<double, 2> nelder_mead_minimise_2d(
    const std::function<std::array<double, 2>(const std::array<double, 2>&)>& func,
    double x0,
    double x1,
    double step0 = 0.1,
    double step1 = 0.01,
    double tol = config::nelder_mead_tol,
    int max_iter = config::nelder_mead_max_its);

std::pair<double, double> get_mean_sd(std::vector<double>& vec);

#endif // INCLUDE_MATHS_HPP_H