// tests.hpp
#ifndef INCLUDE_TESTS_HPP_H
#define INCLUDE_TESTS_HPP_H

/**
 * @file tests.hpp
 * @brief Contains test functions and plotting utilities for various components of the codebase.
*/


#include <vector>
#include <string>
// #include <matplotlibcpp.h>

/// @brief Runs all tests
void test_all();

void test_PowerSpec();

void example_Kin_Spec(const std::string& filename="example_kin_spectrum");
void example_GW_Spec(const std::string& filename="example_GW_spectrum");

void test_cubic_spline();

void test_SiCi_spline();


void test_FluidProfile_params();

void test_rk4_coupled_odes();
void test_rk4_solver();

/// @brief Tests arithmetic and operations on the vec<T> class.
void test_vec();

/// @brief Tests the FluidProfile class for correct behavior and integration.
void test_FluidProfile();

/**
 
 * 
 * @param x Input value.
 * @param params Pointer to parameters needed for the function.
 * 
 * @return Computed value.
*/


/// @brief Tests numerical integration using Boost.Math.
void test_boostmath();



/**
 * @brief Tests the 1D interpolator with a specified type.
 * 
 * @param type Type of interpolator to test (e.g., "linear", "cubic").
*/
void test_interpolator();

void test_prof_ints();
void test_Apsq();

void test_SiCi();

// void test_dlt_SSM();

void test_simpson_integrate();

void test_integrators();

/// @brief Namespace for performance testing alternatives to std::pow
namespace PowerTest {

/**
 * @brief Computes the mean of a list of time durations.
 * 
 * @param times Vector of measured time values.
 * @return Mean of the time values.
 */
double mean(const std::vector<double>& times);

/**
 * @brief Computes the standard deviation of time measurements.
 * 
 * @param times Vector of measured time values.
 * @param mean_time Mean of the time values.
 * @return Standard deviation of the time values.
 */
double st_dev(const std::vector<double>& times, double mean_time);

/**
 * @brief Run a performance test on PowerSpec and report timing statistics.
 * 
 * @param num_runs Number of iterations for benchmarking (default 1000).
 */
void test_power(int num_runs = 1000);

} // namespace PowerTest

/**
 * @brief Plots a user-defined function over a given interval.
 * 
 * @tparam Func Type of the callable object.
 * @param f Function to plot.
 * @param xmin Minimum x-value of the plot.
 * @param xmax Maximum x-value of the plot.
 * @param points Number of sample points (default: 100).
 */
template<typename Func> // only takes in function of 1 variable at this stage
// void plot_func(Func f, double xmin, double xmax, int points=100) {
//     namespace plt = matplotlibcpp;
    
//     std::vector<double> xs, ys;
//     for (int i = 0; i <= points; ++i) {
//         double x = xmin + i * (xmax - xmin) / points;
//         xs.push_back(x);
//         ys.push_back(f(x));
//     }
    
//     plt::plot(xs, ys);
//     plt::xlabel("x");
//     plt::ylabel("f(x)");
//     plt::grid(true);
//     // plt::show();
//     plt::save("plot.png");
// }

// move this to a more appropriate file
/**
 * @brief Reads a file containing point data and generates a plot.
 * 
 * @param filename Path to the data file (CSV or plain text).
 */
void plot_pts(std::string& filename);

#endif // INCLUDE_TESTS_HPP_H
