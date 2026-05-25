// maths_ops.cpp
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <functional>
#include <cassert>
#include <numeric>
#include <memory>
// #include <matplotlibcpp.h>

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

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

std::vector<double> logspace(double log_start, double log_end, std::size_t num) {
    std::vector<double> result;
    // double log_start = std::log10(start);  // Log of the start value
    // double log_stop = std::log10(stop);    // Log of the stop value

    double step = (log_end - log_start) / (num - 1);  // Step size in log space

    for (size_t i = 0; i < num; ++i) {
        double log_val = log_start + i * step;     // Calculate log value at step i
        result.push_back(std::pow(10, log_val));   // Convert back to linear space
    }

    return result;
}

/*
    The Filon qadrature method below was created using Claude Sonnet 4.5.
    We do not claim ownership of this code, but we have verified it works
    as intended. We caution that unexpected results may occur.

    TODO - rename to FilonIntegrator
*/

FilonQuadrature::FilonQuadrature(int n_points) : n_points_(n_points) {
    if (n_points < 3) {
        throw std::invalid_argument("FilonQuadrature requires at least 3 points");
    }
}

double FilonQuadrature::filon_integrate_interval(
    const std::vector<double>& f_vals,
    const std::vector<double>& x_vals,
    double omega,
    double a,
    double b,
    bool use_sin)
{
    // Filon-type quadrature: approximate f(x) with piecewise linear/quadratic
    // and integrate exactly for polynomial * sin(omega*x) or polynomial * cos(omega*x)
    
    const int n = f_vals.size();
    if (n < 2) return 0.0;
    
    double result = 0.0;
    
    // Use Simpson-type approach: trapezoidal for 2 points, Simpson for 3+
    if (n == 2) {
        // Trapezoidal with exact oscillatory integration
        // ∫[a,b] (f0 + (f1-f0)*(x-a)/(b-a)) * sin(omega*x) dx
        double h = b - a;
        double f0 = f_vals[0];
        double f1 = f_vals[1];
        
        if (use_sin) {
            // ∫[a,b] f0*sin(ωx) dx = -f0/ω * [cos(ωb) - cos(ωa)]
            double s0 = -f0 / omega * (std::cos(omega * b) - std::cos(omega * a));
            // ∫[a,b] (x-a)*sin(ωx) dx
            double s1 = ((f1 - f0) / h) * 
                        ((std::sin(omega * b) - std::sin(omega * a)) / omega - 
                         a * (std::cos(omega * b) - std::cos(omega * a)) / omega);
            result = s0 + s1;
        } else {
            // ∫[a,b] f0*cos(ωx) dx = f0/ω * [sin(ωb) - sin(ωa)]
            double c0 = f0 / omega * (std::sin(omega * b) - std::sin(omega * a));
            double c1 = ((f1 - f0) / h) * 
                        (-(std::cos(omega * b) - std::cos(omega * a)) / omega - 
                         a * (std::sin(omega * b) - std::sin(omega * a)) / omega);
            result = c0 + c1;
        }
    } else {
        // Use composite integration over subintervals
        for (int i = 0; i + 1 < n; ++i) {
            double x0 = x_vals[i];
            double x1 = x_vals[i + 1];
            double f0 = f_vals[i];
            double f1 = f_vals[i + 1];
            double h = x1 - x0;
            
            if (std::abs(h) < 1e-15) continue;
            
            if (use_sin) {
                // Linear interpolation: f(x) = f0 + (f1-f0)*(x-x0)/h
                // ∫[x0,x1] f(x)*sin(ωx) dx
                double wa = omega * x0;
                double wb = omega * x1;
                double ca = std::cos(wa);
                double cb = std::cos(wb);
                double sa = std::sin(wa);
                double sb = std::sin(wb);
                
                result += -f0 * (cb - ca) / omega + 
                         (f1 - f0) / (h * omega * omega) * 
                         (omega * h * (cb + ca) / 2.0 - (sb - sa));
            } else {
                // ∫[x0,x1] f(x)*cos(ωx) dx
                double wa = omega * x0;
                double wb = omega * x1;
                double ca = std::cos(wa);
                double cb = std::cos(wb);
                double sa = std::sin(wa);
                double sb = std::sin(wb);
                
                result += f0 * (sb - sa) / omega + 
                         (f1 - f0) / (h * omega * omega) * 
                         (omega * h * (sb + sa) / 2.0 + (cb - ca));
            }
        }
    }
    
    return result;
}

double FilonQuadrature::integrate_sin(
    const std::function<double(double)>& f,
    double omega,
    double a,
    double b)
{
    if (std::abs(omega) < 10.0) {
        // For small omega, use standard quadrature
        auto integrand = [&](double x) { return f(x) * std::sin(omega * x); };
        boost::math::quadrature::gauss<double, 30> integrator;
        return integrator.integrate(integrand, a, b);
    }
    
    // Adaptive subdivision based on oscillation count
    const double interval = b - a;
    const double period = 2.0 * M_PI / omega;
    const int n_oscillations = static_cast<int>(std::ceil(interval / period));
    
    // Target: 2-4 oscillations per subinterval for stability
    const int target_osc_per_interval = 3;
    const int n_subdivisions = std::max(1, (n_oscillations + target_osc_per_interval - 1) / target_osc_per_interval);
    
    const double h_sub = interval / n_subdivisions;
    double total = 0.0;
    
    for (int k = 0; k < n_subdivisions; ++k) {
        double a_sub = a + k * h_sub;
        double b_sub = a + (k + 1) * h_sub;
        
        // Sample function at n_points_ locations in this subinterval
        std::vector<double> x_vals(n_points_);
        std::vector<double> f_vals(n_points_);
        
        for (int i = 0; i < n_points_; ++i) {
            double t = static_cast<double>(i) / (n_points_ - 1);  // uniform in [0,1]
            x_vals[i] = a_sub + t * (b_sub - a_sub);
            f_vals[i] = f(x_vals[i]);
        }
        
        // Apply Filon-type quadrature on this subinterval
        total += filon_integrate_interval(f_vals, x_vals, omega, a_sub, b_sub, true);
    }
    
    return total;
}

double FilonQuadrature::integrate_cos(
    const std::function<double(double)>& f,
    double omega,
    double a,
    double b)
{
    if (std::abs(omega) < 10.0) {
        // For small omega, use standard quadrature
        auto integrand = [&](double x) { return f(x) * std::cos(omega * x); };
        boost::math::quadrature::gauss<double, 30> integrator;
        return integrator.integrate(integrand, a, b);
    }
    
    // Adaptive subdivision based on oscillation count
    const double interval = b - a;
    const double period = 2.0 * M_PI / omega;
    const int n_oscillations = static_cast<int>(std::ceil(interval / period));
    
    // Target: 2-4 oscillations per subinterval for stability
    const int target_osc_per_interval = 3;
    const int n_subdivisions = std::max(1, (n_oscillations + target_osc_per_interval - 1) / target_osc_per_interval);
    
    const double h_sub = interval / n_subdivisions;
    double total = 0.0;
    
    for (int k = 0; k < n_subdivisions; ++k) {
        double a_sub = a + k * h_sub;
        double b_sub = a + (k + 1) * h_sub;
        
        // Sample function at n_points_ locations in this subinterval
        std::vector<double> x_vals(n_points_);
        std::vector<double> f_vals(n_points_);
        
        for (int i = 0; i < n_points_; ++i) {
            double t = static_cast<double>(i) / (n_points_ - 1);  // uniform in [0,1]
            x_vals[i] = a_sub + t * (b_sub - a_sub);
            f_vals[i] = f(x_vals[i]);
        }
        
        // Apply Filon-type quadrature on this subinterval
        total += filon_integrate_interval(f_vals, x_vals, omega, a_sub, b_sub, false);
    }
    
    return total;
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

    // --- 4. Check if a turning point (local minimum) exists --------------
    // A turning point exists if the minimum is strictly interior to the segment
    // and the function decreases then increases (has lower neighbors on both sides)
    bool has_turning_point = false;
    
    if (idx_min > S.i0 && idx_min < S.i1) {
        // Minimum is interior to segment
        // Check if it's a local minimum (lower than both neighbors)
        if (ys[idx_min] < ys[idx_min - 1] && ys[idx_min] < ys[idx_min + 1]) {
            has_turning_point = true;
        }
    }
    
    // Alternative check: scan for any interior local minimum in the segment
    if (!has_turning_point) {
        // for (int j = S.i0 + 1; j < S.i1; ++j) {
        for (int j = S.i0; j < S.i1; ++j) {
            if (ys[j] < ys[j - 1] && ys[j] < ys[j + 1]) {
                has_turning_point = true;
                break;
            }
        }
    }
    
    if (!has_turning_point) {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
    }

    // --- 5. Construct bracket from that segment --------------------------
    double dx  = (b - a) / NSAMPLES;
    double left  = std::max(a, xs[S.i0] - PAD_MULT * dx);
    double right = std::min(b, xs[S.i1] + PAD_MULT * dx);

    return {left, right};
}

double wasserstein_distance_1d(std::vector<double> u_values, 
                                std::vector<double> u_weights,
                                std::vector<double> v_values, 
                                std::vector<double> v_weights) {
    
    // Check inputs
    if (u_values.size() != u_weights.size() || v_values.size() != v_weights.size()) {
        throw std::invalid_argument("Values and weights must be the same size!");
    }

    // Check for sufficient points
    if (u_values.size() < 3 || v_values.size() < 3) {
        throw std::invalid_argument("At least 3 points required in each distribution for Wasserstein distance calculation.");
    }
    
    // Normalize weights
    double u_sum = 0.0, v_sum = 0.0;
    for (double w : u_weights) {
        if (w < 0.0) throw std::invalid_argument("Weights must be non-negative!");
        u_sum += w;
    }
    for (double w : v_weights) {
        if (w < 0.0) throw std::invalid_argument("Weights must be non-negative!");
        v_sum += w;
    }
    for (double& w : u_weights) w /= u_sum;
    for (double& w : v_weights) w /= v_sum;
    
    // Create (value, weight) pairs and sort
    std::vector<std::pair<double, double>> u_dist, v_dist;
    for (size_t i = 0; i < u_values.size(); i++) {
        u_dist.push_back({u_values[i], u_weights[i]});
    }
    for (size_t i = 0; i < v_values.size(); i++) {
        v_dist.push_back({v_values[i], v_weights[i]});
    }
    
    std::sort(u_dist.begin(), u_dist.end());
    std::sort(v_dist.begin(), v_dist.end());
    
    // Merge the two sorted arrays and compute distance
    size_t i = 0, j = 0;
    double u_cdf = 0.0, v_cdf = 0.0;
    double distance = 0.0;
    double prev_x = std::min(u_dist[0].first, v_dist[0].first);
    
    while (i < u_dist.size() || j < v_dist.size()) {
        // Find next position
        double u_x = (i < u_dist.size()) ? u_dist[i].first : INFINITY;
        double v_x = (j < v_dist.size()) ? v_dist[j].first : INFINITY;
        double curr_x = std::min(u_x, v_x);
        
        // Add area: |CDF_u - CDF_v| * width
        distance += std::abs(u_cdf - v_cdf) * (curr_x - prev_x);
        
        // Update CDFs at this position
        if (std::abs(curr_x - u_x) < 1e-15 && i < u_dist.size()) {
            u_cdf += u_dist[i].second;
            i++;
        }
        if (std::abs(curr_x - v_x) < 1e-15 && j < v_dist.size()) {
            v_cdf += v_dist[j].second;
            j++;
        }
        
        prev_x = curr_x;
    }
    
    return distance;
}

double L2_norm(const std::vector<double>& x_values,
               const std::vector<double>& dist1,
               const std::vector<double>& dist2) {
    
    if (x_values.size() != dist1.size() || x_values.size() != dist2.size()) {
        throw std::invalid_argument(
            "L2_norm: x_values, dist1, and dist2 must have the same size. "
            "Got sizes: " + std::to_string(x_values.size()) + ", " +
            std::to_string(dist1.size()) + ", " + std::to_string(dist2.size())
        );
    }
    
    if (x_values.empty()) {
        throw std::invalid_argument("L2_norm: Input vectors cannot be empty");
    }
    
    if (x_values.size() < 3) {
        throw std::invalid_argument("L2_norm: Need at least 3 points for Simpson integration");
    }

    // Normalise distributions
    std::vector<double> d1 = dist1;
    std::vector<double> d2 = dist2;

    double norm1 = simpson_integrate(x_values, d1);
    double norm2 = simpson_integrate(x_values, d2);
    
    if (norm1 <= 0.0 || norm2 <= 0.0) {
        throw std::runtime_error("L1_norm: Cannot normalize - integral is non-positive");
    }
    
    for (size_t i = 0; i < d1.size(); ++i) {
        d1[i] /= norm1;
        d2[i] /= norm2;
    }
    
    // Calculate (P1 - P2)²
    std::vector<double> squared_diff(x_values.size());
    for (size_t i = 0; i < x_values.size(); ++i) {
        double diff = d1[i] - d2[i];
        squared_diff[i] = diff * diff;
    }
    
    // Integrate
    double integral = simpson_integrate(x_values, squared_diff);
    
    // Return sqrt of integral
    return std::sqrt(std::max(0.0, integral));  // max protects against numerical noise giving negative
}

/**
 * Calculate L1 norm (Manhattan distance) between two distributions. For log(x) values, normalises distributions such that
 * ∫ P(log(x)) dlog(x) = 1
 * 
 * L1 = ∫ |P1(x) - P2(x)| dx
 * 
 * @param x_values Positions/coordinates
 * @param dist1 First distribution values
 * @param dist2 Second distribution values
 * @return L1 norm between the two distributions
 */
double L1_norm(const std::vector<double>& x_values,
               const std::vector<double>& dist1,
               const std::vector<double>& dist2) {
    
    if (x_values.size() != dist1.size() || x_values.size() != dist2.size()) {
        throw std::invalid_argument(
            "L1_norm: x_values, dist1, and dist2 must have the same size"
        );
    }
    
    if (x_values.empty()) {
        throw std::invalid_argument("L1_norm: Input vectors cannot be empty");
    }
    
    if (x_values.size() < 3) {
        throw std::invalid_argument("L1_norm: Need at least 3 points for Simpson integration");
    }
    
    // Normalise distributions
    std::vector<double> d1 = dist1;
    std::vector<double> d2 = dist2;

    double norm1 = simpson_integrate(x_values, d1);
    double norm2 = simpson_integrate(x_values, d2);
    
    if (norm1 <= 0.0 || norm2 <= 0.0) {
        throw std::runtime_error("L1_norm: Cannot normalize - integral is non-positive");
    }
    
    for (size_t i = 0; i < d1.size(); ++i) {
        d1[i] /= norm1;
        d2[i] /= norm2;
    }
    
    // Calculate |P1 - P2|
    std::vector<double> abs_diff(x_values.size());
    for (size_t i = 0; i < x_values.size(); ++i) {
        abs_diff[i] = std::abs(d1[i] - d2[i]);
    }
    
    // Integrate
    return simpson_integrate(x_values, abs_diff);
}

std::array<double, 2> grid_search_2d(std::function<std::array<double, 2>(std::array<double, 2>)> F, 
                                           const std::array<double, 2>& bounds_min, 
                                           const std::array<double, 2>& bounds_max, 
                                           const int n1, const int n2, 
                                           const bool write_search, const std::string& filename)
{
    double best_residual_norm = std::numeric_limits<double>::infinity();
    std::array<double, 2> best_x = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};

    const auto x1_min = bounds_min[0];
    const auto x1_max = bounds_max[0];
    const auto x2_min = bounds_min[1];
    const auto x2_max = bounds_max[1];

    std::unique_ptr<std::ofstream> out;
    if (write_search) {
        out = std::make_unique<std::ofstream>(filename);
        *out << std::setprecision(12);
        *out << "x1,x2,eq1,eq2,norm\n";
    }
    
    for (int i = 0; i < n1; ++i) {
        double x1 = x1_min + i * (x1_max - x1_min) / n1;
        for (int j = 0; j < n2; ++j) {
            double x2 = x2_min + j * (x2_max - x2_min) / n2;
            
            try {
                auto res = F({x1, x2});
                double res_norm = std::sqrt(res[0]*res[0] + res[1]*res[1]);

                if (out) {
                    *out << x1 << "," << x2 << "," << res[0] << "," << res[1] 
                         << "," << res_norm << "\n";
                }
                
                if (res_norm < best_residual_norm) {
                    best_residual_norm = res_norm;
                    best_x = {x1, x2};
                }
            } catch (...) {
                continue;
            }
        }
    }

    out.reset();

    return best_x;
}

// nelder-mead 2d minimiser
double nelder_mead_2d_objective(const gsl_vector* x, void* params) {
    auto* p = static_cast<NelderMead2DParams*>(params);

    const std::array<double, 2> args = {
        gsl_vector_get(x, 0),
        gsl_vector_get(x, 1)
    };

    try {
        const auto r = p->func(args);
        return std::sqrt(r[0] + r[1]); // squared residual
    } catch (const std::exception&) {
        return 1e10;
    }
}

std::array<double, 2> nelder_mead_minimise_2d(
    const std::function<std::array<double, 2>(const std::array<double, 2>&)>& func,
    double x0,
    double x1,
    double step0,
    double step1,
    double tol,
    int max_iter)
{
    NelderMead2DParams params{func};

    gsl_multimin_function gsl_func;
    gsl_func.n      = 2;
    gsl_func.f      = &nelder_mead_2d_objective;
    gsl_func.params = &params;

    gsl_vector* x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, x0);
    gsl_vector_set(x, 1, x1);

    gsl_vector* step_size = gsl_vector_alloc(2);
    gsl_vector_set(step_size, 0, step0);
    gsl_vector_set(step_size, 1, step1);

    gsl_multimin_fminimizer* minimizer = 
        gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, 2);
    gsl_multimin_fminimizer_set(minimizer, &gsl_func, x, step_size);

    int status = GSL_CONTINUE;
    int iter   = 0;

    while (status == GSL_CONTINUE && iter < max_iter) {
        ++iter;
        status = gsl_multimin_fminimizer_iterate(minimizer);
        if (status) break;
        const double simplex_size = gsl_multimin_fminimizer_size(minimizer);
        status = gsl_multimin_test_size(simplex_size, tol);
    }

    const double x0_sol = gsl_vector_get(minimizer->x, 0);
    const double x1_sol = gsl_vector_get(minimizer->x, 1);
    const double fval   = minimizer->fval;

    gsl_multimin_fminimizer_free(minimizer);
    gsl_vector_free(x);
    gsl_vector_free(step_size);

    if (status != GSL_SUCCESS) {
        throw std::runtime_error(
            "Nelder-Mead minimisation failed to converge after " +
            std::to_string(iter) + " iterations. Final residual: " +
            std::to_string(fval));
    }

    return {x0_sol, x1_sol};
}