/**
 * @file CubicSpline.tpp
 * @brief Template implementation of the CubicSpline class.
 *
 * Implements natural cubic spline interpolation using the standard
 * tridiagonal system formulation.  The spline passes exactly through
 * all supplied data points and satisfies the natural boundary condition
 * \f$S''(x_0) = S''(x_{n-1}) = 0\f$.
 *
 * This file is included directly by `maths_ops.hpp` and should not be
 * compiled separately.
 */
#ifndef INCLUDE_CUBICSPLINE_TPP_H
#define INCLUDE_CUBICSPLINE_TPP_H

#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cmath>

/**
 * @brief Default constructor. Leaves the spline uninitialised.
 *
 * @c build() must be called before the spline can be evaluated.
 */
template <typename T>
CubicSpline<T>::CubicSpline() : initialised_(false) {}

/**
 * @brief Construct and immediately build the spline from data.
 *
 * Equivalent to default-constructing and then calling @c build(x, y).
 *
 * @param x Strictly monotonic sample points.
 * @param y Function values at each sample point.
 */
template <typename T>
CubicSpline<T>::CubicSpline(const std::vector<T>& x, const std::vector<T>& y) {
    build(x, y);
}

/**
 * @brief Build the spline coefficients from data vectors.
 *
 * Computes the piecewise cubic polynomial coefficients \f$a_i, b_i, c_i, d_i\f$
 * by solving the tridiagonal system arising from the natural cubic spline
 * conditions.  If @p x is strictly decreasing the arrays are reversed
 * internally so the underlying algorithm always operates on an increasing
 * sequence.
 *
 * @param x Strictly monotonic (increasing or decreasing) sample points.
 *          Must contain at least two elements.
 * @param y Function values corresponding to @p x.  Must be the same size as @p x.
 *
 * @throws std::invalid_argument If sizes differ, fewer than two points are
 *         supplied, any x-value is NaN, or the sequence is not strictly monotonic.
 */
template <typename T>
void CubicSpline<T>::build(const std::vector<T>& x, const std::vector<T>& y) {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::invalid_argument("CubicSpline: x and y must be the same size and contain at least two points.");
    }

    const size_t n = x.size();

    // might not be the best way to implement strictly increasing - do in solve_profile() instead?
    // Checks monotonicity
    if (!is_strictly_monotonic(x)) {
        throw std::invalid_argument("CubicSpline: x-values must be strictly increasing or decreasing");
    }

    // Cubic spline needs x strictly increasing
    bool is_increasing = x[0] < x[1];
    auto x_copy = x;
    auto y_copy = y;

    if (!is_increasing) {
        std::reverse(x_copy.begin(), x_copy.end());
        std::reverse(y_copy.begin(), y_copy.end());
    }

    x_ = x_copy;
    y_ = y_copy;
    h_.resize(n - 1);
    a_ = y_copy;
    c_.resize(n);

    // Step 1: Compute h[i]
    for (size_t i = 0; i < n - 1; ++i) {
        h_[i] = x_[i + 1] - x_[i];
    }

    // Step 2: Compute alpha
    std::vector<T> alpha(n - 1, 0);
    for (size_t i = 1; i < n - 1; ++i) {
        alpha[i] = (3.0 / h_[i]) * (a_[i + 1] - a_[i]) - (3.0 / h_[i - 1]) * (a_[i] - a_[i - 1]);
    }

    // Step 3: Solve tridiagonal system
    std::vector<T> l(n), mu(n), z(n);
    l[0] = 1;
    mu[0] = z[0] = 0;

    for (size_t i = 1; i < n - 1; ++i) {
        l[i] = 2 * (x_[i + 1] - x_[i - 1]) - h_[i - 1] * mu[i - 1];
        mu[i] = h_[i] / l[i];
        z[i] = (alpha[i] - h_[i - 1] * z[i - 1]) / l[i];
    }

    l[n - 1] = 1;
    z[n - 1] = c_[n - 1] = 0;

    b_.resize(n - 1);
    d_.resize(n - 1);

    for (int j = static_cast<int>(n) - 2; j >= 0; --j) {
        c_[j] = z[j] - mu[j] * c_[j + 1];
        b_[j] = (a_[j + 1] - a_[j]) / h_[j] - h_[j] * (c_[j + 1] + 2 * c_[j]) / 3;
        d_[j] = (c_[j + 1] - c_[j]) / (3 * h_[j]);
    }

    initialised_ = true;

    check_convergence();
    // std::cout << "CubicSpline has been initialised and built." << std::endl;

    return;
}

/**
 * @brief Check whether a vector is strictly monotonic (increasing or decreasing).
 *
 * A sequence of fewer than two elements is considered trivially monotonic.
 * Any NaN value triggers an exception.
 *
 * @param x Sequence to test.
 * @return @c true if @p x is strictly increasing or strictly decreasing.
 * @throws std::invalid_argument If a NaN value is detected.
 */
template <typename T>
bool CubicSpline<T>::is_strictly_monotonic(const std::vector<T>& x) const {
    if (x.size() < 2) return true; // Trivially monotonic

    bool increasing = true;
    bool decreasing = true;

    for (size_t i = 1; i < x.size(); ++i) {
        if (std::isnan(x[i]))
            throw std::invalid_argument("CubicSpline: NaN x-value detected.");
        if (x[i] <= x[i - 1]) increasing = false;
        if (x[i] >= x[i - 1]) decreasing = false;
    }

    return increasing || decreasing;
}

/**
 * @brief Verify that the spline reproduces the input data within tolerance.
 *
 * Evaluates the spline at every knot and compares the result to the
 * stored y-values.  A warning is printed to @c stderr for any point
 * whose absolute error exceeds @c tol_.  This is called automatically
 * at the end of @c build().
 */
template <typename T>
void CubicSpline<T>::check_convergence() const {
    bool test_passed = true;

    for (std::size_t i = 0; i < x_.size(); ++i) {
        T interp = (*this)(x_[i]);
        T error = std::abs(interp - y_[i]);

        if (error > tol_) {
            test_passed = false;
            std::cerr << "CubicSpline convergence warning at x = " << x_[i]
                      << ": interpolated = " << interp
                      << ", expected = " << y_[i]
                      << ", error = " << error << '\n';
            // exit on failed convergence?
            // initialised_ = false;
        }
    }

    if (test_passed)
        // std::cout << "CubicSpline convergence test passed (tol=" << tol_ << ")!" <<std::endl;
    
    return;
}

/**
 * @brief Evaluate the spline at a point @p xi.
 *
 * Uses binary search to locate the interval containing @p xi, then
 * evaluates the cubic polynomial for that interval.
 *
 * @param xi Query point.  Must lie within [x.front(), x.back()].
 * @return Interpolated value \f$S(\xi)\f$.
 * @throws std::runtime_error   If the spline has not been initialised.
 * @throws std::invalid_argument If @p xi is outside the data domain.
 */
template <typename T>
T CubicSpline<T>::operator()(T xi) const {
    if (!initialised_) {
        throw std::runtime_error("CubicSpline has not been initialised.");
    }

    // Extrapolate flat for out-of-bounds input
    // if (xi <= x_.front()) return y_.front();
    // if (xi >= x_.back()) return y_.back();

    // Gives error for out-of-bounds input
    if (xi < x_.front() || xi > x_.back()) {
        std::cout << "CubicSpline failed at xi=" << xi << "\n";
        throw std::invalid_argument("CubicSpline called outside of domain bounds!");
    }

    // Binary search to find correct interval
    auto it = std::upper_bound(x_.begin(), x_.end(), xi);
    size_t i = std::distance(x_.begin(), it) - 1;

    T dx = xi - x_[i];
    return a_[i] + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
}

/**
 * @brief Return a new spline whose y-values are each increased by @p scalar.
 * @param scalar Additive offset.
 */
template <typename T>
CubicSpline<T> CubicSpline<T>::operator+(T scalar) const { // spline + scalar
    std::vector<T> y_new = y_;
    for (T& val : y_new) val += scalar;
    return CubicSpline<T>(x_, y_new);
}

/**
 * @brief Return a new spline whose y-values are each decreased by @p scalar.
 * @param scalar Subtractive offset.
 */
template <typename T>
CubicSpline<T> CubicSpline<T>::operator-(T scalar) const { // spline - scalar
    return *this + (-scalar);
}

/**
 * @brief Return a new spline whose y-values are each multiplied by @p scalar.
 * @param scalar Multiplicative factor.
 */
template <typename T>
CubicSpline<T> CubicSpline<T>::operator*(T scalar) const { // spline * scalar
    std::vector<T> y_new = y_;
    for (T& val : y_new) val *= scalar;
    return CubicSpline<T>(x_, y_new);
}

/**
 * @brief Return a new spline whose y-values are each divided by @p scalar.
 * @param scalar Divisor. Must be non-zero.
 * @throws std::invalid_argument If @p scalar is zero.
 */
template <typename T>
CubicSpline<T> CubicSpline<T>::operator/(T scalar) const { // spline / scalar
    if (scalar == T(0)) throw std::invalid_argument("Division by zero.");
    return *this * (T(1) / scalar);
}

/**
 * @brief Return a new spline with @p scalar added to each y-value (scalar + spline).
 * @param scalar Additive offset.
 */
// scalar operations (from left)
template <typename T>
CubicSpline<T> operator+(T scalar, const CubicSpline<T>& spline) { // scalar + spline
    return spline + scalar;
}

/**
 * @brief Return a new spline equal to @p scalar minus each y-value (scalar - spline).
 * @param scalar Value to subtract from.
 */
template <typename T>
CubicSpline<T> operator-(T scalar, const CubicSpline<T>& spline) { // scalar - spline
    return spline * T(-1) + scalar;
}

/**
 * @brief Return a new spline with each y-value multiplied by @p scalar (scalar * spline).
 * @param scalar Multiplicative factor.
 */
template <typename T>
CubicSpline<T> operator*(T scalar, const CubicSpline<T>& spline) { // scalar * spline
    return spline * scalar;
}

#endif // INCLUDE_CUBICSPLINE_TPP_H