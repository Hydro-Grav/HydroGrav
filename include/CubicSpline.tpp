#ifndef INCLUDE_CUBICSPLINE_TPP_H
#define INCLUDE_CUBICSPLINE_TPP_H

#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cmath>

template <typename T>
CubicSpline<T>::CubicSpline() : initialised_(false) {}

template <typename T>
CubicSpline<T>::CubicSpline(const std::vector<T>& x, const std::vector<T>& y) {
    build(x, y);
}

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

// scalar operations
template <typename T>
CubicSpline<T> CubicSpline<T>::operator+(T scalar) const { // spline + scalar
    std::vector<T> y_new = y_;
    for (T& val : y_new) val += scalar;
    return CubicSpline<T>(x_, y_new);
}

template <typename T>
CubicSpline<T> CubicSpline<T>::operator-(T scalar) const { // spline - scalar
    return *this + (-scalar);
}

template <typename T>
CubicSpline<T> CubicSpline<T>::operator*(T scalar) const { // spline * scalar
    std::vector<T> y_new = y_;
    for (T& val : y_new) val *= scalar;
    return CubicSpline<T>(x_, y_new);
}

template <typename T>
CubicSpline<T> CubicSpline<T>::operator/(T scalar) const { // spline / scalar
    if (scalar == T(0)) throw std::invalid_argument("Division by zero.");
    return *this * (T(1) / scalar);
}

// scalar operations (from left)
template <typename T>
CubicSpline<T> operator+(T scalar, const CubicSpline<T>& spline) { // scalar + spline
    return spline + scalar;
}

template <typename T>
CubicSpline<T> operator-(T scalar, const CubicSpline<T>& spline) { // scalar - spline
    return spline * T(-1) + scalar;
}

template <typename T>
CubicSpline<T> operator*(T scalar, const CubicSpline<T>& spline) { // scalar * spline
    return spline * scalar;
}

#endif // INCLUDE_CUBICSPLINE_TPP_H
