/**
 * @file ssm_kernel.hpp
 * @brief Fast evaluation of the sound-shell model time-correlation kernel.
 *
 * The double-time integral that enters @c Spectrum::dlt_SSM reduces, for each of the four
 * sign combinations, to
 *
 * \f[
 *   G(x) \;=\; \tfrac14\Big|\int_{\tau_s}^{\tau_{\rm fin}} e^{i x\tau}\,\frac{d\tau}{\tau}\Big|^2
 *          \;=\; \tfrac14\big[(\mathrm{Ci}(x\tau_{\rm fin})-\mathrm{Ci}(x\tau_s))^2
 *                            +(\mathrm{Si}(x\tau_{\rm fin})-\mathrm{Si}(x\tau_s))^2\big].
 * \f]
 *
 * Evaluating that literally costs two sine/cosine-integral pairs per sign combination
 * (eight per @c dlt_SSM call), and subtracts two nearly-equal special-function values.
 * Rescaling \f$\tau=\tau_s t\f$ removes both problems: the common phase drops out of the
 * modulus, leaving
 *
 * \f[
 *   G(x) = \tfrac14\Big|\int_1^{\rho} e^{i a t}\,\frac{dt}{t}\Big|^2, \qquad
 *   a = x\tau_s,\quad \rho = \tau_{\rm fin}/\tau_s .
 * \f]
 *
 * Mapping \f$t\f$ to \f$[-1,1]\f$ about the interval centre and expanding \f$1/t\f$ in the
 * small parameter \f$\epsilon = \Delta\tau/(\tau_s+\tau_{\rm fin})\f$ gives
 *
 * \f[
 *   G(x) = \tfrac14\,\epsilon^2\,\Big|\sum_{n\ge 0}(-\epsilon)^n M_n(w)\Big|^2, \qquad
 *   M_n(w)=\int_{-1}^{1}\! s^n e^{iws}\,ds, \qquad w = \tfrac12 x\,\Delta\tau .
 * \f]
 *
 * The series is truncated at \f$n=4\f$, so the relative error is \f$O(\epsilon^5)\f$ —
 * verified below @c 1e-6 of the peak for \f$\Delta\tau \le 0.3\,\tau_s\f$.  The leading term
 * alone reproduces the familiar \f$G\simeq(\Delta\tau^2/4\tau_s^2)\,\mathrm{sinc}^2(x\Delta\tau/2)\f$.
 *
 * Two properties are exploited by the caller: \f$G\f$ is even, and it depends on
 * \f$x\f$ only through \f$w=x\Delta\tau/2\f$ — so the four sign combinations needed by
 * @c dlt_SSM share a pair of sine/cosine evaluations via angle addition.
 */

#ifndef INCLUDE_SSM_KERNEL_HPP_H
#define INCLUDE_SSM_KERNEL_HPP_H

#include <cmath>

#include "specialfunctions.h"

namespace Spectrum {

/**
 * @class SoundShellKernel
 * @brief Evaluates \f$G(x)\f$ for a fixed sound-wave period \f$[\tau_s,\tau_{\rm fin}]\f$.
 *
 * Construct once per spectrum; @c operator() is then a few tens of flops.  If the
 * expansion parameter \f$\epsilon\f$ is too large for the truncated series (sound waves
 * lasting an appreciable fraction of a Hubble time) the class silently falls back to a
 * direct ALGLIB Si/Ci evaluation, which is slower but exact.
 */
class SoundShellKernel {
  public:
    /// Largest expansion parameter for which the fast series is used.
    static constexpr double eps_max = 0.15;
    /// Below this |w| the moments are evaluated by their power series rather than in closed
    /// form, avoiding the 1/w^5 cancellation of the closed-form expressions.
    static constexpr double small_w = 1.0;

    SoundShellKernel(double tau_s, double tau_fin)
      : tau_s_(tau_s), tau_fin_(tau_fin),
        half_dtau_(0.5 * (tau_fin - tau_s)),
        eps_((tau_fin - tau_s) / (tau_fin + tau_s)),
        exact_(!(std::abs((tau_fin - tau_s) / (tau_fin + tau_s)) <= eps_max)) {}

    /// Half the sound-wave duration; @c w = x * half_dtau().
    double half_dtau() const { return half_dtau_; }
    /// True when the exact (slow) ALGLIB path is in use.
    bool uses_exact() const { return exact_; }
    /// Expansion parameter @f$\epsilon@f$.
    double eps() const { return eps_; }

    /// @brief \f$G(x)\f$.
    double operator()(double x) const {
        if (exact_) return exact_G(x);
        const double w = x * half_dtau_;
        return from_w(w, std::sin(w), std::cos(w));
    }

    /**
     * @brief \f$G\f$ from a precomputed @p w and its sine/cosine.
     *
     * Lets the caller obtain @c sin/cos for several @p w values by angle addition from a
     * single pair, which is where most of the saving in @c dlt_SSM comes from.
     */
    double from_w(double w, double sin_w, double cos_w) const {
        double m0, m1, m2, m3, m4;
        moments(w, sin_w, cos_w, m0, m1, m2, m3, m4);

        const double e2 = eps_ * eps_;
        const double re = m0 + e2 * (m2 + e2 * m4);
        const double im = eps_ * (m1 + e2 * m3);   // Q = re - i*im

        return 0.25 * e2 * (re * re + im * im);
    }

  private:
    /**
     * @brief Moments \f$M_n(w)=\int_{-1}^1 s^n e^{iws}ds\f$ for n = 0..4.
     *
     * Even @c n are real, odd @c n are purely imaginary; @p m1 and @p m3 are the imaginary
     * parts, i.e. \f$M_1 = i\,m_1\f$.
     */
    static void moments(double w, double sin_w, double cos_w,
                        double& m0, double& m1, double& m2, double& m3, double& m4) {
        if (std::abs(w) <= small_w) {
            // M_n  = 2 sum_m (-1)^m w^2m     / ((2m)!   (n+2m+1))     n even
            // Mt_n = 2 sum_m (-1)^m w^(2m+1) / ((2m+1)! (n+2m+2))     n odd
            const double w2 = w * w;
            double even = 1.0;   // w^(2m) / (2m)!
            double odd  = w;     // w^(2m+1) / (2m+1)!
            double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0;
            for (int m = 0; m < 8; ++m) {
                const double sgn = (m & 1) ? -1.0 : 1.0;
                const double e = sgn * even;
                const double o = sgn * odd;
                s0 += e / (2 * m + 1);
                s1 += o / (2 * m + 3);
                s2 += e / (2 * m + 3);
                s3 += o / (2 * m + 5);
                s4 += e / (2 * m + 5);
                even = odd * w / (2.0 * m + 2.0);          // w^(2m+2)/(2m+2)!
                odd  = even * w / (2.0 * m + 3.0);         // w^(2m+3)/(2m+3)!
            }
            m0 = 2.0 * s0; m1 = 2.0 * s1; m2 = 2.0 * s2; m3 = 2.0 * s3; m4 = 2.0 * s4;
        } else {
            // Closed forms, from M_n = T_n + (i n / w) M_{n-1}.
            const double iw  = 1.0 / w;
            const double iw2 = iw * iw,  iw3 = iw2 * iw;
            const double iw4 = iw2 * iw2, iw5 = iw4 * iw;
            const double S = sin_w, C = cos_w;

            m0 = 2.0 * S * iw;
            m1 = 2.0 * S * iw2 -  2.0 * C * iw;
            m2 = 2.0 * S * iw  -  4.0 * S * iw3 +  4.0 * C * iw2;
            m3 = 6.0 * S * iw2 - 12.0 * S * iw4 -  2.0 * C * iw  + 12.0 * C * iw3;
            m4 = 2.0 * S * iw  - 24.0 * S * iw3 + 48.0 * S * iw5
               + 8.0 * C * iw2 - 48.0 * C * iw4;
        }
    }

    /// Direct Si/Ci evaluation, used when @c eps_ is outside the series' range of validity.
    double exact_G(double x) const {
        double si_f, ci_f, si_s, ci_s;
        alglib::sinecosineintegrals(x * tau_fin_, si_f, ci_f);
        alglib::sinecosineintegrals(x * tau_s_,   si_s, ci_s);
        const double dCi = ci_f - ci_s;
        const double dSi = si_f - si_s;
        return 0.25 * (dCi * dCi + dSi * dSi);
    }

    double tau_s_, tau_fin_, half_dtau_, eps_;
    bool exact_;
};

} // namespace Spectrum

#endif // INCLUDE_SSM_KERNEL_HPP_H
