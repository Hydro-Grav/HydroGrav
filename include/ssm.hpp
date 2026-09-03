/**
 * @file ssm.hpp
 * @brief Gravitational wave calculator using the sound shell model (see arXiv:1608.04735)
 */

#ifndef INCLUDE_SSM_HPP_H
#define INCLUDE_SSM_HPP_H

#include <vector>
#include <string>

#include "phasetransition.hpp"
#include "profile.hpp"
#include "maths.hpp"
#include "ssm_kernel.hpp"

namespace Hydrodynamics {

  struct

/**
 * @brief Lifetime distribution function of bubbles
 *
 * Supported values of @p nuc_type are:
 *   - "exp" – simple exponential decay: \f$e^{-\tilde T}\f$
 *   - "sim" – simulated form \f$0.5\tilde T^2e^{-\tilde T^3/6}\f$
 */
std::function<double(double)> lifetime_distribution_function(const std::string& nuc_type);

/**
 * @brief Evaluate fluid profile integrals \f$f'(\chi)\f$ and \f$l(\chi)\f$
 *
 * Given a vector of momentum arguments @p chi_vals and a
 * @c FluidProfile object, this function computes the integrals over the fluid profiles.
 */
std::pair<std::vector<double>, std::vector<double>>
fluid_profile_integrals(const std::vector<double>& chi_vals, const FluidProfile& prof);

/**
 * @brief Compute squared amplitude \f$|A_+|^2\f$ for a profile.
 *
 * This helper forms the combination
 * \f$0.25(f^2 + c_s^2 l^2)\f$, where \f$c_s^2\f$ is the background sound
 * speed stored in the profile and \f$f',l\f$ originate from
 * @c fluid_profile_integrals().
 */
std::vector<double> Ap_sq(const std::vector<double>& chi_vals, const FluidProfile& prof);

/**
 * @class ApsqSpline
 * @brief Cubic spline of \f$|A_+|^2(\chi)\f$ over the @c config chi grid.
 *
 * Building it evaluates the oscillatory fluid-profile integrals at every grid point, which
 * dominates the cost of @c Spectrum::Ekin().  A single @c GWSpec() call needs @c Ekin() for
 * three different momentum grids of the same profile, so it builds one of these up front and
 * passes it to each, instead of repeating the profile integrals.
 */
class ApsqSpline {
  public:
    explicit ApsqSpline(const FluidProfile& prof);

    /// Interpolated \f$|A_+|^2\f$ at momentum @p chi.
    double operator()(double chi) const { return alglib::spline1dcalc(spline_, chi); }

  private:
    alglib::spline1dinterpolant spline_;
};

} // namespace Hydrodynamics

/**
 * @namespace Spectrum
 * @brief Sound shell model implementation of gravitational wave spectra
 */
namespace Spectrum {

/**
 * @class PowerSpec
 * @brief Encapsulates a power spectrum as a function of momentum.
 *
 * The class stores lists of momentum values @c K and corresponding spectrum
 * values @c P together with derived frequency values.  It provides arithmetic
 * operators, file I/O, and accessors for peak frequency/amplitude properties.
 */
class PowerSpec {
  public:
    /**
     * @brief Construct from momentum and spectrum vectors.
     *
     * @param K_vals   Momentum samples (dimensionless K=k*Rs).
     * @param P_vals   Spectrum values at each momentum point.
     * @param profile  Fluid profile used to generate the spectrum.
     * @param dtau     Sound–shell duration.
     *
     * Throws std::invalid_argument if @p K_vals and @p P_vals differ in size.
     */
    PowerSpec(const std::vector<double>& K_vals, std::vector<double>& P_vals, const Hydrodynamics::FluidProfile& profile, double dtau=std::numeric_limits<double>::quiet_NaN());

    /// Frequency values corresponding to the stored momenta.
    const std::vector<double>& freq() const { return freq_vals_; };
    /// Momentum (kRs) values.
    const std::vector<double>& K() const { return K_vals_; };
    /// Power spectrum values.
    const std::vector<double>& P() const { return P_vals_; };

    /// Fluid profile used to compute this spectrum.
    const Hydrodynamics::FluidProfile& profile() const { return profile_; };
    /// Phase transition parameters pointer (owned externally).
    const PhaseTransition::PTParams* params() const { return params_; };

    /// Sound‑wave duration parameter.
    double dtau() const { return dtau_; };

    /**
     * @brief Return peak frequency and amplitude of the spectrum.
     *
     * @return Pair <f_peak, A_peak>.
     */
    std::pair<double, double> peak_vals() const;

    /**
     * @brief Write spectrum to a CSV file with header information.
     *
     * @param filename Path to output file (default "spectrum.csv").
     */
    void write(const std::string& filename="spectrum.csv", const bool write_header=false) const;

    #ifdef ENABLE_MATPLOTLIB
    /**
     * @brief Plot the spectrum using matplotlibcpp.
     *
     * @param filename Image filename (e.g. "spectrum.png").
     */
    void plot(const std::string& filename="spectrum.png") const;
    #endif

    // Scalar arithmetic operators ------------------------------------------------
    friend PowerSpec operator*(const PowerSpec &spec, double scalar);
    friend PowerSpec operator*(double scalar, const PowerSpec &spec);
    PowerSpec &operator*=(double scalar);

    friend PowerSpec operator/(const PowerSpec& spec, double scalar);
    PowerSpec &operator/=(double scalar);

  private:
    std::vector<double> freq_vals_, K_vals_, P_vals_;
    double dtau_;
    const Hydrodynamics::FluidProfile profile_;
    const PhaseTransition::PTParams* params_;
};

/// Helper that computes \f$\tilde p=\sqrt{k^2-2kpz+p^2}\f$
inline double ptilde(double k, double p, double z) {
  const auto arg = k*k - 2.0 * k * p * z + p*p;

  if (std::abs(arg) < 1e-10)
    return 0.0; // avoids numerical precision issues giving arg < 0

  return std::sqrt(arg);
}

/**
 * @brief Find minimum \f$p_t\f$ over a grid of k and p values.
 */
double find_min_pt(const std::vector<double>& k_vals, const std::vector<double>& p_vals);

/**
 * @brief Computes the sound-shell model time-correlation function \f$\Delta_{\rm SSM}\f$.
 *
 * Evaluates the double-time integral of the velocity–velocity correlator
 * that enters the GW power spectrum integrand.
 *
 * @param k       Wavenumber of emitted GW.
 * @param p       Intermediate momentum.
 * @param pt      Transferred momentum \f$|\mathbf{k}-\mathbf{p}|\f$.
 * @param cs      Background speed of sound.
 * @param tau_s   Conformal time at the start of the phase transition.
 * @param tau_fin Conformal time at the end of the sound-wave period.
 * @return Value of \f$\Delta(k,p,\tilde p)\f$.
 */
double dlt_SSM(double k, double p, double pt, const double cs, const double tau_s, const double tau_fin);

/**
 * @brief Hot-loop overload of @c dlt_SSM taking a prebuilt @c SoundShellKernel.
 *
 * @param kernel  Kernel built once from (tau_s, tau_fin).
 * @param sin_wk  sin(k * kernel.half_dtau()), precomputed per k.
 * @param cos_wk  cos(k * kernel.half_dtau()), precomputed per k.
 */
double dlt_SSM(double k, double p, double pt, const double cs,
               const SoundShellKernel& kernel, double sin_wk, double cos_wk);

/**
 * @brief Calculates kinetic (velocity) power spectrum.
 *
 * Overloaded versions accept either a @c FluidProfile or a @c PTParams object
 */
PowerSpec Ekin(const std::vector<double>& k_vec, const Hydrodynamics::FluidProfile& prof);
PowerSpec Ekin(const std::vector<double>& k_vec, const PhaseTransition::PTParams& params);
/// Overload reusing a prebuilt \f$|A_+|^2\f$ spline; see @c Hydrodynamics::ApsqSpline.
PowerSpec Ekin(const std::vector<double>& k_vec, const Hydrodynamics::FluidProfile& prof,
               const Hydrodynamics::ApsqSpline& apsq);

/**
 * @brief Normalise a kinetic spectrum to unit integral.
 */
PowerSpec norm_spec(const PowerSpec& spec);

/**
 * @brief Normalised kinetic spectrum.
 * 
 * Overloaded versions accept either a @c FluidProfile or a @c PTParams object
 */
PowerSpec zetaKin(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& prof);
PowerSpec zetaKin(const std::vector<double>& kRs_vals, const PhaseTransition::PTParams& params);
/// Overload reusing a prebuilt \f$|A_+|^2\f$ spline; see @c Hydrodynamics::ApsqSpline.
PowerSpec zetaKin(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& prof,
                  const Hydrodynamics::ApsqSpline& apsq);

/**
 * @brief Calculates gravitational wave spectrum.
 */
PowerSpec GWSpec(const std::vector<double>& kRs_vals, const PhaseTransition::PTParams& params, double dtau=0.0);

/**
 * @brief Builds spline interpolating functions for kinetic spectrum
 */
void build_kinetic_spectrum_spline(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& profile, alglib::spline1dinterpolant& log_zk_spline);
/// Overload reusing a prebuilt \f$|A_+|^2\f$ spline; see @c Hydrodynamics::ApsqSpline.
void build_kinetic_spectrum_spline(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& profile,
                                   const Hydrodynamics::ApsqSpline& apsq, alglib::spline1dinterpolant& log_zk_spline);

/**
 * @brief Builds an O(1) interpolant of the normalised kinetic spectrum.
 *
 * @p kRs_vals must be logarithmically spaced.  Unlike the ALGLIB overloads this stores
 * \f$\zeta_{\rm kin}\f$ itself rather than its logarithm, so evaluation costs one @c log and
 * no @c exp -- the form @c GWSpec needs in its innermost loop.
 */
LogGridInterpolant kinetic_spectrum_interpolant(const std::vector<double>& kRs_vals,
                                                const Hydrodynamics::FluidProfile& profile,
                                                const Hydrodynamics::ApsqSpline& apsq);

/**
 * @brief Calculates the non-linear timescale of the phase transition using the time for the plasma to develop turbulence
 */
double get_nl_timescale(const Hydrodynamics::FluidProfile& prof);
/// Overload reusing a prebuilt \f$|A_+|^2\f$ spline; see @c Hydrodynamics::ApsqSpline.
double get_nl_timescale(const Hydrodynamics::FluidProfile& prof, const Hydrodynamics::ApsqSpline& apsq);

/// Approximation used for dtau in arXiv:2308.12943.
// double dtau_approx(const PhaseTransition::PTParams& params);

/**
 * @brief Prefactor for gravitational-wave spectrum (Eq. 93 in arXiv:2308.12943).
 *
 * This overload takes explicit cosmological parameters.
 */
double gw_prefac(double Ekin_max, double Rs, double wNeN_rat, double T0, double Ts, double H0, double Hs, double g0, double gs);

/**
 * @brief Prefactor for gravitational-wave spectrum using a profile.
 */
double gw_prefac(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& profile);
/// Overload reusing a prebuilt \f$|A_+|^2\f$ spline; see @c Hydrodynamics::ApsqSpline.
double gw_prefac(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& profile,
                 const Hydrodynamics::ApsqSpline& apsq);

#ifdef ENABLE_MATPLOTLIB
/**
 * @brief Plot multiple GW spectra on a single log-log axes.
 */
void plot_spectra(const PowerSpec& gw_spec_bag, const PowerSpec& gw_spec_munu, const PowerSpec& gw_spec_veff, const std::string& filename="gw_spec_combined", const double f_min=1e-6, const double f_max=1e+0);
#endif

} // namespace Spectrum

#endif // INCLUDE_SSM_HPP_H
