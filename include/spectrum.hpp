// spectrum.hpp
#ifndef INCLUDE_SPECTRUM_HPP_H
#define INCLUDE_SPECTRUM_HPP_H

#include <vector>
#include <string>
#include <variant>
#include <algorithm>

#include "phasetransition.hpp"
#include "hydrodynamics.hpp"
#include "maths_ops.hpp"

/*
TO DO:
- use template in PowerSpec ctor -> only need one ctor
- destructor for PowerSpec to save memory somewhere?
- update operator overload to use type traits (for int, unsigned int, etc. scalar)
- PowerSpec class documentation
- PowerSpec dflt ctor create it's own instance of params and profile so i don't have to keep calling it
- store k_vals adn P_vals in PowerSpec rather than calling with K(), P()? might be quicker
- simplify powerspec class to just take std::vector (not called for single k,P and these can just be stored in a vec anyway)
*/

namespace Spectrum {

/**
 * @class PowerSpec
 * @brief Represents a power spectrum, either as a scalar or a vector function of momentum.
 *
 * The PowerSpec class handles scalar and vector power spectra and provides arithmetic
 * operations and access to physical properties such as the momentum vector and maximum spectrum value.
 */
class PowerSpec {
  public:
    // ctors - pass in momentum (k) and spectrum (P) since P is calculated differently for kinetic/GW
    // PowerSpec(const std::vector<double>& K_vals, std::vector<double>& P_vals, const PhaseTransition::PTParams& params);
    PowerSpec(const std::vector<double>& K_vals, std::vector<double>& P_vals, const Hydrodynamics::FluidProfile& profile);

    const std::vector<double>& freq() const { return freq_vals_; }; // Frequency
    const std::vector<double>& K() const { return K_vals_; }; // Momentum
    const std::vector<double>& P() const { return P_vals_; }; // Power spectrum

    const Hydrodynamics::FluidProfile profile() const { return profile_; }; // fluid profile
    const PhaseTransition::PTParams* params() const { return params_; }; // PT parameters

    double max() const; // Max value of power spectrum

    void write(const std::string& filename="spectrum.csv") const;

    #ifdef ENABLE_MATPLOTLIB
    void plot(const std::string& filename="spectrum.png") const;
    #endif

    CubicSpline<double> interpolate() const; // generate cubic spline interpolation of P vals

    // Scalar arithmetic
    friend PowerSpec operator*(const PowerSpec &spec, double scalar);
    friend PowerSpec operator*(double scalar, const PowerSpec &spec);
    PowerSpec &operator*=(double scalar);

    friend PowerSpec operator/(const PowerSpec& spec, double scalar);
    PowerSpec &operator/=(double scalar);  

  private:
    std::vector<double> freq_vals_, K_vals_, P_vals_;
    const Hydrodynamics::FluidProfile profile_; // fluid profile
    const PhaseTransition::PTParams* params_;

};

inline double ptilde(double k, double p, double z) {
    const auto arg = k*k - 2.0 * k * p * z + p*p;

    if (arg < 0.0)
        throw std::runtime_error("arg<0 in Spectrum::ptilde"); // avoids numerical precision issues giving arg < 0

    if (std::abs(arg) < 1e-10)
        return 0.0; // avoids numerical precision issues giving arg < 0

    return std::sqrt(arg);
  }

double find_min_pt(const std::vector<double>& k_vals, const std::vector<double>& p_vals);

double ff(double tau_m, double kcs);

std::vector<std::vector<std::vector<double>>> dlt(const int nt, const std::vector<double>& k_vals, const std::vector<double>& p_vals, const std::vector<double>& z_vals, const PhaseTransition::PTParams& params);
std::vector<double> dlt_SSM(const std::vector<double>& kRs_vals, const std::vector<double>& pRs_vals, const std::vector<double>& z_vals, const PhaseTransition::PTParams& params);
double dlt_SSM2(double k, double p, double pt, const double cs, const double tau_s, const double tau_fin);

/**
 * @brief Calculates kinetic (velocity) power spectrum
 *
 * @param k Momentum
 * @param beta Inverse PT duration
 * @param Rs Characteristic length scale R_*
 * @param nuc_type Bubble nucleation method (exp or sim)
 *
 * @return Kinetic power spectrum
 */
PowerSpec Ekin(const std::vector<double>& k_vec, const Hydrodynamics::FluidProfile& prof);
PowerSpec Ekin(const std::vector<double>& k_vec, const PhaseTransition::PTParams& params);

/*
 * @brief Calculates normalised kinetic power spectrum from Ekin(k)
 *
 * @param Ekin kinetic power spectrum
 *
 * @return Normalised kinetic power spectrum
 */
PowerSpec norm_spec(const PowerSpec& spec);
PowerSpec zetaKin(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& prof);
PowerSpec zetaKin(const std::vector<double>& kRs_vals, const PhaseTransition::PTParams& params);


PowerSpec GWSpec(const std::vector<double>& kRs_vals, const PhaseTransition::PTParams& params);
PowerSpec GWSpec2(const std::vector<double>& kRs_vals, const PhaseTransition::PTParams& params, const size_t n_fp=5000, const size_t np=1000, const size_t nz=1000);

/**
 * @brief Calculates prefactor for GW power spectrum $\Omega_{GW}$
 *
 * @param csq The square of the speed of sound (in the broken phase)
 * @param T0 Temperature of the universe today
 * @param H0 Hubble constant today
 * @param g0 Number of dof today
 * @param gs Number of dof at beginning of FOPT
 *
 * @return Constant prefactor of eq (93) in arXiv:2308.12943
 */
double gw_prefac(double Ekin_max, double Rs, double wNeN_rat, double T0, double Ts, double H0, double Hs, double g0, double gs);

/**
 * @brief Calculates prefactor for GW power spectrum $\Omega_{GW}$
 *
 * @param csq The square of the speed of sound (in the broken phase)
 * @param u Universe parameters at start of PT and present day
 *
 * @return Constant prefactor of eq (93) in arXiv:2308.12943
 */
double gw_prefac(const std::vector<double>& kRs_vals, const Hydrodynamics::FluidProfile& profile);

#ifdef ENABLE_MATPLOTLIB
void plot_spectra(const PowerSpec& gw_spec_bag, const PowerSpec& gw_spec_munu, const PowerSpec& gw_spec_veff, const std::string& filename="gw_spec_combined", const double f_min=1e-6, const double f_max=1e+0);
#endif

} // namespace Spectrum

#endif // INCLUDE_SPECTRUM_HPP_H
