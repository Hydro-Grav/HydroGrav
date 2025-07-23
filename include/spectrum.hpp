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
    // using Spectrum = std::pair<double, double>;
    using Spectrum = std::pair<std::vector<double>, std::vector<double>>;

    // ctors - pass in momentum (k) and spectrum (P) since P is calculated differently for kinetic/GW
    // PowerSpec(double k, double P);
    PowerSpec(const std::vector<double>& kvec, std::vector<double>& Pvec, const PhaseTransition::PTParams& params);

    const std::vector<double>& K() const { return data_.first; }; // Momentum
    const std::vector<double>& P() const { return data_.second; }; // Power spectrum
    const PhaseTransition::PTParams params() const { return params_; }; // PT parameters

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
    Spectrum data_;
    const PhaseTransition::PTParams params_;

};

double ptilde(double k, double p, double z);

double ff(double tau_m, double kcs);
double dtau_fin(double tau_fin, double tau_s);

std::vector<std::vector<std::vector<double>>> dlt(const int nt, const std::vector<double>& k_vals, const std::vector<double>& p_vals, const std::vector<double>& z_vals, const PhaseTransition::PTParams& params);
std::vector<std::vector<std::vector<double>>> dlt_SSM(const std::vector<double>& k_vals, const std::vector<double>& p_vals, const std::vector<double>& z_vals, const PhaseTransition::PTParams& params);
std::vector<double> dlt_SSM2(const std::vector<double>& kRs_vals, const std::vector<double>& pRs_vals, const std::vector<double>& z_vals, const PhaseTransition::PTParams& params);


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

} // namespace Spectrum

#endif // INCLUDE_SPECTRUM_HPP_H
