/**
 * @file wallvelocity.hpp
 * @brief Bubble-wall velocity in local thermal equilibrium (LTE) for a generic EoS.
 *
 * In LTE the particle distributions are in equilibrium, so the entropy current
 * \f$s u^\mu\f$ is conserved across the wall (Ai, Garbrecht & Tamarit,
 * arXiv:2109.13710). Combined with energy-flux conservation and \f$w = Ts\f$ this
 * gives a third junction condition (Ai, Laurent & van de Vis, arXiv:2303.10171):
 *
 * \f[ \gamma_+ T_+ = \gamma_- T_- , \f]
 *
 * which, together with the energy-momentum junction conditions, the shock
 * conditions and the fluid equations of motion already solved by
 * FluidProfile, fixes the wall velocity. The only model input is the
 * equation of state \f$p_{s,b}(T)\f$, \f$e_{s,b}(T)\f$ held by PTParams_Veff.
 *
 * The solver scans the residual
 * \f$\Delta(v_w) = \gamma(v_+)T_+/(\gamma(v_-)T_-) - 1\f$
 * over \f$v_w \in [v_{\rm min}, v_J - \delta]\f$ and refines sign changes with
 * Brent's method. Detonations are not searched: in LTE their friction decreases
 * with \f$v_w\f$, so they are unstable (arXiv:2303.10171).
 *
 * Caveats: LTE is the minimum-friction limit, so the result is an upper bound
 * on the true wall velocity (arXiv:2411.13641); close to \f$v_J\f$ the
 * steady-state assumption itself is questionable (arXiv:2402.15408).
 */
#ifndef INCLUDE_WALLVELOCITY_HPP_H
#define INCLUDE_WALLVELOCITY_HPP_H

#include <memory>
#include <string>
#include <vector>

#include "config.hpp"
#include "phasetransition.hpp"
#include "profile.hpp"

namespace Hydrodynamics {

/**
 * @brief Outcome of an LTE wall-velocity calculation.
 */
enum class LTEStatus {
    Converged,   ///< A steady-state deflagration/hybrid solution was found.
    Runaway,     ///< The driving force beats the maximal hydrodynamic obstruction (at v_J): no steady state below v_J.
    NoExpansion, ///< The driving force is too weak for any expanding steady state in the scanned range.
    Failed       ///< The residual could not be evaluated or bracketed reliably.
};

/** @brief Name of an LTEStatus value. */
std::string to_string(LTEStatus status);

/**
 * @brief Numerical options for LTEWallSolver. Defaults are taken from config.hpp.
 */
struct LTEOptions {
    double vw_min    = config::lte_vw_min;        ///< Smallest wall velocity in the scan.
    double vJ_margin = config::lte_vJ_margin;     ///< Scan stops at v_J - vJ_margin.
    double root_tol  = config::lte_root_tol;      ///< Absolute tolerance on the returned v_w.
    int scan_points  = config::lte_scan_points;   ///< Number of scan points.
    size_t profile_n = config::lte_profile_steps; ///< Integration steps for each trial profile.
};

/**
 * @brief Result of LTEWallSolver::solve().
 */
struct LTEResult {
    LTEStatus status = LTEStatus::Failed; ///< Outcome.
    double vw = std::numeric_limits<double>::quiet_NaN(); ///< LTE wall velocity (NaN unless Converged).
    int mode = -1;                        ///< Hydrodynamic mode at the solution (0 = deflagration, 1 = hybrid).
    double vJ = std::numeric_limits<double>::quiet_NaN(); ///< Jouguet velocity (upper end of the search).
    WallState wall;                       ///< Fluid state either side of the wall at the solution.
    std::vector<double> vw_scan;          ///< Scanned wall velocities.
    std::vector<double> resi_scan;        ///< \f$\Delta(v_w)\f$ at each scanned velocity (NaN where evaluation failed).

    /**
     * @brief Copy of the input parameters with the solved wall velocity.
     *
     * Rs (and Rbar, if supplied) are rescaled with v_w unless solve() was told
     * otherwise. Null unless status is Converged. FluidProfile and PowerSpec
     * keep a raw pointer to their parameters, so this object must outlive any
     * profile or spectrum built from it.
     */
    std::shared_ptr<PhaseTransition::PTParams_Veff> params;
};

/**
 * @class LTEWallSolver
 * @brief Computes the bubble-wall velocity in local thermal equilibrium from a tabulated EoS.
 *
 * @code
 * PhaseTransition::PTParams_Veff base(0.5, alN, TN, beta, Rs, Rbar, nuc_type, un, eos); // vw is a placeholder
 * Hydrodynamics::LTEWallSolver lte(base);
 * const auto res = lte.solve();
 * if (res.status == Hydrodynamics::LTEStatus::Converged) {
 *     Hydrodynamics::FluidProfile profile(*res.params);
 * }
 * @endcode
 */
class LTEWallSolver {
  public:
    /**
     * @param base    Veff parameters. Their wall velocity is ignored (any valid value will do).
     * @param options Numerical options.
     */
    explicit LTEWallSolver(const PhaseTransition::PTParams_Veff& base, LTEOptions options = {});

    /**
     * @brief Solve for the LTE wall velocity.
     *
     * @param rescale_lengths If true, the returned parameters have Rs (and Rbar,
     *                        if supplied) rescaled to the solved wall velocity.
     */
    LTEResult solve(bool rescale_lengths = true) const;

    /**
     * @brief Entropy-matching residual \f$\Delta(v_w) = \gamma(v_+)T_+/(\gamma(v_-)T_-) - 1\f$.
     *
     * Solves the fluid profile for the trial wall velocity (energy-momentum and
     * shock conditions only) and evaluates the LTE condition on its wall state.
     *
     * The residual is only as precise as the inner shooting on \f$T_-\f$ in
     * FluidProfile: it is piecewise smooth in \f$v_w\f$ with jumps of
     * \f$O(10^{-5})\f$. In tests against the reference code of arXiv:2303.10171
     * the resulting wall velocity agrees to better than \f$10^{-4}\f$.
     *
     * @return The residual, or NaN if the profile could not be solved
     *         (this happens for slow walls, \f$v_w \lesssim 0.1\f$).
     */
    double entropy_residual(double vw) const;

    /** @brief Jouguet velocity of the EoS (computed once at construction). */
    double jouguet_velocity() const { return vJ_; }

    /**
     * @brief Largest \f$|w/(T\,dp/dT) - 1|\f$ over the EoS table, in either phase.
     *
     * The LTE condition relies on \f$w = Ts\f$ with \f$s = dp/dT\f$.
     */
    double max_entropy_inconsistency() const;

    /** @brief Write the scanned \f$\Delta(v_w)\f$ curve to a CSV file. */
    static void write_scan(const LTEResult& result, const std::string& filename = "lte_scan.csv");

  private:
    /// Trial-profile evaluation at one wall velocity.
    struct Evaluation {
        double resi = std::numeric_limits<double>::quiet_NaN();
        WallState wall;
        int mode = -1;
    };

    Evaluation evaluate(double vw) const;

    const PhaseTransition::PTParams_Veff& base_;
    LTEOptions opts_;
    double vJ_;
};

} // namespace Hydrodynamics

#endif // INCLUDE_WALLVELOCITY_HPP_H
