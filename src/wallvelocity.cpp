// wallvelocity.cpp
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <stdexcept>

#include "wallvelocity.hpp"
#include "logger.hpp"
#include "maths.hpp"
#include "physics.hpp"

namespace Hydrodynamics {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

/*
 * Sign of Delta(vw) = gamma_+ T_+ / (gamma_- T_-) - 1 when the driving force
 * exceeds the LTE (hydrodynamic) obstruction, i.e. when the wall would keep
 * accelerating. Fixed empirically on the template EoS (see
 * unit_tests/test_wallvelocity.cpp): Delta > 0 at small vw and changes sign
 * just below vJ, where the obstruction peaks; a transition too strong for a
 * steady state keeps Delta > 0 up to vJ, one too weak has Delta < 0 throughout.
 * The physical (stable) root is therefore the + -> - crossing.
 */
constexpr double kDrivingSign = +1.0;

bool same_sign(double a, double b) { return (a > 0.0) == (b > 0.0); }

} // namespace

std::string to_string(LTEStatus status) {
    switch (status) {
        case LTEStatus::Converged:   return "converged";
        case LTEStatus::Runaway:     return "runaway";
        case LTEStatus::NoExpansion: return "no expansion";
        case LTEStatus::Failed:      return "failed";
    }
    return "unknown";
}

LTEWallSolver::LTEWallSolver(const PhaseTransition::PTParams_Veff& base, LTEOptions options)
    : base_(base), opts_(options), vJ_(jouguet_velocity_veff(base)) {
    if (opts_.scan_points < 2) throw std::invalid_argument("LTEWallSolver: scan_points must be >= 2");
    if (!(opts_.vw_min > 0.0)) throw std::invalid_argument("LTEWallSolver: vw_min must be > 0");
    if (opts_.vw_min >= vJ_ - opts_.vJ_margin) {
        throw std::invalid_argument("LTEWallSolver: vw_min must lie below vJ - vJ_margin (vJ=" + std::to_string(vJ_) + ")");
    }
}

LTEWallSolver::Evaluation LTEWallSolver::evaluate(double vw) const {
    Evaluation ev;

    // Trial profiles report fallbacks as warnings; they are expected during the scan.
    logging::ScopedLevel quiet(logging::Level::Error);
    try {
        const FluidProfile fp(base_, vw, opts_.profile_n, false);
        const auto& ws = fp.wall_state();

        ev.wall = ws;
        ev.mode = fp.mode();
        // gamma_+ T_+ / (gamma_- T_-) - 1, with gamma^2 = 1/(1 - v^2)
        ev.resi = std::sqrt(gammaSq(ws.vp) / gammaSq(ws.vm)) * ws.TpTN / ws.TmTN - 1.0;
        if (!std::isfinite(ev.resi)) ev.resi = kNaN;
    } catch (const std::exception&) {
        ev.resi = kNaN;
    }
    return ev;
}

double LTEWallSolver::entropy_residual(double vw) const {
    return evaluate(vw).resi;
}

double LTEWallSolver::max_entropy_inconsistency() const {
    // Skip the outermost 5% of the table, where the smoothing spline's derivative is least reliable.
    const auto lo = base_.TTN_min();
    const auto hi = base_.TTN_max();
    const auto pad = 0.05 * (hi - lo);
    double worst = 0.0;
    for (const auto TTN : linspace(lo + pad, hi - pad, 100)) {
        worst = std::max(worst, std::abs(base_.entropy_consistency(TTN, false)));
        worst = std::max(worst, std::abs(base_.entropy_consistency(TTN, true)));
    }
    return worst;
}

LTEResult LTEWallSolver::solve(bool rescale_lengths) const {
    LTEResult res;
    res.vJ = vJ_;

    HG_LOG(Info) << "Solving for the LTE wall velocity (alN_munu=" << base_.alN_munu()
                 << ", psiN=" << base_.psiN() << ", vJ=" << vJ_ << ")...";

    const auto inconsistency = max_entropy_inconsistency();
    if (inconsistency > config::lte_entropy_consistency_tol) {
        HG_LOG(Warn) << "EoS is not thermodynamically consistent: max |w/(T dp/dT) - 1| = " << inconsistency
                     << ". The LTE condition gamma+ T+ = gamma- T- assumes w = T dp/dT, so the wall velocity may be unreliable.";
    }

    // Every trial profile costs a full shooting solve, so never evaluate the same vw twice:
    // the Brent bracket ends come from the scan and the root's wall state from Brent's last step.
    std::map<double, Evaluation> cache;
    auto eval = [this, &cache](double v) -> const Evaluation& {
        auto it = cache.find(v);
        if (it == cache.end()) it = cache.emplace(v, evaluate(v)).first;
        return it->second;
    };

    // 1. scan Delta(vw) below the Jouguet velocity
    res.vw_scan = linspace(opts_.vw_min, vJ_ - opts_.vJ_margin, static_cast<size_t>(opts_.scan_points));
    res.resi_scan.resize(res.vw_scan.size());
    for (size_t i = 0; i < res.vw_scan.size(); ++i) {
        res.resi_scan[i] = eval(res.vw_scan[i]).resi;
        HG_LOG(Debug) << "LTE scan: vw=" << res.vw_scan[i] << ", Delta=" << res.resi_scan[i];
    }

    // 2. bracket sign changes between neighbouring finite points
    std::vector<std::pair<double, double>> brackets;
    size_t n_finite = 0;
    double last_vw = kNaN, last_resi = kNaN;
    for (size_t i = 0; i < res.vw_scan.size(); ++i) {
        const auto r = res.resi_scan[i];
        if (!std::isfinite(r)) continue;
        ++n_finite;
        if (std::isfinite(last_resi) && !same_sign(last_resi, r)) {
            brackets.emplace_back(last_vw, res.vw_scan[i]);
        }
        last_vw = res.vw_scan[i];
        last_resi = r;
    }

    if (n_finite < 2) {
        HG_LOG(Error) << "LTE wall velocity: the fluid profile could not be solved over the scan range.";
        res.status = LTEStatus::Failed;
        return res;
    }

    // 3. no sign change: classify by the sign of Delta
    if (brackets.empty()) {
        double sum = 0.0;
        for (const auto r : res.resi_scan) if (std::isfinite(r)) sum += r;
        res.status = (sum * kDrivingSign > 0.0) ? LTEStatus::Runaway : LTEStatus::NoExpansion;
        if (res.status == LTEStatus::Runaway) {
            HG_LOG(Warn) << "No LTE steady state below vJ=" << vJ_
                         << ": the driving force exceeds the maximal hydrodynamic obstruction (runaway in LTE).";
        } else {
            HG_LOG(Warn) << "No LTE steady state above vw=" << opts_.vw_min
                         << ": the driving force is too weak to sustain an expanding wall.";
        }
        return res;
    }

    if (brackets.size() > 1) {
        HG_LOG(Warn) << "LTE residual changes sign " << brackets.size()
                     << " times below vJ; taking the lowest root. Inspect the scan (write_scan) to check.";
    }

    // 4. refine the lowest root
    const auto [a, b] = brackets.front();
    double vw = kNaN;
    try {
        auto f = [&eval](double v) { return eval(v).resi; };
        vw = brent_root(f, a, b, opts_.root_tol);
    } catch (const std::exception& e) {
        // a NaN inside the bracket breaks Brent's method: fall back to linear interpolation of the scan
        const auto fa = eval(a).resi;
        const auto fb = eval(b).resi;
        vw = a - fa * (b - a) / (fb - fa);
        HG_LOG(Warn) << "Brent refinement of the LTE root failed (" << e.what()
                     << "); using linear interpolation between scan points, vw=" << vw;
    }

    const auto& ev = eval(vw);
    HG_LOG(Debug) << "LTE wall velocity: " << cache.size() << " trial profiles solved";
    if (!std::isfinite(ev.resi)) {
        HG_LOG(Error) << "LTE wall velocity: fluid profile failed at the refined root vw=" << vw;
        res.status = LTEStatus::Failed;
        return res;
    }

    res.status = LTEStatus::Converged;
    res.vw = vw;
    res.mode = ev.mode;
    res.wall = ev.wall;
    res.params = std::shared_ptr<PhaseTransition::PTParams_Veff>(base_.with_vw(vw, rescale_lengths));

    HG_LOG(Info) << "LTE wall velocity: vw=" << vw << " (" << (ev.mode == 0 ? "deflagration" : "hybrid")
                 << ", vw/vJ=" << vw / vJ_ << ", residual=" << ev.resi << ")";
    if (vw / vJ_ > config::lte_near_vJ_warn) {
        HG_LOG(Warn) << "LTE wall velocity is close to the Jouguet velocity (vw/vJ=" << vw / vJ_
                     << "). The steady-state assumption may not hold here (see arXiv:2402.15408).";
    }

    return res;
}

void LTEWallSolver::write_scan(const LTEResult& result, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) throw std::runtime_error("LTEWallSolver::write_scan: could not open " + filename);

    file << "# LTE wall velocity scan: Delta = gamma+ T+ / (gamma- T-) - 1\n"
         << "# status = " << to_string(result.status) << ", vw = " << result.vw << ", vJ = " << result.vJ << "\n"
         << "vw,Delta\n";
    file << std::setprecision(10);
    for (size_t i = 0; i < result.vw_scan.size(); ++i) {
        file << result.vw_scan[i] << "," << result.resi_scan[i] << "\n";
    }

    HG_LOG(Info) << "Wrote LTE wall-velocity scan to " << filename;
}

} // namespace Hydrodynamics
