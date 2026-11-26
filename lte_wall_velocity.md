# Design: LTE bubble-wall velocity solver for HydroGrav (Veff EoS)

*Branch: `feature/vw`. Revised from a review of `feature/inverse_transitions`: see "Changes for feature/vw" below.*

## Context
HydroGrav takes the wall velocity `vw` as an input (`PTParams::vw_`) and builds the self-similar fluid profile from it ([FluidProfile](include/profile.hpp#L120)). The goal is to **compute** `vw` from the generic tabulated EoS `p_{s,b}(T), e_{s,b}(T)` (`PTParams_Veff`), using **only the local-thermal-equilibrium (LTE) approximation**. Out-of-equilibrium (Boltzmann) friction is out of scope. Agreed scope: **Veff EoS only**; bag/μν cases are tested by feeding a tabulated template EoS.

Main conclusion: in LTE, `vw` follows from the EoS through **one extra junction condition at the wall (entropy-flux conservation)**. HydroGrav already solves everything else that is needed (energy-momentum junction conditions, shock matching, profile ODEs with T(ξ)). The feature is therefore a 1-D root-find in `vw` wrapped around the existing `FluidProfile` machinery, plus small refactors to expose the wall state.

### Changes for feature/vw (vs. the first review)
- **Logger** (`include/logger.hpp`): output is silenced during the vw scan with `logging::ScopedLevel(logging::Level::Error)`. A custom stream silencer is no longer needed.
- **`Rbar`** now exists (`Rbar_approx(vw, β) = vw/β`). `dR() = |vw − c_s|` and `Rsound()` read `vw_` directly, and `get_nl_timescale` uses `Rbar()`/`Rsound()` ([ssm.cpp:580](src/ssm.cpp#L580)). **Both Rs and Rbar must be rescaled** with the solved vw. Rbar is rescaled only when it was supplied, since `-1` is the "not given" sentinel at [phasetransition.hpp:259](include/phasetransition.hpp#L259).
- **Inverse PTs** are gone on this branch (`PTParams` now requires alN > 0). They become a future-work note only.
- **Unit tests** are globbed (`unit_tests/CMakeLists.txt`), so a new `test_*.cpp` needs no CMake edit. Library sources are globbed too. Only the explicit `HEADERS` list in the root `CMakeLists.txt` needs the new header.
- **Still true:** no 1-D root finder exists in `maths.hpp`; the wall state is still only logged (Debug), not stored; v_J is still computed inside `get_mode_veff` and discarded.

---

## 1. New physics

### 1.1 What the code already enforces (unchanged)
At the wall (wall frame, + = symmetric/front, − = broken/behind), energy-momentum conservation gives
$$w_+\gamma_+^2 v_+ = w_-\gamma_-^2 v_-, \qquad w_+\gamma_+^2 v_+^2 + p_+ = w_-\gamma_-^2 v_-^2 + p_-,$$
implemented in combined form in [`matching_eqs_wall`](src/profile.cpp#L808). That gives two equations in four wall unknowns (v±, T±). Fixing `vw` closes the system: v₋ = vw for deflagrations, v₋ = c₋ for hybrids, T₊ = T_N and v₊ = vw for detonations, plus the shock conditions. This is why `vw` is currently an input. [Espinosa et al. 2010; Giese, Konstandin, van de Vis 2020]

### 1.2 The LTE condition: entropy conservation across the wall
In LTE (δf = 0) the scalar EoM is ∂²φ + V'_eff(φ,T) = 0, and the entropy current S^μ = s u^μ is conserved [Ai, Garbrecht, Tamarit, arXiv:2109.13710]. Integrated across a stationary wall:
$$s_+\gamma_+ v_+ = s_-\gamma_- v_- .$$
Divide the energy-flux condition by this and use w = Ts to get [Ai, Laurent, van de Vis, arXiv:2303.10171, Eqs. (14)–(15)]
$$\boxed{\;\gamma_+ T_+ = \gamma_- T_-\;}$$
This third junction condition removes the freedom in `vw`. Energy-momentum + entropy + shock conditions + fluid EoMs then (generically) have a unique deflagration/hybrid solution, which gives the LTE `vw`. The friction is **hydrodynamic obstruction**: the plasma heats in front of the wall, and the inhomogeneous temperature produces a backreaction pressure [Konstandin & No, arXiv:1011.3735; Balaji, Spannowsky, Tamarit, arXiv:2010.08013]. Only p(T), e(T) in each phase are needed, which is exactly what `PTParams_Veff` holds.

**Assumption needed:** thermodynamic consistency, s = dp/dT and w = Ts. The existing ODE `dT/dv = Tγ²μ` ([profile.cpp:62](src/profile.cpp#L62)) already assumes it. The γT form is used instead of s±(T) from spline derivatives because it needs no noisy derivative and T_N cancels.

### 1.3 Mode structure in LTE
From 2303.10171:
- **Deflagrations/hybrids** are the physical steady states. The hybrid has v₋ = c₋, and the deflagration→hybrid transition is continuous. HydroGrav's c₋ = c_b(T₋) ([profile.cpp:1492](src/profile.cpp#L1492)) is consistent with this.
- **Detonations** in LTE have friction that decreases with vw. They are unstable and "not physical", and whenever one exists a deflagration/hybrid also exists. ⇒ **search vw ∈ (0, v_J) only.**
- Obstruction is **maximal at v_J** [Ai, Nagels, Vanvlasselaer, arXiv:2401.05911]. If no root lies below v_J, the result is **runaway** in LTE.
- For α_N below α_min^def (template: max[(1−Ψ_N)/3, (μ−ν)/(3μ)]) there is **no expanding solution**.

### 1.4 Template parameters (diagnostics and validation)
For c_s = c_b = 1/√3, LTE vw depends only on α_N and Ψ_N = w_b(T_N)/w_s(T_N), with α = Δθ̄/(3w_s) and θ̄ = e − p/c_b². The paper's fit (Eqs. 32–34) is
ξ^low = √[(3α_N+Ψ_N−1)/(2(2−3Ψ_N+Ψ_N³))], ξ^high = ξ_J(1 − a(1−Ψ_N)^b/α_N), ξ^fit = ‖(ξ^low, ξ^high)‖_p,
with a = 0.2233, b = 1.704, p = −3.433. It is used only as a test oracle.

### 1.5 Caveats to state in docs and logs
- LTE is the minimum-friction limit, so the LTE vw is an **upper bound** (arXiv:2411.13641).
- The steady state near v_J is questionable (arXiv:2402.15408, 2606.30740). Warn when vw/v_J > 0.95.
- Entropy is produced at the shock. The LTE condition is imposed **only at the wall**.
- Future work: inverse PTs obey the same condition (arXiv:2406.01596).

---

## 2. Code-review findings relevant to this feature (feature/vw)
1. The wall state (v±, T±, w±) is only logged at Debug ([profile.cpp:1554](src/profile.cpp#L1554), [1596](src/profile.cpp#L1596), [1654](src/profile.cpp#L1654)) → store it.
2. `vw_` is taken from params in the ctor ([profile.cpp:207](src/profile.cpp#L207)), and `ssm.cpp:94` reads `prof.params()->vw()` → trial profiles need a vw override, and the final profile needs params carrying the solved vw.
3. v_J is computed inside [`get_mode_veff`](src/profile.cpp#L712) and discarded → extract it.
4. Rs, Rbar, dR and Rsound depend on vw → rescale them on copy.
5. Doc bug: [profile.hpp:392](include/profile.hpp#L392) describes `matching_eqs_wall` as "p₊ = p₋, T₊v₊ = T₋v₋". It actually imposes the energy-momentum conditions.
6. No check of w = T dp/dT. The smoothing spline ([phasetransition.cpp:348](src/phasetransition.cpp#L348)) and the wN/eN warning ([:415](src/phasetransition.cpp#L415)) show the EoS can be off.
7. Minor, pre-existing: the ctor tests `mode_ == 0 || mode_ == 1` at [profile.cpp:226](src/profile.cpp#L226) before `mode_` is assigned. It is benign (α_min evaluates to 0) and noted only.
8. Bag-path T(w) relations hard-code a₊/a₋ = 1. This is out of scope (Veff only).

---

## 3. Implementation plan (order of work)

**Step 0.** Save this plan into the repo as `lte_wall_velocity.md` (root, alongside `optimisation.md`).

**Step 1: maths.** Add `double brent_root(const std::function<double(double)>& f, double a, double b, double tol, int max_iter)` to [maths.hpp](include/maths.hpp)/[maths.cpp](src/maths.cpp), wrapping `gsl_root_fsolver_brent` (GSL is already linked). Add a unit test in `test_maths.cpp`.

**Step 2: config.** Add to [config.hpp](include/config.hpp):
- `lte_vw_min = 1e-2`
- `lte_vJ_margin = 1e-3`
- `lte_scan_points = 40`
- `lte_root_tol = 1e-6`
- `lte_profile_steps = 1000`
- `lte_entropy_consistency_tol = 1e-2`
- `lte_near_vJ_warn = 0.95`

**Step 3: PTParams** ([phasetransition.hpp](include/phasetransition.hpp)/[.cpp](src/phasetransition.cpp)).
- `PTParams`: add a protected `void set_vw(double vw, bool rescale_lengths)`. It validates 0 < vw < 1. When rescaling it sets `Rs_ = Rs_approx(vw, β)` and, if `Rbar_ != -1`, `Rbar_ = Rbar_approx(vw, β)`.
- `PTParams_Veff`:
  - `std::unique_ptr<PTParams_Veff> with_vw(double vw, bool rescale_lengths = true) const` (copy-construct, then `set_vw`).
  - `double alN_munu() const` and `double psiN() const`, using the commented θ̄ code at [phasetransition.cpp:426](src/phasetransition.cpp#L426).
  - `double entropy_consistency(double TTN, bool broken) const`, returning w/(T_N·TTN·dp/dT) − 1 via `spline1ddiff`.

**Step 4: profile** ([profile.hpp](include/profile.hpp)/[profile.cpp](src/profile.cpp)).
- Add `struct WallState { double vp, vm, TpTN, TmTN, wpwN, wmwN, xi_sh; };` (wall-frame velocities), a member `wall_state_`, and a public `wall_state()`. Fill it in `solve_profile_veff` for all three modes from values already in scope. Bag path: fill v± and w± only, and set T to NaN.
- Turn the public ctor into a delegating ctor that calls a new **private** `FluidProfile(const PTParams&, double vw, size_t n, bool dev_log)`, which holds the current body and initialises `vw_` from its argument. Declare `friend class LTEWallSolver;`.
- Move the Veff wall/shock matching and Jouguet logic into **free functions** taking `const PTParams_Veff&`: `veff_matching_eqs_wall(...)` and `double jouguet_velocity_veff(const PTParams_Veff&)`. The existing members become thin wrappers, and `get_mode_veff` calls `jouguet_velocity_veff`. This lets the solver get v_J before building any profile.
- Fix the doc comment (finding 5).

**Step 5: new `include/wallvelocity.hpp`, `src/wallvelocity.cpp`** (namespace `Hydrodynamics`). Add the header to `hydrograv.hpp` and to the `HEADERS` list in `CMakeLists.txt`.
```cpp
enum class LTEStatus { Converged, Runaway, NoExpansion, Failed };
struct LTEOptions { double vw_min, vJ_margin, root_tol; int scan_points; size_t profile_n; }; // defaults from config
struct LTEResult {
  LTEStatus status; double vw; int mode; double vJ; WallState wall;
  std::vector<double> vw_scan, resi_scan;                    // Δ(vw) curve
  std::shared_ptr<PhaseTransition::PTParams_Veff> params;    // input copy with solved vw (+ rescaled Rs/Rbar)
};
class LTEWallSolver {
 public:
  explicit LTEWallSolver(const PhaseTransition::PTParams_Veff& base, LTEOptions = {});
  LTEResult solve(bool rescale_lengths = true) const;
  double entropy_residual(double vw) const;   // Δ = γ(v+)T+/(γ(v-)T-) − 1, NaN on failure
  double jouguet_velocity() const;
  static void write_scan(const LTEResult&, const std::string& file = "lte_scan.csv");
  static std::string status_str(LTEStatus);
};
std::string to_string(LTEStatus);
```
`solve()`:
1. Run the EoS consistency check over the table range. Warn if it fails.
2. Compute v_J.
3. Scan Δ on a grid in [vw_min, v_J − margin] under `logging::ScopedLevel(Error)`, using the private ctor with n = `profile_n` and catching exceptions as NaN.
4. Bracket sign changes over adjacent finite points and refine each with `brent_root`.
5. Classify: no sign change → `Runaway` or `NoExpansion` according to the sign of Δ, determined empirically on the template EoS (expected: Δ keeps the "driving" sign up to v_J → Runaway). Several roots → take the lowest one and warn.
6. Build `params` via `with_vw`, then log a summary at Info, plus a warning if vw/v_J > 0.95.

`shared_ptr` is used so `LTEResult` stays copyable. `FluidProfile` holds a raw pointer, so `result.params` must outlive any profile or spectrum built from it.

**Step 6: test** `unit_tests/test_wallvelocity.cpp`, tagged `[lte]`. Mock template EoS with c_s = c_b = 1/√3:
- p_s = a_s T⁴/3 − ε, e_s = a_s T⁴ + ε, p_b = a_b T⁴/3, e_b = a_b T⁴, with Ψ_N = a_b/a_s and ε = α_N a_s T_N⁴ (so that α_N = Δθ̄/(3w_s)). Use T ∈ [0.5, 2]·T_N.
- Assertions:
  - `psiN()`/`alN_munu()` reproduce the inputs.
  - `entropy_consistency` ≈ 0.
  - A few (α_N, Ψ_N) points give `Converged`, with vw matching reference values to ~1e-3 and the fit to a few %. Reference values come from running the authors' Python snippet (2303.10171 App. C) in the scratchpad if Python + SciPy are available; otherwise compare to the fit only.
  - |Δ(vw)| < tol at the solution, and `shock_flag()` is true.
  - α_N < (1−Ψ_N)/3 → `NoExpansion`; large α_N → `Runaway`.
  - `with_vw` rescales Rs/Rbar.

**Step 7: example** `example/run_lte_wall_velocity.cpp`: mock EoS → `LTEWallSolver` → `write_scan` → `FluidProfile` / `GWSpec` with `*res.params`. Register it in [example/CMakeLists.txt](example/CMakeLists.txt) and add it to the README run list.

---

## 4. Verification
1. Configure and build: `cmake -S . -B build && cmake --build build -j`, with no new warnings (the `ENABLE_COMPILER_WARNINGS` flags).
2. `./bin/unit_tests` passes in full (regression: public `FluidProfile` behaviour unchanged), and `./bin/unit_tests "[lte]"` passes.
3. `./bin/run_lte_wall_velocity` runs. Inspect `lte_scan.csv`: Δ should be smooth with a single zero crossing, and vw < v_J.
4. Report results, including any tolerance that had to be loosened, faithfully. No commits unless asked.

## 5. References
- Espinosa, Konstandin, No, Servant, arXiv:1004.4187. Junction conditions and modes.
- Giese, Konstandin, van de Vis, arXiv:2004.06995. Template/μν model.
- Konstandin, No, *Hydrodynamic obstruction to bubble expansion*, arXiv:1011.3735.
- Balaji, Spannowsky, Tamarit, *Cosmological bubble friction in local equilibrium*, arXiv:2010.08013.
- Ai, Garbrecht, Tamarit, *Bubble wall velocities in local thermal equilibrium*, arXiv:2109.13710. Entropy matching condition.
- **Ai, Laurent, van de Vis, *Model-independent bubble wall velocities in LTE*, JCAP (2023), arXiv:2303.10171.** Primary reference.
- Ai, Nagels, Vanvlasselaer, *Criterion for ultra-fast bubble walls*, JCAP 03 (2024) 037, arXiv:2401.05911.
- *Bubble-wall velocity in LTE: hydrodynamical simulations vs analytical treatment*, arXiv:2402.15408.
- *Bounds on the bubble wall velocity*, JHEP 02 (2025) 119, arXiv:2411.13641.
- Laurent, Vanvlasselaer, *Dynamical evolution of the pressure on the bubble wall*, arXiv:2606.30740.
- Ekstedt et al., *WallGo*, arXiv:2411.04970. Independent cross-check with an LTE mode.
- Barni, Blasi, Vanvlasselaer, *The hydrodynamics of inverse phase transitions*, arXiv:2406.01596. Future extension.

---

## 6. Implementation notes (feature/vw)

Implemented as planned. The new code is `include/wallvelocity.hpp`, `src/wallvelocity.cpp`, `unit_tests/test_wallvelocity.cpp` and `example/run_lte_wall_velocity.cpp`.

**Validation.** Wall velocities were compared against the Python code in Appendix C of arXiv:2303.10171 (run verbatim, except that `simps` is imported as `scipy.integrate.simpson` for SciPy ≥ 1.14), with c_s² = c_b² = 1/3:

| α_N | Ψ_N | reference v_w | HydroGrav v_w | mode |
|---|---|---|---|---|
| 0.04 | 0.90 | 0.513568 | agrees within 1e-3 (unit test) | deflagration |
| 0.05 | 0.90 | 0.597418 | 0.597390 | hybrid |
| 0.10 | 0.85 | 0.653129 | 0.653126 | hybrid |
| 0.02 | 0.90 | "alN too small" | `NoExpansion` | – |
| 0.30 | 0.90 | "alN too large" | `Runaway` | – |

`jouguet_velocity_veff` reproduces the analytic bag-model v_J.

**Sign convention.** Δ(v_w) = γ₊T₊/(γ₋T₋) − 1 is positive when the driving force wins. The physical root is the + → − crossing below v_J.

**Deviations from the plan.**
- `EquationOfState::smoothing` was added. **Pre-existing issue found:** the smoothing spline in `PTParams_Veff` (λ = 1e-2) biases smooth EoS tables. For a T⁴ table on [0.5, 2] T_N it shifts p(T_N) by about 22% and breaks w = T dp/dT by up to about 200% near the table edges, while c_s² looks fine because the bias cancels in dp/de. The default is unchanged (1e-2), so existing results are unchanged. Set `smoothing = 0` for exact or noise-free tables. `LTEWallSolver` warns when max|w/(T dp/dT) − 1| > 1e-2. **This bias also affects all existing Veff profiles and spectra, so the default deserves a decision.**
- `WallState` is filled for the Veff EoS only; it is NaN for bag/μν.
- `LTEResult::params` is a `shared_ptr` so that `LTEResult` stays copyable.

**Known limitations.**
- The existing profile solver fails for v_w ≲ 0.1. Those scan points are NaN and skipped, so very slow LTE walls are not found.
- Δ(v_w) has jumps of O(1e-5) from the inner T₋ shooting (`find_TmTN_veff`). This limits v_w to roughly 1e-4 precision; observed agreement is 3e-6 to 3e-5.
- The final profile at the solved v_w (n = 5000) may still use the existing μν shock fallback near v_J, as it does for any input v_w.
- Runtime (single core, Release `-O3`): 1.6–3.7 s per `solve()`. A trial profile costs ~0.2 s, almost all in the inner T₋ shooting (`find_TmTN_veff`: a 200-point bracket scan plus golden-section search, ~240 shots). `solve()` caches every trial so no vw is solved twice. It uses an 8-point scan plus Brent (`lte_root_tol` = 1e-5, which matches the residual's precision): 12–19 profiles per convergent solve, and 8 for runaway/no-expansion. Further speed-up would come from making the inner shooting itself faster: a Brent solve on the signed residual, warm starts, adaptive integration.

---

## 7. Planned (not implemented): faster inner shooting in `FluidProfile` (stage 3)

**Status:** proposed only. Stages 1 (Release build) and 2 (cached evaluations, 8-point scan, `lte_root_tol` = 1e-5) are done, and one `LTEWallSolver::solve()` now takes 1.6–3.7 s on one core. All of the remaining cost is in the trial profiles, at about 0.2 s each. These changes would speed up *every* Veff `FluidProfile` (for deflagrations and hybrids), not only the LTE solver, which matters for parameter scans.

### 7.1 Where the time goes (current code, `src/profile.cpp`)
- **`find_TmTN_veff`** finds T₋ by minimising `|v1_residual_veff|`. It first runs `find_bracket`, a 200-point scan over the *whole* EoS table range, then `golden_section_minimize` (about 40 more evaluations). That is about 240 shooting solves per profile.
- **Each shot (`v1_residual_veff` → `deflagration_profile_veff`)** costs about 0.8 ms:
  - `wall_matching_veff`: a `grid_search_2d` (about 0.2 ms) to seed `newton_solve_2d_bounded`, with a Nelder–Mead fallback;
  - a fixed 1000-step RK4 from v₊ᵁᶠ down to v = 1e-10. Each step evaluates the `csq_s`/`csq_b` splines four times, even though the solution is truncated at the shock afterwards;
  - a linear scan for dξ/dv = 0, then `find_shock_idx_veff`, a linear scan over the grid for the minimum of |shock eq. 1|.
- The `profile_n` passed to `FluidProfile` does **not** reach the shooting: `v1_residual_veff` uses its default n = 1000. That is why n = 1000 and n = 5000 profiles cost the same.

### 7.2 Numerical problems the same code causes
- **The shock position is quantised to the RK4 grid.** `find_shock_idx_veff` returns a grid index, not a root. The shooting residual (shock eq. 2) is therefore piecewise constant or jumpy in T₋, with jumps of O(1e-5). This is the most likely origin of the jumps in the LTE Δ(v_w) noted in §6, which limit v_w to about 1e-4.
- **Minimising |f| instead of root-finding on f** halves the attainable precision and needs many more evaluations than Brent's method on a signed residual.
- **Probable cause of frequent μν shock fallbacks.** `check_shock_convergence` requires both shock equations below `minimiser_tol` = 1e-9 at a *grid point*, which quantisation rarely allows. The "Shock residual above tolerance. Using mu-nu fallback" warning then fires and the whole solve is repeated with the approximate condition. The fallback was seen in practice, e.g. for the solved v_w at α_N = 0.05, Ψ_N = 0.9.

### 7.3 Proposed changes, in order
1. **Locate the shock continuously.** Find the sign change of shock eq. 1 between grid points, then refine it by root-finding on a local interpolant or by a short re-integration. That makes the residual a smooth function of T₋. It should remove the O(1e-5) jumps and most μν fallbacks.
2. **Root-find on the signed residual.** Have `v1_residual_veff` return the signed value, and replace `find_bracket` + golden section with a coarse bracket scan (about 10–20 points, or a narrow bracket around a hint) plus `brent_root`. Keep the current minimiser as a fallback for cases with no sign change. This takes about 240 shots down to about 15–25.
3. **Warm starts.**
   - *Within a profile:* seed `newton_solve_2d_bounded` in `wall_matching_veff` from the previous shot's (v₊, T₊), and run `grid_search_2d` only if Newton fails.
   - *Across the LTE scan:* pass a T₋ hint from the neighbouring v_w through the private `FluidProfile` constructor, so the bracket starts narrow.
4. **Adaptive integration with early stop.** Replace the fixed-step RK4 in the shooting with an adaptive integrator such as GSL `odeiv2` (`rk8pd`/`rkf45`). Stop at the first of: the shock condition, dξ/dv = 0, or v → 0. This also allows step 1 to be done properly. The final profile can still be sampled at `n` points for the SSM integrals.
5. **Pass `n` through to the shooting**, or document that the shooting resolution is separate from the output resolution.

**Expected effect (estimates, to be measured):** about 10× per Veff profile from steps 2 and 3, and more with step 4. Together with the LTE solver's 12–19 profiles per solve, that would bring one LTE solve to roughly 0.2–0.5 s on one core. The precision of v_w should improve from about 1e-4 to the root tolerance.

### 7.4 Risks
- This changes the code path behind **every existing Veff result**. Profiles, κ and Ω_GW will move slightly, by about the size of the current quantisation error, and should become more accurate.
- The current minimise-|f| approach copes with residuals that only touch zero or have no clean sign change. The Brent path needs a robust fallback, namely the current method.
- Hybrid profiles near v_J are where the existing solver is most fragile (`vmUF < 0` TODO in `solve_profile_veff`). They need explicit test coverage.

### 7.5 Verification plan
- Before changing anything, record reference outputs from the current code: v(ξ), w(ξ), T(ξ), ξ_sh, wall state, K and Ω_GW. Cover the mock template EoS at several (α_N, Ψ_N, v_w) across deflagration, hybrid and detonation, plus at least one realistic Veff table, e.g. a benchmark from `run_bps`.
- After the change: existing unit tests pass, and the profiles agree with the references to within the old quantisation error. The `shock_flag()` fallback should trigger less often.
- LTE: `[lte]` tests pass, v_w matches the Appendix C reference code to ≤ 1e-5, and Δ(v_w) is smooth in a fine scan around the root.
- Timing: one trial profile and one `solve()` single-core in Release, against the current 0.2 s and 1.6–3.7 s.
