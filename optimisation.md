# Optimising `Spectrum::GWSpec`

Notes on the September 2026 optimisation pass over the sound-shell model GW spectrum
calculation. Two things came out of it: `GWSpec` got roughly **14× faster**, and the
quadrature settings it shipped with turned out to be **6.7% wrong at large `kRs`**.

All timings are the `GWSpec Timer` line from `bin/run_gw_spectrum` (default parameters,
`kRs = logspace(-3, 3, 100)`, `dtau = 10 R*`) on a 12-thread Intel i7-1355U (2 P-cores +
8 E-cores).

| Step | Time | Max error vs over-resolved reference |
|---|---|---|
| Baseline, as found (`build/` was configured `Debug`, i.e. `-O0`) | 2.23 s | 6.7% |
| `-O3` | 1.86 s | 6.7% |
| `schedule(dynamic)` on the k-loop | 1.23 s | 6.7% |
| Analytic `Δ_SSM` kernel | 0.50 s | 6.7% |
| Shared `\|A₊\|²` spline | 0.42 s | 6.7% |
| O(1) kinetic-spectrum interpolant | 0.33 s | 6.7% |
| Resonance-aligned `p̃` quadrature | 0.21 s | **0.07%** |
| `Ap_sq` scheduling | **0.16 s** | 0.07% |

---

## 1. The accuracy problem

The shipped `z_samples = 31`, `z_max_refinements = 5`, `z_tolerance = 1e-6` under-resolved the
inner angular integral at large `kRs`. At `kRs = 1000` the spectrum was **6.7% low**:

| | `Ω(kRs = 1000)` |
|---|---|
| Old scheme, shipped settings | 1.30127e-19 |
| Old scheme, `z_max_refinements = 15`, `z_tolerance = 1e-10` | 1.39528e-19 |
| Independent brute-force integration | 1.395283e-19 |
| New resonance-aligned scheme | 1.39419e-19 |

Three independent routes agree on 1.3953e-19. Low and mid `kRs` were unaffected — the error
grows with `kRs` because the feature being missed shrinks relative to the integration range.

The over-resolved run of the *old* scheme took about ten minutes; the new scheme reaches the
same answer in 0.16 s.

## 2. Why: the integrand has one narrow, analytically locatable resonance

Writing the inner integral over `p̃` instead of the angle `z`, the kinematic range is
`|k − p| ≤ p̃ ≤ k + p`. Of the four sign branches in `Δ_SSM`, **only `cs(p + p̃) = k` has a
solution inside that range** — the `m = −1` branches would need `cs > 1`. So there is exactly
one resonance, at

```
p̃* = k/cs − p,     present for   k(1/cs − 1)/2 ≤ p ≤ k(1/cs + 1)/2
```

(for `cs = 0.567`, that is `0.38k ≤ p ≤ 1.38k`), and around it `Δ_SSM` is a `sinc²` with

```
HWHM = 2.783 / (cs · Δτ)
```

Verified numerically: scanning `Δ_SSM` over `p̃` on a dense grid puts the peak at `k/cs − p`
and gives HWHM = 0.4908 against a predicted 0.4908, at `kRs` = 1, 30 and 300.

With `Δτ = 10 R*` the half-width is ≈ 0.49 in `p̃R*`, while the integration range at
`kRs = 1000` is ~2000 wide. **The adaptive rule was being asked to find a feature 10⁻³ of the
domain width**, and at large `kRs` it stepped over it. Away from the resonance `Δ_SSM` only
decays as `1/(p̃ − p̃*)²`, so the tails cannot simply be truncated either.

## 3. Why: `Δ_SSM` is four evaluations of one 1-D function

The double-time integral reduces, for each sign combination, to

```
G(x) = ¼ |∫_{τs}^{τfin} e^{ixτ} dτ/τ|²
     = ¼ [ (Ci(x·τfin) − Ci(x·τs))² + (Si(x·τfin) − Si(x·τs))² ]
```

so `Δ_SSM = Σ_{m,n} G(cs·p + m·cs·p̃ + n·k)`, with `G` even and fixed for a whole spectrum.
Substituting `τ = τs·t` pulls the common phase out of the modulus:

```
G(x) = ¼ |∫_1^ρ e^{iat} dt/t|²,     a = x·τs,   ρ = τfin/τs
```

Mapping `t` onto `[−1, 1]` about the interval centre and expanding `1/t` in

```
ε = Δτ/(τs + τfin)
```

gives `G(x) = ¼ ε² |Σ_n (−ε)ⁿ Mₙ(w)|²` with `Mₙ(w) = ∫_{−1}^1 sⁿ e^{iws} ds` and
**`w = x·Δτ/2`**. Truncating at `n = 4` leaves an `O(ε⁵)` relative error. The leading term
alone reproduces the familiar

```
G ≈ (Δτ²/4τs²) · sinc²(x·Δτ/2)
```

Two consequences used in the implementation:

- **There is no fast `1/τs` oscillation in `G`** — the modulus kills it. `G` varies only on the
  scale `1/Δτ`. This is what makes the `sinc²` picture in §2 exact rather than approximate.
- `G` depends on `x` only through `w`, so the four `(m, n)` branches share sine/cosine
  evaluations by angle addition: `w = (cs(p ± p̃) ± k)·Δτ/2`, and `k·Δτ/2` is fixed per `k`.

Result: **eight ALGLIB `sinecosineintegrals` calls per quadrature node became two `sin`/`cos`
calls.** Measured agreement with the original implementation, over 4×10⁵ random `(k, p, z)`:

| `Δτ/τs` | ε | worst \|Δ − Δ_ref\| relative to peak `Δ` |
|---|---|---|
| 0.011 (default) | 0.0055 | 7.7e-14 |
| 0.05 | 0.024 | 6.0e-11 |
| 0.1 | 0.048 | 3.3e-09 |
| 0.3 | 0.130 | 1.4e-06 |

Above `ε = 0.15` (`Δτ ≳ 0.35 τs`, i.e. sound waves lasting an appreciable fraction of a Hubble
time) the class falls back to exact Si/Ci automatically. For the default parameters
`τs = 622.5 R*` and `ε = 0.008`, so the fast path is used with a wide margin.

## 4. What was changed

### `CMakeLists.txt`
- Defaults `CMAKE_BUILD_TYPE` to `Release` when the caller has not set one. The repo's
  `build/` was configured `Debug`, so everything was compiled `-O0`.
- Adds opt-in `-DHYDROGRAV_NATIVE_ARCH=ON` (`-march=native -ffp-contract=fast`).
  `-ffast-math` is deliberately **not** used: the code relies on defined NaN/Inf behaviour.

> Worth noting this mattered much less than expected — 2.23 s → 1.86 s, not the usual 5–20×,
> because the hot code was inside the precompiled `libalglib`, not in this project's objects.
> The real wins were algorithmic.

### `include/ssm_kernel.hpp` (new)
`SoundShellKernel` — evaluates `G(x)` per §3. Constructed once per spectrum from
`(τs, τfin)`; exposes `from_w(w, sin w, cos w)` for the angle-addition path and
`operator()(x)` for general use, plus `uses_exact()` for the fallback.

### `src/ssm.cpp`
- `dlt_SSM` reimplemented over `SoundShellKernel`. The original six-argument signature is kept
  and delegates; a hot-loop overload takes a prebuilt kernel and the per-`k` `sin`/`cos`.
- `GWSpec`'s inner integral changed from `z` to `p̃`
  (`dz = −p̃ dp̃/(kp)`), split into fixed-order Gauss–Legendre panels whose widths grow
  geometrically away from `p̃*` (`resonance_panels()`). The peak is resolved regardless of how
  small it is relative to the range, and the `1/x²` tails cost only logarithmically many
  panels. No adaptivity, so the outer integrand is smooth and deterministic.
- The outer `p` range is split at `k(1/cs ∓ 1)/2` — the band where the resonance exists —
  rather than leaving the adaptive rule to discover it.
- `Hydrodynamics::ApsqSpline` added: `|A₊|²(χ)` was being rebuilt 2–3 times per `GWSpec` call
  (via `get_nl_timescale`, `build_kinetic_spectrum_spline` and `gw_prefac`). It is now built
  once and passed to each; new overloads of `Ekin`, `zetaKin`, `build_kinetic_spectrum_spline`,
  `get_nl_timescale` and `gw_prefac` accept it. Results are unchanged.
- `kinetic_spectrum_interpolant()` replaces the `exp(spline1dcalc(log_zk_spline, ·))` pattern
  in the innermost loop with a `LogGridInterpolant` storing ζ directly — no ALGLIB state, no
  binary search, no `exp`.
- `schedule(dynamic, 1)` on the k-loop, `schedule(dynamic, 16)` on the χ loop of
  `fluid_profile_integrals`. Both had `schedule(static)` over work that varies by orders of
  magnitude across the loop range; on a hybrid P/E-core CPU this was leaving most threads idle.
  The χ fix alone took `Ap_sq` from 0.080 s to 0.030 s.

### `include/maths.hpp` / `src/maths.cpp`
- `LogGridInterpolant` — O(1) Catmull-Rom interpolation on a `logspace()` grid. Index
  arithmetic instead of a binary search; `at_log()` for callers that already hold `ln x`.
  Third-order accurate; ~3e-6 on the production grid.
- `FilonQuadrature::filon_integrate_interval` replaced by `filon_subdivided`, which walks the
  nodes in a rolling window: one `f` evaluation and one `sin`/`cos` pair per node, and no
  per-subinterval `std::vector` allocations. Bit-identical results.

### `include/ssm.hpp`
New declarations for the above; `PowerSpec::profile()` now returns `const&` instead of
returning a 5000-point profile **by value** on every call.

### `include/config.hpp`
`z_samples` / `z_max_refinements` / `z_tolerance` are gone, replaced by the panel parameters:

```cpp
pt_gauss_legendre_samples = 15   // GL points per panel
pt_panel_inner_width      = 2.0  // innermost panel width, in resonance half-widths
pt_panel_growth           = 3.0  // geometric growth away from the resonance
pt_max_panels             = 32   // cap per side
sinc_sq_hwhm              = 1.39155737825151
```

`pRs_max_refinements` 5 → 4 and `pRs_tolerance` 1e-6 → 1e-5. With the inner integral
deterministic and the `p` range pre-split, the outer rule converges quickly and the residual
error of the spectrum is set by the panel resolution, not by these.

These were chosen by sweeping against the over-resolved reference:

| Setting | k-loop time | Max error |
|---|---|---|
| GL=15, growth 2, width 1, refine 5, tol 1e-6 | 0.220 s | 0.123% |
| **GL=15, growth 3, width 2, refine 4, tol 1e-5** | **0.069 s** | **0.069%** |
| GL=15, growth 3, width 2, refine 3, tol 1e-4 | 0.038 s | 0.066% |
| GL=15, growth 4, width 2, refine 3, tol 1e-4 | 0.037 s | 0.658% |
| GL=10, growth 3, width 2, refine 3, tol 1e-4 | 0.025 s | 0.868% |

`refine = 4` / `tol = 1e-5` was taken over the marginally faster `refine = 3` for margin in
other parameter regimes; the panel resolution dominates the error either way.

## 5. Validation

Cross-checked against a dense independent evaluation of the same integral (uniform `p̃` grid
at 8 points per resonance half-width, 6000 log-`p` points):

| Case | mode | `Δτ/R*` | ε | max error |
|---|---|---|---|---|
| deflagration, `dtau = 10 R*` | deflagration | 10.0 | 0.0080 | 0.062% |
| hybrid, `dtau = 10 R*` | hybrid | 10.0 | 0.0080 | 0.052% |
| detonation, `dtau = 10 R*` | detonation | 10.0 | 0.0080 | 0.063% |
| bag EoS, `dtau = 10 R*` | hybrid | 10.0 | 0.0080 | 0.059% |
| hybrid, `dtau = 2 R*` | hybrid | 2.0 | 0.0016 | 0.053% |
| hybrid, `dtau = 100 R*` | hybrid | 100.0 | 0.0744 | 0.053% |
| hybrid, `dtau` from nl timescale | hybrid | 10.6 | 0.0085 | 0.047% |
| detonation, `dtau` from nl timescale | detonation | 12.1 | 0.0097 | 0.065% |

Whole-spectrum deviation from the over-resolved reference: **6.9e-4**, worst at
`kRs = 35.1`. LISA SNR is unchanged at 0.271878.

New tests (38 test cases, 291,989 assertions, all passing):
- `[ssm_kernel]` — `dlt_SSM` against a literal Si/Ci reference across `(k, p, z)` and several
  `Δτ/τs`; the exact-fallback path; evenness of `G`.
- `[loggrid]` — knot reproduction, interpolation accuracy, third-order convergence under grid
  refinement, out-of-range behaviour, malformed input.
- `[gwspec]` — spectrum pinned at five `kRs` against the brute-force values, to 2e-3.

## 6. Known gotcha (pre-existing, not changed)

`gw_prefac()` normalises by the peak of `Ekin()` **over the momentum grid it is handed**, so
the absolute amplitude of `GWSpec` depends on the `kRs` grid you pass. A coarse grid samples
the `Ekin` peak poorly and shifts the whole spectrum. This bit the regression test — five
`kRs` values give visibly different amplitudes from the same five points taken out of a
100-point grid.

This was left alone deliberately: fixing it changes published numbers and is a physics
decision, not a performance one. But it is worth resolving before comparing spectra computed
on different grids.

## 7. Reproducing

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release   # now the default
cmake --build build -j12
./bin/unit_tests
./bin/run_gw_spectrum
```

To regenerate the over-resolved reference for the old scheme, revert `src/ssm.cpp` and set
`z_max_refinements = 15`, `z_tolerance = 1e-10`, `pRs_max_refinements = 12`,
`pRs_tolerance = 1e-9` — expect ~10 minutes.

## 8. Not done

- **Tabulating `G`.** `G` could be precomputed once per `GWSpec` call on a grid in
  `y = x·Δτ`, reducing the innermost cost to four interpolations. Not needed: after the panel
  change the k-loop is ~0.04 s of the 0.16 s total, and the setup (fluid profile, `|A₊|²`,
  the `Ekin` evaluations) now dominates.
- **The stationary / δ-function limit**, which would collapse the double integral to a single
  one. It is an approximation with `O(Δτ/τs)` error and the exact path is now fast enough.
- Further work on `GWSpec` should target the ~0.10 s of setup, not the k-loop.
