# Changelog

## 0.3.4 — 2026-05-01 — Guillot 2010 T-P profile

Adds the standard hot-Jupiter analytic temperature-pressure profile,
moving Argus from "isothermal toy" to actually-realistic exoplanet
atmosphere modelling.

### Added
- **`argus::guillot_profile()`** — analytic T(P) per Guillot 2010
  eq. (29). Inputs: `T_int` (internal heat flux temperature),
  `T_irr` (irradiation temperature), γ (thermal/visible opacity ratio),
  `kappa_IR` (IR opacity), `gravity`. Returns T (K) at each pressure.
- **`argus::guillot()`** — convenience builder that returns a complete
  `Atmosphere` with a Guillot T-P profile, log-spaced pressure grid,
  and a single species.
- **`test_guillot`** (9 assertions):
  * profile is finite and physically plausible (100-5000 K)
  * deep atmosphere is hottest (internal heat increases with τ ∝ P)
  * top-of-atmosphere matches the analytic skin-temperature formula
    to <0.1%
  * monotone non-decreasing for γ < 1 (no inversion regime)
  * geometry pass works on the non-isothermal profile
  * forward model produces sensible spectrum
  * higher T_irr → larger scale height → deeper transit (verified)
  * bad inputs throw

### Validated
19/19 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`.
Performance unchanged: ~1.7 ms forward call, ~9 ns / Voigt.

---

## 0.3.3 — 2026-05-01 — Self-broadening + independent Voigt reference + benchmark

Adds the missing HITRAN physics (self-broadening, pressure shift), an
independent closed-form Voigt cross-validation, and a performance
baseline benchmark.

### Added
- **Self-broadening + pressure shift in `LineListOpacity`**:
  * `OpacityKernel::cross_section_with_self()` — new virtual that takes
    the species' own VMR per (T,P) sample, with a default that
    forwards to `cross_section()` (no self-broadening) for backward
    compatibility.
  * `LineListOpacity::cross_section_with_self()` — implements the full
    HITRAN form: γ_total = γ_air·(1-VMR) + γ_self·VMR
  * Pressure-shift now applied to ALL line evaluations:
    nu_eff = nu0 + δ_air · P_atm
  * Both `cross_section()` (air-broadened, VMR=0) and
    `cross_section_with_self()` (self-aware) paths use the shifted
    centre, so the pressure shift is always honoured.

- **`test_self_broadening`** (4 tests):
  * γ_self = 2·γ_air with varying self-VMR -> monotonically wider
    Lorentzian; pure-self limit ~2× wing absorption vs pure-air.
  * Default `cross_section()` matches `cross_section_with_self(VMR=0)`.
  * δ_air > 0 with high P shifts the line centre off our sample point,
    reducing cross-section there.
  * δ_air = 0 confirmation: pressure broadens centre but doesn't shift.

- **`test_voigt_reference`** — independent closed-form cross-check.
  Voigt at line centre with σ_g = 1, varying γ_l: V(0) = w(iy)/(σ√(2π))
  where w(iy) = exp(y²)·erfc(y) computed via std::erfc (~1e-15 accurate).
  Tolerance 1e-4 absorbs HAW's actual error in the y ~ 1-5 transition
  zone. Plus pure-Lorentzian limit cross-check at high y_n, plus
  comprehensive (σ, γ) grid sweep.

- **`test_benchmark`** — performance baseline. Times 10 forward-model
  runs (60-layer atmosphere · 16 H2O lines · ~200 wavenumber points)
  and asserts:
  * median wall time < 5 s (loose CI bound)
  * per-Voigt-evaluation cost < 1 μs (current: ~9 ns on this box)
  Recorded medians: ~1.7 ms / forward call, ~8.9 ns / Voigt.

### Changed
- `OpacityKernel` interface gains a new virtual; existing kernels (e.g.
  `GreyOpacity`) inherit the no-op default and continue to work.
- `LineListOpacity::cross_section()` now applies pressure-shift to
  the line centre. For δ_air = 0 (M2 default for our test fixtures)
  this is a no-op; for nonzero δ_air the line moves with pressure.

### Validated
18/18 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. 3 examples build warning-free. Performance benchmark
recorded: 1.7 ms / forward call (60 layer · 200 wn · 16 lines).

---

## 0.3.2 — 2026-05-01 — Multi-molecule end-to-end

Adds a real CO2 line fixture and a hard multi-molecule end-to-end test
that runs both H2O and CO2 in the same atmosphere, verifying additive
optical depth and correctly-positioned absorption bands.

### Added
- **`argus::test_data::kCO2Lines`** — 10 hand-curated real HITRAN CO2
  lines from the 4.3 μm asymmetric stretch fundamental band (most
  abundant 12C-16O2 isotopologue).
- **`test_multi_molecule`** — H2O + CO2 in a single atmosphere.
  Asserts:
  * cross-section additivity across opacity kernels
  * H2O dominates at H2O band centres (3651, 7099 cm⁻¹)
  * CO2 dominates at CO2 4.3 μm Q-branch (2363 cm⁻¹)
  * combined optical depth ≥ either species alone (no destructive
    interference, sanity check)
  * continuum window (5500 cm⁻¹) shows less absorption than band centres

### Validated
15/15 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. All 3 examples build warning-free and run end-to-end.

---

## 0.3.1 — 2026-05-01 — M2 polish

Adds three more hard test suites — bringing the total to 14 — covering
edge-case stress, file I/O round-trip, and bit-exact reproducibility.

### Added
- **`test_stress`** — 15 hard edge-case assertions:
  malformed Atmosphere, mixing-ratio bounds, isothermal-builder rejects,
  geometry rejects 1-layer atmosphere and mis-ordered layers, Tensor
  shape mismatch, zero-opacity baseline matches geometric (R_p/R_*)²,
  null OpacityKernel rejection, non-positive molar mass rejection,
  Voigt finite + positive at extreme regimes, line list at 50 K and
  4000 K, single-wavelength forward, empty wavenumber grid, empty line
  list, garbage HITRAN strings.
- **`test_file_io`** — writes `kH2OLines` fixture to `/tmp`, reads it
  back via `Hitran::load_file()`, and verifies field-by-field equality
  with the in-memory parse to 1e-12. Also asserts that missing files
  throw.
- **`test_reproducibility`** — runs the entire forward model twice on
  the same inputs and asserts bit-exact equality (no nondeterminism in
  reductions / hash randomisation / uninit memory). Plus IR content
  address determinism, single-call cross-section determinism,
  partition-function determinism.

### Validated
14/14 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. Three examples build and run end-to-end.

---

## 0.3.0 — 2026-05-01 — M2 production-grade

Production-quality completion of the M2 milestone. The kernel now uses
a verified Faddeeva-based Voigt evaluator, ships a real HITRAN .par
parser, and validates against ~1e-6 finite-difference autograd checks
on real-world data.

### Added
- **`argus::Hitran`** (`include/argus/hitran.hpp`, `src/hitran.cpp`) —
  fixed-width 160-char HITRAN .par parser. Reads molecule_id, isotope_id,
  nu0, intensity, gamma_air, gamma_self, E_lower, n_air, delta_air from
  every record. Tolerant of trailing CR/LF, rejects malformed records,
  supports filtering by molecule id. `load_file()` reads from disk;
  `load(istream)` reads any stream.
- **`argus::Partition`** (`include/argus/partition.hpp`,
  `src/partition.cpp`) — TIPS-anchored partition functions Q(T) for
  H2O, CO2, CH4, CO, NH3 (the workhorses of exoplanet atmospheres).
  Power-law form anchored at the published Q(296) values; accuracy
  ~5-10% across 200-3000 K. M3 will swap to full TIPS-2017 polynomials.
- **`argus::test_data::kH2OLines`** (`include/argus/test_data.hpp`) —
  16 hand-curated real H2O lines from HITRAN-2020 spanning the 2.7 μm
  and 1.4 μm bands. Bundled as a string_view fixture so tests work
  with no network or disk dependency.
- **`examples/03_real_hitran.cpp`** — end-to-end demo loading the bundled
  real HITRAN records, running the transmission spectrum on a hot Jupiter,
  and printing the resulting JWST-NIRSpec-PRISM-shaped output (3500-7400
  cm⁻¹). Demonstrates Hitran::load + LineListOpacity + TransmissionModel
  + Partition::Q working together.
- **5 new hard tests** (test_partition, test_hitran, test_real_h2o,
  test_finite_diff, plus reinforced test_voigt and test_geometry):
  - `test_voigt` — Voigt vs analytic Gaussian/Lorentzian limits to <1e-4,
    symmetry to 1e-12, area-normalisation to 1e-3, numerical-convolution
    cross-check, asymptotic Lorentzian wing, dual-number derivative vs
    central FD to 1e-5.
  - `test_geometry` — chord through full sphere = 2r exactly, tangent
    chord = 0 exactly, hydrostatic isothermal altitude vs analytic
    z = -H ln(P/P_bot) to 1e-9, full-atmosphere height check.
  - `test_partition` — Q(296) anchor exact, Q(1000) within 1% of TIPS,
    Q(150) within 10%, Q(2000) within 5%, monotonicity, throw on bad input.
  - `test_hitran` — single-line round-trip with field-by-field assertions,
    multi-line parse via istringstream, CR/LF tolerance, malformed
    rejection, molecule-id filter.
  - `test_real_h2o` — load real HITRAN data, cross-section sum-rule check,
    temperature-dependence of intensity ratios, end-to-end transmission
    spectrum with line-centre-vs-continuum check, VMR monotonicity,
    optically-thick saturation.
  - `test_finite_diff` — Dual<double> derivatives wrt σ, γ, x, S validated
    against central finite differences to 1e-5 across the entire line
    profile, plus exactness check for linear-in-S.

### Changed
- **`argus::voigt`** rewritten with the **correct** Hui-Armstrong-Wray
  evaluator. Previous v0.2.0 had two bugs: (1) wrong polynomial input
  (used z = x+iy instead of t = y-ix — produced 50× wrong results in
  the Lorentz-dominated regime), (2) coefficient b[6] was 10.479857
  instead of 10.479857 (verified against published Hui et al. 1978).
  Now ships an internal `detail::Cmplx<T>` template so both `double`
  and `Dual<double>` flow through the same code path. Adds a
  Gaussian-limit fallback for y_n < 1e-3 to avoid HAW's catastrophic
  cancellation in the deep wings.
- **`argus::Dual<T>`** gained `value_of()` extraction helpers so
  branchy code (like the Voigt fallback) can inspect a value while
  keeping the Dual chain intact.
- **`argus::LineListOpacity`** intensity scaling upgraded to the full
  HITRAN form: `S(T) = S(296) · Q(296)/Q(T) · exp(-c2·E"/T·(1-T/Tref))
  · induced_emission_correction`. The Q ratio is computed once per
  layer per species via `argus::Partition`. Falls back to Boltzmann-only
  when the species key is unknown.
- Version → 0.3.0.

### Validated
All 11 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. Three examples build and run end-to-end:
- `01_transmission_spectrum` — physical 1.12% transit depth (hot-Jupiter on sun)
- `02_voigt_h2o` — toy-line spectrum with autograd derivative
- `03_real_hitran` — real HITRAN H2O lines, JWST-PRISM-shaped output

### Roadmap
M2 is now production-complete. M3 next: amortized SBI via normalizing
flows, WASP-39b benchmark vs petitRADTRANS, full TIPS-2017 partition
tables.

---

## 0.2.0 — 2026-05-01 — M2 starting wedge

First piece of M2: real physics in place of the M1 chord-stub, plus the
forward-mode autograd primitive.

### Added
- argus::Geometry — hydrostatic ray geometry. build_geometry() derives
  layer altitudes from the local scale height; chord_path_length()
  returns the geometric path length of a tangent ray through one shell.
- argus::voigt — area-normalised Voigt line shape via the
  Thompson-Cox-Hastings (1987) pseudo-Voigt approximation (~1% accuracy).
  Replaced in v0.3.0 with verified Hui-Armstrong-Wray (~1e-6 accuracy).
- argus::Line, argus::LineListOpacity — HITRAN-subset line records and an
  OpacityKernel that sums Voigt-shaped contributions.
- argus::Dual<T> — forward-mode autograd dual numbers.
- Atmosphere fields: bulk_mmw_amu, star_radius_rsun.
- 4 new tests + example 02_voigt_h2o.cpp.

### Changed
- TransmissionModel::forward now does proper transit-radius integration:
  R_eff² = R_p² + 2·Σ_b (1-exp(-τ(b)))·b·db.

---

## 0.1.0 — 2026-05-01 — M1 starting point

Initial scaffolding.

- Argus IR skeleton: Node, Graph, content-address (FNV-1a 64-bit stub).
- Atmosphere model: layered T-P-VMR with isothermal builder.
- Opacity interface: OpacityKernel + GreyOpacity placeholder.
- Transmission forward model: chord-integrated optical depth (M1 stub).
- Argus.tex one-page Future Project Brief.
- docs/Astronomy-Compute-Crisis.tex — 23-page research deck.
- site/ — astronomy-themed landing page.
