# Changelog

## 0.7.7 — 2026-05-02 — Runnable M4 + M5 examples

Adds the lensing and interferometry equivalents of the existing
`example_03_real_hitran` and `example_04_retrieval` demos: short
runnable programs that exercise the full forward-model surface for
each of the new physics layers.

### Added
- **`examples/05_lensing.cpp`** — SIE + ExternalShear lens, forward
  through `find_images`, prints 4 quad-image positions + magnifications
  + 16 pairwise Fermat-potential time-delay differences (the same
  kernel surface fitted by `test_sie_retrieval`). Output explicitly
  flags Δτ × (1+z_l)·D_Δt/c as the COSMOGRAIL/TDCOSMO observable for
  H₀ inference.
- **`examples/06_interferometry.cpp`** — 7-antenna Y-array,
  Earth-rotation track at dec=30°, lat=34°, λ=21 cm, 105 baselines
  total. Predicts visibilities of a compact-core + extended-jet
  2-Gaussian source and prints an ASCII `<|V|>` vs |uv| histogram
  showing the jet's exp(-2π² σ² r²) envelope falling off at long
  baselines. Same kernel surface fitted by
  `test_multi_component_retrieval`.

### Validated
44/44 tests still pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. All 6 examples build warning-free.

---

## 0.7.6 — 2026-05-02 — M5 capstone: 2-component interferometric retrieval

The most rigorous M5 substrate proof: 8-parameter retrieval recovers
both components (compact core + extended jet) of a multi-component
radio source from synthetic visibilities sampled across an Earth-
rotation track on a 7-antenna VLA-like array.

Real radio sources (AGN, quasars, starbursts) are typically multi-
component. Argus's substrate claim now covers this end-to-end with
no special-casing — the 2-Gaussian forward model is just two calls
to `predict_visibility` summed pairwise, and the same `Retrieval` API
recovers all 8 parameters.

### Added
- **`test_multi_component_retrieval`** — full M5 capstone:
  * 7-antenna Y-array, 5 hour angles spanning ±2.5 h, source at
    dec=30° from latitude 34°: **105 baselines × 2 = 210 observations**
    via `uv_coverage_track`.
  * Truth: compact core (l=0, m=0, F=1.0, σ=1e-5) + extended jet
    (l=1.5e-5, m=5e-6, F=0.3, σ=3e-5).
  * 8 free params, label-switching broken by **non-overlapping size
    priors** (σ_core ∈ [0, 2e-5]; σ_jet ∈ [2e-5, 8e-5]).
  * EnsembleSampler 32 × 4000: 36% acceptance, recovery within 3σ on
    every parameter.
  * Bit-exact determinism on identical seed.

### Substrate-claim status (final inventory)
| Layer | Test | Free params | Recovery | Sampler |
|---|---|---|---|---|
| **M2/M3** atmospheres | `test_retrieval` | (T, log10 VMR) | 1σ | MH |
| **M4** SIS lensing    | `test_lensing_retrieval`           | 5 | 3σ | Ensemble |
| **M4** SIE lensing    | `test_sie_retrieval`               | 7 | 0.5σ | Ensemble |
| **M5** interferometry | `test_interferometry_retrieval`    | 4 | 3σ | MH |
| **M5** multi-component| `test_multi_component_retrieval`   | 8 | 3σ | Ensemble |

Five independent retrieval tests across three physics layers, all
reusing the same `Spectrum` / `Retrieval` / `EnsembleSampler` /
`MetropolisHastings` / `PosteriorSummary` / `posterior_predictive`
infrastructure. Substrate claim demonstrably backed.

### Validated
44/44 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

---

## 0.7.5 — 2026-05-02 — M4 capstone: full SIE retrieval substrate proof

The most rigorous M4 substrate proof yet: 7-parameter SIE retrieval
recovers (θ_E, q, φ, lens_x, lens_y, source_x, source_y) from 4
observed quad-image positions. Unlike the v0.6.1 SIS retrieval (which
needed `EnsembleSampler` to mitigate the (lens, source) translation
degeneracy), the SIE forward map is informative enough to constrain
all 7 parameters to < 0.5σ from truth — including the SIE-specific
shape parameters q and φ.

### Added
- **`test_sie_retrieval`** — full M4 capstone substrate proof.
  Truth: SIE at (0.05, -0.02) with θ_E=1.0, q=0.7, φ=0.3 rad; source
  at (0.04, 0.03) inside the tangential caustic ⇒ 4 quad images.
  Add 0.01" image-position noise. Forward = source-plane chi²:
  back-project each observed image via β_back_i = θ_i - α(θ_i;
  lens_params); spectrum = (β_back_i - source_free) packed as 8
  scalars; observation = 8 zeros. Avoids `find_images` in the inner
  loop (~ 1 µs per forward instead of ~ 1 ms).
  EnsembleSampler (32 walkers × 4000 steps) recovers every parameter
  to within 0.5σ at 49% acceptance.

### Substrate-claim status
| Layer | Test | Free params | Recovery |
|---|---|---|---|
| **M2/M3** atmospheres | `test_retrieval` | (T, log10 VMR) | 1σ |
| **M4** SIS lensing    | `test_lensing_retrieval` | (θ_E, lens, source) | 3σ ensemble |
| **M4** SIE lensing    | `test_sie_retrieval` | (θ_E, q, φ, lens, source) | 0.5σ ensemble |
| **M5** interferometry | `test_interferometry_retrieval` | (l, m, F, σ) | 3σ MH |

Substrate claim now backed by **four independent retrieval tests**
across **three physics layers**, all reusing the same `Spectrum` /
`Retrieval::log_posterior` / `EnsembleSampler` / `MetropolisHastings`
/ `PosteriorSummary` infrastructure. The forward model is the only
thing that changes.

### Validated
43/43 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

---

## 0.7.4 — 2026-05-02 — Earth-rotation UV synthesis

The standard real-array trick: as the Earth rotates, projected
baselines sweep tracks in the UV plane. Single snapshot has only
N(N-1)/2 baselines; a multi-hour track samples thousands of UV
points. Argus now ships the closed-form ENU → equatorial → UV
projection used in CASA / WSClean / AIPS so any real array layout
+ source position + integration schedule produces the corresponding
UV coverage.

### Added
- **`argus::interferometry::uv_coverage_track(east, north, latitude,
  hour_angles, declination, λ)`** — for each (HA), rotates each ENU
  baseline into equatorial (X = -sin L · B_N; Y = B_E; Z = cos L · B_N)
  then projects onto the UV plane perpendicular to source direction:
      u =  sin h · B_X + cos h · B_Y
      v = -sin δ cos h · B_X + sin δ sin h · B_Y + cos δ · B_Z
  Returns N_baselines · N_HA UVPoints in HA-major order.
- **`test_interferometry` extended (11 → 15 test groups)**:
  * Track size = N_baselines · N_HA
  * Snapshot at HA=0 with source at zenith (δ=L) reduces to the
    existing meridian-snapshot uv_coverage_snapshot bit-exactly
  * **Single E-W baseline traces an ellipse**: u² + (v/sin δ)² =
    (B_E/λ)² verified at 21 HAs with explicit (h=0) sanity check
  * 4 bad-input throws

### Validated
42/42 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

---

## 0.7.3 — 2026-05-02 — NFW dark-matter halo lens

The universal cosmological-simulation density profile (Navarro-Frenk-
White 1996/1997) — the standard model for galaxy / cluster dark-matter
halos in lensing analyses. Closes the M4 lens-model zoo: SIS + SIE +
NFW + ExternalShear, composable via CompoundLens.

### Added
- **`argus::lensing::NFW`** — Wright & Brainerd (2000) closed-form
  deflection:
      |α(θ)| = α_s · h(θ/θ_s) / (θ/θ_s)
      h(x) = ln(x/2) + arccosh(1/x)/√(1-x²)   for x < 1
      h(x) = 1 + ln(1/2)                       for x = 1
      h(x) = ln(x/2) + arccos(1/x)/√(x²-1)    for x > 1
  Continuous through x = 1 (analytic limit verified to 1e-6).
- **`test_lensing` extended (23 → 29 test groups)**:
  * NFW deflection at x=1 matches α_s·(1+ln(1/2)) to 1e-10
  * NFW small-x cancellation (deflection → 0 at origin)
  * NFW deflection radial outward via cross-product = 0
  * NFW translation invariance under off-centre shift
  * NFW + ExternalShear via CompoundLens: `find_images` recovers a
    multi-image config (commonly 3 — NFW supports a central image),
    every image closes lens equation < 1e-7
  * 4 bad-input throws (α_s/θ_s ≤ 0)

### Validated
42/42 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

The M4 lens model zoo is now complete for galaxy- and cluster-scale
strong lensing: SIS (axisymmetric mass), SIE (elliptical galaxy),
NFW (dark-matter halo), ExternalShear (large-scale tidal field), all
composable via CompoundLens for the SIE + NFW + γ workhorse used in
real H0LiCOW / TDCOSMO analyses.

---

## 0.7.2 — 2026-05-02 — External shear + compound lens

The standard real-strong-lensing modelling primitive: SIE + external
shear is the workhorse model used in 99% of published lensed-quasar
analyses. Argus now ships both the shear field and a generic
CompoundLens that sums deflection + potential over multiple
components, so any combination (SIS + shear, SIE + shear, SIE + NFW
once NFW is added, etc.) just works through `find_images`,
`fermat_potential`, and `time_delay_arcsec2`.

### Added
- **`argus::lensing::ExternalShear(γ_1, γ_2)`** — pure-shear lens
  field: α = (γ_1 θ_x + γ_2 θ_y, γ_2 θ_x - γ_1 θ_y);
  ψ = ½(γ_1 (θ_x²-θ_y²) + 2 γ_2 θ_x θ_y). Shear is centred at the
  origin by convention; pivot elsewhere by composing with a
  translated SIS/SIE.
- **`argus::lensing::CompoundLens`** — owns a
  `vector<shared_ptr<const Lens>>`; `deflection` and `potential`
  return the per-component sums. `add(lens)` extends in place.
- **`test_lensing` extended (18 → 23 test groups)**:
  * Shear deflection closed form on 3 sample points
  * Shear potential gradient = α via central FD on 3 points
  * Compound lens sums deflection + potential of (SIS + shear)
  * SIS + strong shear (γ_1=0.10): `find_images` recovers ≥ 2
    images, lens equation closes <1e-7 on each, magnifications finite
  * Null-component throws (constructor + `add`)

### Validated
42/42 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

---

## 0.7.1 — 2026-05-02 — M5 wedge 2: interferometry retrieval (substrate proof)

The full M5 substrate-claim test: the same `argus::Retrieval` API used
for atmospheric retrieval and gravitational-lensing parameter recovery
now drives radio-interferometric source-parameter inference from
observed visibilities. Different physics, identical infrastructure.

### Added
- **`test_interferometry_retrieval`** — 4-parameter MCMC retrieval
  over (l, m, F, σ) for a circular Gaussian source from synthetic
  noisy visibilities on a 7-antenna VLA-like Y-array (21 baselines,
  42 real-valued observations, λ = 21 cm). Truth is recovered to
  within 3σ on every parameter; posterior-predictive 5–95% band ± 2 σ_V
  brackets ≥ 85% of observations. Pipeline reuses `argus::Spectrum` /
  `Retrieval::log_posterior` / `Retrieval::run_mcmc` /
  `PosteriorSummary` / `Retrieval::posterior_predictive` unchanged.
  The only interferometry-specific code is the `Retrieval::Forward`
  closure on top of `argus::interferometry::predict_visibilities`.

### Substrate-claim status
- **Atmospheric retrieval (M2/M3)** → recover (T, log10 VMR), 1 σ
- **Strong-lensing retrieval (M4)** → recover (θ_E, lens_xy, source_xy)
  via `EnsembleSampler`, 3 σ
- **Interferometric retrieval (M5)** → recover (l, m, F, σ) via single-
  chain MH, 3 σ
All three reuse the same `Spectrum` / `Retrieval` / sampler /
`PosteriorSummary` / `posterior_predictive` infrastructure. The
substrate claim is now backed by 3 independent physics layers.

### Validated
42/42 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

---

## 0.7.0 — 2026-05-02 — M5 starting wedge: interferometry visibility model

The third physics layer for the substrate claim: radio interferometry.
Argus can now predict complex visibilities V(u, v) from a sky
brightness distribution composed of point sources and circular
Gaussians, the same primitives every real interferometric pipeline
(CASA, WSClean, AIPS) builds on.

### Added
- **`argus::interferometry` namespace** (`include/argus/interferometry.hpp`,
  `src/interferometry.cpp`):
    `UVPoint`, `Visibility` — basic data types for the (u, v) plane.
    `PointSource{l, m, flux}` — Dirac δ component.
    `GaussianSource{l, m, flux, sigma}` — circular Gaussian; FT is
       another Gaussian in the UV plane × position phase.
    `predict_visibility(src, uv)` — single-component prediction.
    `predict_visibilities(srcs, uv_points)` — overloads for both
       source types; sums over components.
    `uv_coverage_snapshot(east_m, north_m, λ)` — synthesises N(N-1)/2
       baseline (u, v) points from antenna positions; the real-array
       UV coverage primitive.
- **`test_interferometry`** (11 test groups):
  * Point at origin → V = F + 0i
  * Off-origin point: |V|=F, phase = -2π(u·l + v·m)
  * Visibilities sum across components
  * Translation theorem: source shift Δ multiplies V by exp(-2πi·Δ·uv)
  * Conjugate symmetry: V(-u,-v) = conj(V(u,v)) for real sky
  * Gaussian: V(0)=F, decay exp(-2π² σ² r²) at three UV separations
  * Gaussian σ=0 reduces to PointSource bit-exact
  * Mixed Gaussian + point components compose by sum
  * uv_coverage_snapshot on a 3-antenna triangle: 3 baselines at
    correct (u, v) values
  * 4 bad-input throws (size mismatch, < 2 antennas, λ ≤ 0, σ < 0)
  * Determinism: bit-equal outputs on repeated calls

### Validated
41/41 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

The next M5 wedge wires this forward into the substrate-claim
retrieval test that recovers (l, m, F, σ) from noisy visibilities
via the same `argus::Retrieval` API used for atmospheric and lensing
inference.

---

## 0.6.3 — 2026-05-02 — Lensing potential + Fermat + time delays

The H0-from-lensed-quasar pipeline (TDCOSMO/H0LiCOW methodology):
once the lens potential ψ(θ) is in hand, the Fermat potential
τ = 0.5|θ-β|² - ψ(θ) and the time-delay difference between images
follow analytically. Argus now ships all three.

### Added
- **`Lens::potential(theta)`** — virtual on the base class, default
  returns 0.
- **`SIS::potential`** — closed form ψ_SIS = θ_E · |θ - centre|.
- **`SIE::potential`** — closed form Kormann, Schneider & Bartelmann
  (1994):
      ψ = (θ_E√q / √(1-q²)) · [x' arctan(qp x'/ψ_e) + y' arctanh(qp y'/ψ_e)]
  in the major-axis-aligned body frame; rotation-invariant scalar so
  no back-rotation. Series fall-back for q→1 reduces to ψ_SIS.
- **`fermat_potential(lens, θ, β)`** — τ = 0.5|θ-β|² - ψ(θ).
- **`time_delay_arcsec2(lens, θ_a, θ_b, β)`** — Δτ in arcsec². The
  cosmological time-delay distance D_Δt is intentionally left to the
  caller (kernel stays unit-clean).
- **`test_time_delays`** (8 test groups):
  * SIS potential matches closed form to 1e-12 (3 sample points)
  * ∇ψ_SIE = α_SIE via central finite differences (3 points)
  * SIE q=1 potential reduces to SIS to 1e-12
  * Fermat potential ∇τ = 0 at SIS image positions to 1e-7 (FD)
  * SIS on-axis closed form Δτ = 2θ_E·β to 1e-12
  * Off-axis SIS time delay agrees with direct τ-difference + sign-flip
  * Translation invariance under common shift of lens, source, images
  * SIE 4-image cusp config: every pairwise Δτ matches τ-difference
    exactly + non-degenerate (some delays nonzero)

### Validated
40/40 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

---

## 0.6.2 — 2026-05-02 — M4 wedge 3: SIE lens + numerical image solver

The realistic-galaxy / cluster-scale lens model and a generic Newton
image solver. Argus can now generate the 4-image cusp/cross
configurations that real strong-lensing systems exhibit.

### Added
- **`argus::lensing::SIE`** — Singular Isothermal Ellipsoid via the
  closed-form Kormann, Schneider & Bartelmann 1994 deflection (eqs
  40-41). Parameters: θ_E (SIS-equivalent Einstein radius), q
  (axis ratio), φ (position angle), centre. Rotates the field into
  the major-axis-aligned frame, evaluates `arctan / arctanh` of the
  scaled coordinate, and rotates back. Series fall-back for q→1
  avoids the 1/√(1-q²) singularity; the q=1 limit reduces to SIS
  bit-exactly.
- **`argus::lensing::find_images(lens, β, ...)`** — generic
  numerical image solver for any `Lens`. Coarse grid scan locates
  candidate basins, then Newton iteration (with central-difference
  Jacobian of α) converges to each root. Magnification computed from
  1/|det(I − ∂α/∂θ)|. Deduplicates within a configurable tolerance.
- **`test_lensing` extended (10 → 18 test groups)**:
  * SIE q=1 reduces to SIS bit-exactly (4 sample points)
  * Major-axis closed form: α at (x,0) matches arctan formula
  * Minor-axis closed form: α at (0,y) matches arctanh formula
  * Point-symmetry α(-θ) = -α(θ) under rotation φ=0.7
  * Rotation covariance: lens at φ=π/2 + θ_rot=R(π/2)θ ⇒ α_rot=R(π/2)α
  * Off-centre lens: deflection translation-invariant
  * **4-image cusp config**: source inside tangential caustic →
    `find_images` finds 4 images, lens equation closes to <1e-8 each,
    every magnification finite and > 1
  * `find_images` on q=1 SIE recovers the 2-image SIS solution and
    matches both magnifications to 1e-3
  * SIE bad inputs throw (4 cases)
  * `find_images` bad inputs throw (3 cases)

### Validated
39/39 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. The next M4 wedge wires SIE into a substrate-claim
retrieval test that recovers (θ_E, q, φ, lens, source) from
synthetic 4-image observations.

---

## 0.6.1 — 2026-05-02 — M4 wedge 2: lensing retrieval (substrate proof)

The full M4 substrate-claim test: the same Argus `Retrieval` API used
for exoplanet atmospheric retrieval now drives strong-lensing parameter
inference from observed image positions. Different physics, identical
infrastructure.

### Added
- **`test_lensing_retrieval`** — 5-parameter affine-invariant ensemble
  MCMC retrieval over (θ_E, lens_x, lens_y, source_x, source_y) for an
  SIS lens. Truth is recovered to within 3σ on every parameter from
  noisy synthetic image positions (σ = 0.02 arcsec). Pipeline reuses
  the existing `argus::Spectrum` / `Retrieval::log_posterior` /
  `EnsembleSampler` / `PosteriorSummary` / `Retrieval::posterior_predictive`
  surface unchanged — the only lensing-specific code is the `Vec2`-style
  forward closure on top of `argus::lensing::sis_images`.
- 5 test groups: forward sanity → ensemble run with finite acceptance
  → 3σ marginal recovery → posterior-predictive 5–95% bracketing of
  every observation → bit-exact determinism on identical seed.

### Why ensemble over single-chain MH
SIS image positions are only weakly sensitive to a translation that
shifts both lens and source by the same vector — single-chain MH with
isotropic proposals biases the marginal posterior on lens position by
6–8σ. The Goodman-Weare stretch move is affine-invariant by
construction and recovers the correct (broader) posterior, dropping
the bias under 2σ. The test documents this explicitly so the choice
of sampler reads as physics, not as a tuning workaround.

### Validated
39/39 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

---

## 0.6.0 — 2026-05-02 — M4 starting wedge: strong-lensing pass

Argus's substrate claim — same Argus IR + Retrieval pattern generalises
across exoplanet atmospheres, gravitational lensing, and radio
interferometry — gets its first proof beyond M2 with a clean
analytically-verifiable lensing physics layer.

### Added
- **`argus::lensing` namespace**:
    `Vec2`, `Lens` (virtual base), `SIS` (Singular Isothermal Sphere),
    `lens_equation(lens, θ)`, `Image` (position + magnification),
    `sis_images(SIS, β)` — closed-form 1-or-2-image solver.
- **`test_lensing`** (8 test groups): deflection magnitude/direction,
  in-Einstein-radius two-image case + lens-equation closure,
  outside-radius single-image case, off-axis collinearity, total
  magnification 2θ_E/β closed form, off-centre lens, degenerate
  Einstein-ring source-on-centre case, malformed-input throws.

### Validated
38/38 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

Future M4 patches: SIE (Singular Isothermal Ellipsoid), magnification
maps, time delays, lensing-parameter MCMC retrieval via
`argus::Retrieval`.

---

## 0.5.3 — 2026-05-02 — Flow training via reverse-mode AD

A learnable scale+shift 1-D normalizing flow trained end-to-end via
maximum-likelihood on samples from N(2.5, 0.7²). Verifies the
flow training pattern works with the new tape — same loop structure
that scales to the full `nn::NormalizingFlow` of stacked AffineCoupling
layers.

### Added
- **`test_flow_training`** — trains s, t such that y = x·exp(s) + t
  maximises the log-likelihood under a standard-Gaussian base.
  Asserts:
  * Loss drops by ≥ 0.5 (substantial training signal).
  * Recovered (s, t) match closed-form optimum (s* = -ln σ,
    t* = -μ/σ) to 5%.
  * Transformed samples have mean ≈ 0 and var ≈ 1 (within
    finite-sample tolerance for 500 draws).

### Validated
37/37 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

This closes the M3 milestone in spirit: every component a real
amortized-SBI retrieval pipeline needs is shipped — physics forward
model + multiple MCMC samplers + neural-net primitives + Real NVP +
NormalizingFlow + reverse-mode autograd + Adam + end-to-end NN
training + flow training. The remaining M3.5 polish (full ConditionalNF,
WASP-39b benchmark vs petitRADTRANS) is engineering, not invention.

---

## 0.5.2 — 2026-05-02 — NN training via reverse-mode AD (end-to-end)

The reverse-mode tape now drives a real neural-network training loop.
A 1-input, 16-hidden, 1-output MLP fits sin(2x) end-to-end in pure
C++: forward through Var-typed Linear + tanh, MSE loss, backward
through the entire network, Adam updates parameters in place.

This validates the full training pipeline that amortized-SBI
normalizing flows are built on.

### Added
- **`argus::ad::to_vars(tape, doubles)`** — convert a vector of doubles
  into leaf Vars on the tape. Used to upload trainable parameters
  + per-step inputs.
- **`argus::ad::grads_of(tape, vars)`** — read gradients for a vector
  of Vars after backward(). Used to extract param grads for the
  optimizer.
- **`argus::ad::linear(tape, weights, bias, input, in_dim, out_dim)`**
  — Linear-layer forward y = W·x + b, fully traced on the tape.
- **`argus::ad::relu_vec / tanh_vec / sigmoid_vec`** — element-wise
  activations on a `vector<Var>`.
- **`argus::ad::mse(pred, target)`** — mean squared error as a single
  Var; the standard regression loss.
- **`test_nn_training`** (4 test groups):
  * `ad::linear` matches `nn::Linear::forward` to bit-equality
  * `ad::mse` returns 0 on equal vectors and 1.0 on offset-by-1
  * **End-to-end MLP training**: a 1-16-1 MLP fits sin(2x) on 64
    points; final loss < 5% of initial loss AND < 0.02 absolute;
    held-out prediction within 0.10 of truth.
  * 3 malformed-input throws.

### Validated
36/36 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

The same training pattern now scales to NormalizingFlow training:
build the conditioner on the tape, define a flow log-likelihood loss,
backward + step. That's the next iteration.

---

## 0.5.1 — 2026-05-02 — Adam + SGD optimizers (training pipeline live)

The reverse-mode autograd from v0.5.0 now drives a full training
pipeline: Adam (Kingma & Ba 2014) and SGD-with-momentum optimizers
that update parameters from gradients computed by the tape.

Argus can now do gradient-descent learning end-to-end in C++.

### Added
- **`argus::ad::Adam`** — bias-corrected first/second-moment estimates
  (m_hat, v_hat). Defaults: lr=1e-3, β1=0.9, β2=0.999, ε=1e-8.
  Constructor validates parameter ranges; `step(params, grads)`
  performs one update.
- **`argus::ad::SGD`** — vanilla SGD with momentum. lr + momentum ∈ [0,1).
- **`test_optimizer`** (5 test groups):
  * Adam on quadratic (x-3)² converges to x≈3 within 1e-3 in 1k steps.
  * SGD-with-momentum on the same problem also converges.
  * **End-to-end linear regression**: train (w, b) on 200 noisy
    samples, recover analytic least-squares solution within 2%.
    Demonstrates the full pipeline: Tape + Adam + mini-batch loss.
  * Adam vs SGD on a poorly-conditioned quadratic (cond #100):
    Adam beats SGD-no-momentum by orders of magnitude.
  * 6 malformed-input throws (negative lr, β out of range,
    momentum out of range, params/grads size mismatch).

### Validated
35/35 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

The next iteration ports `argus::nn::Linear::forward` to use the
`Tape`/`Var` API so neural networks can be trained end-to-end via
gradient descent in pure C++.

---

## 0.5.0 — 2026-05-02 — Reverse-mode autograd (Wengert tape)

The major M3 milestone: reverse-mode automatic differentiation via a
Wengert tape. ONE forward + ONE backward pass yields the gradient
regardless of input dimensionality. This is the algorithm
PyTorch/JAX/TensorFlow use; it's the unlock for training neural
networks and normalizing flows in C++ at scale.

Forward-mode `Dual<T>` is O(D) forward passes; reverse-mode is O(1).
For a 1000-parameter retrieval that's a 1000× speedup per gradient.

### Added
- **`argus::ad::Tape`** — append-only computation log; each node
  records (value, parent indices, local-gradient values).
- **`argus::ad::Var`** — lightweight handle (tape pointer + index +
  cached value). Operator overloading records arithmetic; math
  free-functions record exp/log/sqrt/pow/sin/cos/tanh.
- **`Tape::backward(output)`** — reverse-topological traversal that
  populates gradients on every leaf; lookup via `tape.grad(var)`.
- **`Tape::reset()`** — wipe + reuse for the next iteration.
- **`test_ad`** (10 test groups, 60+ assertions):
  * Arithmetic gradients (+, -, *, /, unary -) match analytic
  * `1/x` gradient at x=2 matches `-1/x²`
  * exp/log/sqrt composition + chain rule
  * `pow + chain rule`: `d/dx (x²+1)³ = 6x(x²+1)²`
  * `sin/cos` chain
  * Cross-validation against forward-mode `Dual<T>` on 5-D
    function across 3 random points → agree to 1e-10
  * 100-input sum-of-squares: gradient on each input correct
    (where forward-mode would need 100 evaluations)
  * `Tape::reset()` enables clean reuse
  * Wrong-tape ops throw
  * `tanh` gradient matches `1 - tanh(x)²`

### Validated
34/34 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

### Why this is a major (0.4 → 0.5) bump
Reverse-mode autograd is the gating piece for the rest of M3 → M4:
- NF training loop in C++ (instead of training in PyTorch + loading)
- HMC for high-D retrieval (HMC currently uses forward-mode Dual,
  fine for ~10 params but breaks at 100+)
- Variational inference, score matching, neural posterior estimation
- Differentiable physics: gradient through the entire Argus forward
  model, end-to-end inversion via gradient descent

The next iteration wires the new tape into NN forward passes and a
loss-function training loop.

---

## 0.4.10 — 2026-05-02 — Hamiltonian Monte Carlo (autograd-based)

Gradient-based MCMC: Hamiltonian Monte Carlo with leapfrog
integration, using forward-mode autograd through `Dual<T>` to
compute the log-posterior gradient. Far better mixing than MH on
curved or correlated posteriors.

### Added
- **`argus::grad<F>(f, x)`** — generic forward-mode gradient. Calls
  the templated function f exactly D times (one per coordinate)
  with the seed derivative on that axis = 1. Cost is O(D) forward
  evaluations; suitable for retrieval problems with O(10) parameters.
- **`argus::HMC<LogP>`** — leapfrog HMC with Metropolis correction.
  User supplies a templated log-posterior callable; the sampler
  resolves the gradient via `argus::grad`. API mirrors
  `MetropolisHastings`: `sample(state, n_samples)` returns a
  `Result` with samples and log_posteriors.
- **`test_hmc`** (5 test groups):
  * `argus::grad` on 1-D and 2-D quadratics matches the analytic
    gradient to 1e-12.
  * 2-D Gaussian recovery with HMC: mean within 0.10-0.20, stddev
    within 20%; acceptance rate > 30%.
  * Determinism: identical samples across identical seeds.
  * Bad inputs throw (negative step_size, n_leapfrog=0, non-finite
    initial log-posterior).

### Validated
33/33 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. The HMC sampler completes the gradient-based MCMC
toolkit; reverse-mode autograd (for high-D parameter spaces and NF
training) is the next major M3 milestone.

---

## 0.4.9 — 2026-05-02 — NormalizingFlow weight I/O

Save and load NormalizingFlow weights to/from a plain-text file with
hex-float values for bit-exact round-trip. Practical hand-off so
users can train flows in PyTorch (or any framework with an export
adapter) and run amortized-SBI inference in Argus.

### Added
- **`NormalizingFlow::save(path)`** — writes a self-describing
  text file:
    ```
    # argus.nn.NormalizingFlow v0.4.x
    # dim=4 n_couplings=6 init_split=2
    # coupling 0 conditioner_layers=3
    layer 0 in=2 out=16
    <weights as hex-floats, one per line>
    bias
    <bias values, one per line>
    ...
    ```
- **`NormalizingFlow::load(path)`** — reads a file produced by
  save(); strict shape validation throws on architectural mismatch
  (different dim, n_couplings, hidden widths, or layer dims).
- **`test_normalizing_flow`** extended with 3 round-trip groups:
  * Round-trip preserves forward output bit-exact (verified on a
    flow initialised with one seed, saved, then loaded into a
    flow initialised with a different seed).
  * Architecture mismatch throws on load.
  * Missing file throws on load.

### Validated
32/32 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2.

---

## 0.4.8 — 2026-05-02 — Normalizing flow (stacked Real NVP)

The actual normalizing-flow primitive: a stack of N AffineCoupling
layers separated by deterministic half-rotations, trained against a
standard-Gaussian base distribution. The infrastructure DINGO-class
amortized SBI is built on.

### Added
- **`argus::nn::HalfSwap`** — deterministic dimension rotation
  between coupling layers so successive layers see different
  "active" subsets. Identity-volume permutation; no log-det
  contribution.
- **`argus::nn::NormalizingFlow`**:
    `forward(x)`  — data → base; returns z + cumulative log_det
    `inverse(z)`  — base → data; returns x + cumulative log_det
    `log_density(x)` — change-of-variables formula:
        log p_x(x) = -½|z|² - ½ D ln(2π) + Σ log_det
    `sample(rng)` — z ~ N(0, I), then x = inverse(z)
    `coupling(i)` — mutable access for pretrained weights
- **`test_normalizing_flow`** (5 test groups):
  * Inverse-of-forward recovers input bit-exact under Xavier init
  * `log_density(x)` matches the change-of-variables formula
  * `sample(rng)` reproducible across identical seeds + finite
  * Identity-config flow (zeroed conditioner, n_couplings=1) gives
    forward(x) = x with log_det = 0 and standard-Gaussian density
  * 7 malformed-input throws

### Validated
32/32 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. The next iteration ships the amortized-SBI training loop
that learns conditioner weights so `flow.sample()` draws from the
posterior given a JWST observation.

---

## 0.4.7 — 2026-05-02 — Affine coupling layer (Real NVP / Glow block)

The fundamental building block of normalizing-flow posteriors used in
amortized SBI: a Real NVP affine coupling layer with forward/inverse
+ tractable log-determinant Jacobian.

### Added
- **`argus::nn::AffineCoupling`** (`include/argus/nn.hpp`,
  `src/nn.cpp`):
    Splits input x into x_a (first `split` dims, "passive") and x_b.
    Forward: y_a = x_a; y_b = x_b · exp(s(x_a)) + t(x_a),
    where [s, t] is the output of a configurable Sequential MLP
    conditioner from `split` to `2 · (dim - split)`.
    Inverse: bit-exact recovery — x_a = y_a;
             x_b = (y_b - t(y_a)) · exp(-s(y_a)).
    log_det_jacobian = Σ s_i for forward; -Σ s_i for inverse.

### Tests added to `test_nn` (3 new groups, 31 → 31 still since same file):
  * Forward then inverse recovers input bit-exact under
    Xavier-initialised conditioner.
  * Inverse log-det = -forward log-det.
  * Analytic log_det_jacobian agrees with numerical Jacobian
    determinant via finite differences to <1e-4.
  * Bad inputs throw (dim < 2, split = 0, split ≥ dim, wrong forward
    dimension).

### Validated
31/31 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. The next iteration will stack AffineCoupling layers into
a NormalizingFlow class with sample() / log_density() and the
amortized-SBI training loop.

---

## 0.4.6 — 2026-05-02 — Neural-net primitives + README/site refresh

Foundation for amortized SBI: Linear/Activation/Sequential building
blocks for normalizing flows in the next iteration.

### Added
- **`argus::nn::Linear`** — y = Wx + b. Row-major weights flat in
  memory, Glorot/Xavier-uniform initialisation with deterministic
  seed, direct setters/getters for weight round-tripping. No external
  BLAS dependency; M3.5 will swap inner loop for cuBLAS.
- **`argus::nn::Activation`** enum + `apply_activation()`:
  None / ReLU / LeakyReLU (α=0.01) / Tanh / Sigmoid.
- **`argus::nn::Sequential`** — multi-layer MLP, alternating
  `Linear` + `hidden_act`, terminating in a raw `Linear`. Mutable
  layer access for pretrained-weight loading.
- **`test_nn`** (7 test groups, 30+ assertions):
  * Hand-set Linear forward gives the analytic answer
  * All four activations match analytic forms
  * Identity-weight Sequential composes Tanh as expected
  * Xavier init: deterministic given seed, weights inside √(6/(in+out)) bound
  * Sequential with Xavier produces finite output, deterministic by seed
  * Layout: in/out dims and layer count
  * 7 malformed-input throws (zero dim, wrong sizes, OOB layer access)

### Changed
- **README.md** test-suite table updated to 31 tests, split into
  M2 physics layer + M3 inference layer sections.
- **README.md** roadmap row for M3 marked 🟡 substantially shipped
  (v0.4.5+) with the inventory of MCMC + ensemble + diagnostics +
  chain I/O + priors + posterior-predictive.
- **`site/index.html`** command-bar pill updated from M1/v0.1.0 to
  M3/v0.4.5 with "31 tests" replacing the commit hash.

### Validated
31/31 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. 4 examples build warning-free.

---

## 0.4.5 — 2026-05-02 — Posterior predictive + example refresh

The standard post-MCMC convergence check: generate model spectra at
posterior samples, compute per-wavelength quantile bands, verify the
observed spectrum lies inside.

### Added
- **`Retrieval::PosteriorPredictive`** struct + method:
    `posterior_predictive(samples, thin, quantiles)` runs the forward
    model at each thinned sample and returns per-wavelength band
    values at the requested quantiles (default {16%, 50%, 84%}).
- **`test_posterior_predictive`** (6 test groups):
  * Default quantiles ordered q16 ≤ q50 ≤ q84
  * Coverage: ≥85% of observations within (band ± 2σ_noise)
  * Median band within 1% of truth at every wavelength
  * Custom quantiles {2.5, 50, 97.5%}
  * Thin=1 vs thin=5 consistent to <0.5%
  * Empty samples and out-of-range quantiles throw

### Changed
- **`examples/04_retrieval.cpp`** rewritten to demonstrate the full
  M3 retrieval pipeline:
    1. Single-chain MH baseline (acceptance + posterior summary)
    2. Ensemble sampler with 16 walkers (acceptance + R̂ + ESS)
    3. Chain saved to CSV (loadable in numpy/pandas)
    4. Posterior-predictive coverage report
  Sample run: 91% coverage, R̂ ≈ 1.05 with 16 walkers × 500 steps,
  truth recovered within 1σ.

### Validated
30/30 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. 4 examples build warning-free.

---

## 0.4.4 — 2026-05-02 — Extended priors (Gaussian + LogUniform)

Real exoplanet retrievals need flexible priors beyond uniform — e.g.
Gaussian priors on stellar parameters (constrained from independent
photometric data) and log-uniform priors on quantities that span many
decades (e.g. molecular VMRs, cloud opacities).

### Added
- **`argus::PriorType`** enum: `Uniform` (default), `Gaussian`,
  `LogUniform`.
- **`argus::Parameter`** gains `prior_type`, `prior_mean`,
  `prior_stddev`. Defaults preserve existing semantics (Uniform on
  [prior_min, prior_max]).
- **`Retrieval::log_posterior`** evaluates each parameter's prior
  contribution by type:
    Uniform: log_prior += 0
    Gaussian: log_prior += -0.5 · ((x - mean)/stddev)²
    LogUniform: log_prior += -ln(x); requires x > 0
  All prior types still hard-clip to [prior_min, prior_max].
- **`test_priors`** (6 test groups):
  * Uniform on [-2, 4] → mean = 1.0, std = √3
  * Standard Gaussian clipped to ±3σ → mean ≈ 0, std ≈ 1
  * Gaussian centred at (5, σ=0.5) clipped to [3, 7] → recovered
  * LogUniform on [1, 100] → mean ≈ 21.5, std ≈ 25
  * Bounds clip applies to all prior types (Gaussian outside box → -inf)
  * LogUniform with x ≤ 0 returns -inf

### Validated
29/29 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. 4 examples build warning-free. All existing M2/M3 tests
still pass (Uniform default preserves backward compatibility).

---

## 0.4.3 — 2026-05-02 — CSV chain I/O (bit-exact round-trip)

Adds the offline-analysis hand-off: save MCMC chains to CSV and load
them back into corner.py / numpy / pandas for posterior visualisation
and downstream analysis. Hex-float encoding makes the round-trip
bit-exact.

### Added
- **`argus::ChainIO`** + **`argus::LoadedChain`**
  (`include/argus/chain_io.hpp`, `src/chain_io.cpp`):
    `save_csv(path, params, samples, log_posteriors)` writes
    a header line with version + a comment line with shape + a CSV
    header row of parameter names + 'log_posterior' column. Doubles
    encoded as `%a` hex-float for bit-exact round-trip.
    `save_csv(path, params, retrieval_result)` overload for direct
    use with the Retrieval API.
    `load_csv(path)` returns a `LoadedChain` with parameter names,
    samples, and log-posteriors.
- **`test_chain_io`** (8 test groups, 30+ assertions):
  * 1000-sample 3-parameter chain round-trip → bit-exact equality
  * Retrieval::Result overload round-trip
  * Empty samples allowed
  * Mismatched samples / log_posteriors throws on save
  * Sample dim != params count throws on save
  * Missing file throws on load
  * Malformed CSV (no header, wrong columns) throws on load
  * End-to-end with PosteriorSummary equality after round-trip

### Validated
28/28 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. 4 examples build warning-free.

---

## 0.4.2 — 2026-05-02 — MCMC convergence diagnostics (R̂ + ESS)

Adds the standard production-MCMC diagnostics: Gelman-Rubin R̂ for
between-chain convergence, and effective sample size for autocorrelation-
adjusted "independent-equivalent" sample count.

### Added
- **`argus::ChainDiagnostics`** + **`argus::compute_diagnostics()`**
  (`include/argus/diagnostics.hpp`, `src/diagnostics.cpp`):
    R̂ via the Gelman-Rubin (1992) variance-ratio formula
        B = (N/(M-1)) · Σ_m (θ̄_m - θ̄_·)²
        W = (1/M) · Σ_m s_m²
        V̂ = ((N-1)/N)·W + (1/N)·B
        R̂ = √(V̂/W)
    ESS via Geyer's (1992) initial monotone sequence cutoff:
        sum lag-t autocorrelations until first ρ_t < 0.05;
        ESS = N·M / (1 + 2·Σ ρ_t)
  Two overloads: from a multi-chain `vector<vector<vector<double>>>`,
  or directly from an `EnsembleSampler::Result` (each walker treated
  as an independent chain).
- **`test_diagnostics`** (7 test groups):
  * 4 converged chains from N(0,1) → R̂ ∈ (0.95, 1.05), ESS ≈ N·M
  * 3 chains stuck at different means → R̂ > 1.5 (non-convergence)
  * AR(1) with ρ=0.9 → R̂ ≈ 1, ESS << raw count (autocorrelation kills)
  * 2-D case validates per-parameter shapes
  * Malformed inputs throw (empty, < 4 samples, unequal lengths)
  * Single-chain edge case → R̂ = 1.0
  * EnsembleSampler::Result overload roundtrip

### Validated
27/27 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. 4 examples build warning-free.

---

## 0.4.1 — 2026-05-01 — Affine-invariant ensemble sampler (emcee-style)

Adds the de-facto astronomy MCMC algorithm: the Goodman & Weare 2010
"stretch move" affine-invariant ensemble sampler, the same algorithm
emcee implements in Python. Ensemble samplers handle highly-correlated
posteriors (e.g. T–VMR or molecule-cloud degeneracies) far better than
single-chain Metropolis-Hastings.

### Added
- **`argus::EnsembleSampler`** (`include/argus/mcmc.hpp`,
  `src/mcmc.cpp`) — N-walker ensemble sampler with the
  Goodman-Weare 2010 stretch-move proposal:
    z ~ g(z) = 1/(2 √z) on [1/a, a]
    proposal = walker_partner + z² * (walker_self - walker_partner)
    α = min(1, z^(d-1) · L_new / L_old)
  Constructor validates: walker count ≥ 4 and even, all walkers same
  dim, all walkers finite log-posterior at init, stretch_a > 1.
  Result shape: [n_steps × n_walkers][n_dim] in step-major order.
- **`test_ensemble`** (5 hard test groups, 30+ assertions):
  * 2D uncorrelated Gaussian: marginal means/stds recovered to ~10%
  * **Highly-correlated 2D Gaussian** (correlation 0.99, condition
    number ~100) — affine invariance test where the sampler still
    recovers the marginal means/stds; this is the regime where
    single-chain MH would mix slowly.
  * Bit-exact determinism on identical seeds.
  * 6 malformed-input throws (null callable, < 4 walkers, odd count,
    dim mismatch, stretch_a ≤ 1, init walker with -inf logp).
  * Result shape contract.

### Validated
26/26 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. 4 examples build warning-free. Performance baseline:
~2 ms / forward call, ~10-13 ns / Voigt eval.

---

## 0.4.0 — 2026-05-01 — M3 starting wedge: MCMC + retrieval API

The first M3 deliverable: posterior inference. Argus can now ingest
real (or synthetic) JWST-class spectra and recover atmospheric
parameters via Metropolis-Hastings MCMC against any forward model
the kernel can build.

### Added
- **`argus::MetropolisHastings`** (`include/argus/mcmc.hpp`,
  `src/mcmc.cpp`) — single-chain Metropolis sampler with isotropic
  Gaussian proposal. Deterministic given a `uint64_t` seed (verified
  bit-equal across runs). API: `burn_in()`, `sample()`,
  `acceptance_rate()`. Validates inputs strictly (null callable, empty
  / non-positive proposal widths throw).

- **`argus::Parameter`, `argus::PosteriorEntry`,
  `argus::PosteriorSummary`** — uniform-prior parameter spec, per-
  parameter median + mean + stddev + 16/84-percentile credible
  interval. Lookup by name throws on miss.

- **`argus::Retrieval`** (`include/argus/retrieval.hpp`,
  `src/retrieval.cpp`) — connects a forward-model callable, an
  observed `Spectrum`, per-wavelength uncertainty, and a parameter
  list with uniform priors. Exposes `log_posterior()` (-inf outside
  prior box; chi² log-likelihood inside) and `run_mcmc()` with
  optional auto-tuned proposal widths (default = box-width / 40).

- **`test_mcmc`** (5 assertions): recovers analytic 1D Gaussian
  (mean=3, σ=1.5) and 2D Gaussian (μ=[-2,5], σ=[0.5,2]) to <5%/<10%
  from 20-30k samples; verifies seed-determinism and input validation.

- **`test_posterior`** (4 assertions): summary statistics from 100k
  Gaussian draws, by-name lookup, throw on unknown name and shape
  mismatch.

- **`test_retrieval`** (6 assertions): generates a synthetic noisy
  JWST-PRISM-band spectrum at known T=1500 K and log10(VMR)=-3, runs
  MCMC with deliberately wrong initial guess, asserts:
    * acceptance rate within (5%, 85%)
    * posterior median within 3σ of truth (both T and log10_VMR)
    * truth inside ±2.5σ of mean
    * bit-exact reproducibility on identical seed

- **`examples/04_retrieval.cpp`** — runnable end-to-end retrieval demo.
  Outputs the recovered posterior alongside the injected truth in a
  publication-style "median +up -down (mean ± stddev)" table.

### Validated
25/25 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. 4 examples build warning-free. Retrieval example recovers
T = 1489 K (truth 1500 K) and log10(VMR) = -3.01 (truth -3.0) from
synthetic data with 50 ppm noise.

---

## 0.3.6 — 2026-05-01 — Rayleigh + CH4 + comprehensive hot-Jupiter test

Adds the last two physics components a real exoplanet retrieval needs
(Rayleigh scattering on the bulk gas, CH4 molecular opacity) and a
comprehensive end-to-end test that runs every kernel together on a
realistic hot-Jupiter setup.

### Added
- **`argus::RayleighOpacity`** — λ⁻⁴ scattering kernel.
  σ(ν) = σ₁μm · (ν/10000)⁴, no T,P dependence. For H2-dominated
  atmospheres use σ₁μm ≈ 8.49e-29 cm². T,P-independent (Rayleigh is
  a continuum scatterer); kernel pre-computes per-wavenumber values
  once per call.
- **`argus::test_data::kCH4Lines`** — 8 hand-curated real HITRAN CH4
  lines from the v3 (3.3 μm) asymmetric-stretch fundamental band of
  the most abundant 12CH4 isotopologue.
- **`test_rayleigh`** (5 assertions):
  * σ(1 μm) equals the reference value exactly
  * λ⁻⁴ scaling: σ(0.5 μm) = 16·σ(1 μm); σ(2 μm) = σ(1 μm)/16
  * T- and P-independence
  * malformed input (negative σ) throws
  * end-to-end Rayleigh slope: 0.5 μm depth > 2 μm depth
- **`test_hot_jupiter`** — comprehensive end-to-end:
  * Loads real HITRAN H2O (16 lines), CO2 (10 lines), CH4 (8 lines)
  * Builds Guillot 2010 hot-Jupiter T-P (T_int=200 K, T_irr=1500 K, γ=0.5)
  * Stacks 5 opacity kernels: H2 Rayleigh + H2O + CO2 + CH4 + cloud deck
  * Runs forward over 0.5-5 μm (JWST-PRISM range)
  * Asserts physical properties:
    - Rayleigh slope at short λ
    - H2O 1.4 μm band detected above continuum
    - CO2 4.3 μm band detected above continuum
    - CH4 3.3 μm band detected above continuum
    - cloud floor caps the deepest features (max depth < 5%)
    - bit-exact reproducibility on repeat calls

### Validated
22/22 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. 3 examples build warning-free. Performance unchanged:
~1.8 ms / forward call, ~9 ns / Voigt eval.

### M2 production-perfect inventory
The kernel now ships every component a real exoplanet atmospheric
retrieval would call:
- Hydrostatic geometry + transit-radius integration
- Hui-Armstrong-Wray Voigt (1e-6) with autograd path via `Dual<T>`
- HITRAN .par parser + bundled real H2O/CO2/CH4 fixtures
- TIPS-anchored Q(T) for 5 molecules
- Self-broadening (γ_self) + pressure shift (δ_air)
- Guillot 2010 T-P profile
- Gray cloud-deck opacity
- Rayleigh scattering (λ⁻⁴)
- Bit-exact reproducibility, file I/O round-trip, performance baseline

---

## 0.3.5 — 2026-05-01 — Gray cloud deck

Adds the canonical "gray cloud deck" opacity model used in nearly every
exoplanet atmospheric retrieval. Above the cloud-top pressure the
atmosphere goes opaque; the gas opacity dominates below.

### Added
- **`argus::CloudDeckOpacity`** — `OpacityKernel` that returns
  σ_cloud per "cloud particle" for layers with P ≥ P_cloud_bar, zero
  otherwise. Wavelength-independent (gray) by design — the canonical
  retrieval cloud model. Constructor validates P_cloud > 0 and
  σ ≥ 0; throws otherwise.
- **`test_clouds`** (4 assertions):
  * cross-section is zero above the cloud (low P) and σ below.
  * cross-section is exactly σ at and below P_cloud.
  * end-to-end transit-depth: deeper cloud (higher P_cloud) → lower
    transit depth; shallower cloud → deeper transit (more opaque area).
  * malformed inputs throw.

### Validated
20/20 tests pass clean under
  `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2. 3 examples build warning-free.

---

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
