# Argus

**A differentiable astrophysical inference kernel.**

> "LLVM for the sky — starting with a single planet."

Argus is a C++20 core for differentiable forward modelling and posterior
inference across three astrophysics domains:

* **Exoplanet atmospheric retrieval** — line-by-line radiative transfer
  with HITRAN-grade Voigt opacity (M2/M3, shipped).
* **Strong gravitational lensing** — SIS/SIE/external-shear lens models
  with Fermat-potential time delays (M4, shipped).
* **Radio interferometry** — visibility forward model on synthesised UV
  coverage (M5, shipped).

The substrate claim — *the same `argus::Spectrum` / `Retrieval` API
recovers parameters across all three physics layers* — is demonstrated
end-to-end by **six independent retrieval tests**:

| Test | Layer | Free params | Recovery | Sampler |
|---|---|---|---|---|
| `test_retrieval` | atmospheric isothermal (M2/M3) | 2 | 1σ | MH |
| `test_atmosphere_retrieval_full` | **atmospheric multi-physics** (Guillot + H₂O + CO₂ + CH₄ + cloud + Rayleigh) | 4 | 3σ | MH |
| `test_hmc_atmosphere` | atmospheric — differentiable-physics HMC | 2 | 3σ | HMC + Dual<T> autograd |
| `test_lensing_retrieval` | SIS lensing (M4) | 5 | 3σ | Ensemble |
| `test_sie_retrieval` | SIE lensing (M4) | 7 | 0.5σ | Ensemble |
| `test_interferometry_retrieval` | radio interferometry (M5) | 4 | 3σ | MH |
| `test_multi_component_retrieval` | 2-component interferometry (M5) | 8 | 3σ | Ensemble |

All seven reuse the same `Spectrum` / `Retrieval::log_posterior` /
`EnsembleSampler` / `MetropolisHastings` / `HMC` / `PosteriorSummary` /
`posterior_predictive` infrastructure — only the forward model changes.

## Why this exists

JWST is producing native-resolution transmission and emission spectra that
current retrievals (TauREx, POSEIDON, CHIMERA, petitRADTRANS, ExoJAX) **bin
away** to keep MCMC tractable — losing the very signal that argued for an $11B
telescope. Ariel launches in 2029 expecting ~1000 atmospheres; today's tools
cannot keep up.

There is no `libRT` for exoplanets. There is no shared IR. There is no
GPU-native, autograd-aware production kernel. Argus aims to be that substrate.

## Project layout

```
Argus/
├── Argus.tex / Argus.pdf            one-page Future Project Brief
├── CMakeLists.txt                   C++20 build
├── include/argus/                   public kernel headers
│   ├── argus.hpp                    umbrella include
│   ├── version.hpp
│   ├── tensor.hpp                   dense tensor (CPU now, GPU next)
│   ├── atmosphere.hpp               layered T-P-VMR atmosphere
│   ├── opacity.hpp                  OpacityKernel interface + GreyOpacity stub
│   ├── line_list.hpp                Line + LineListOpacity (HITRAN-subset)
│   ├── hitran.hpp                   fixed-width .par parser
│   ├── partition.hpp                TIPS-anchored Q(T)
│   ├── test_data.hpp                bundled real HITRAN H2O lines for tests
│   ├── voigt.hpp                    Hui-Armstrong-Wray Voigt (~1e-6, templated)
│   ├── geometry.hpp                 hydrostatic ray geometry, chord paths
│   ├── radiative_transfer.hpp       TransmissionModel (proper transit-radius integration)
│   ├── dual.hpp                     forward-mode autograd dual numbers
│   ├── mcmc.hpp                     Metropolis-Hastings + emcee ensemble sampler
│   ├── retrieval.hpp                Retrieval API + Parameter / Prior / PosteriorSummary
│   ├── diagnostics.hpp              Gelman-Rubin R̂ + effective sample size
│   ├── chain_io.hpp                 CSV save/load for MCMC chains (hex-float bit-exact)
│   ├── nn.hpp                       neural-net primitives + Real NVP coupling + NormalizingFlow
│   ├── ad.hpp                       reverse-mode autograd (Wengert tape) + Adam/SGD optimizers
│   ├── lensing.hpp                  SIS / SIE / ExternalShear / CompoundLens + image solver + Fermat time delays
│   ├── interferometry.hpp           visibility forward (Point + Gaussian sources) + UV coverage
│   └── ir.hpp                       Argus IR (typed physics graph + content addressing)
├── src/                             implementations
├── tests/                           assert-based hard tests (49 tests, incl. 2 perf benchmarks)
├── examples/
│   ├── 01_transmission_spectrum.cpp first end-to-end demo (grey opacity)
│   ├── 02_voigt_h2o.cpp             4-line H2O-like spectrum + autograd demo
│   ├── 03_real_hitran.cpp           16 real HITRAN H2O lines, JWST-PRISM-shaped spectrum
│   └── 04_retrieval.cpp             full M3 pipeline: ensemble MCMC + R̂/ESS + CSV + posterior predictive
├── site/                            astronomy-themed project site (port 8767)
└── docs/
    ├── Astronomy-Compute-Crisis.tex landscape research deck (top-10 problems)
    └── Astronomy-Compute-Crisis.pdf compiled PPT-style PDF
```

## Build

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build
./build/examples/example_01_transmission_spectrum
```

Requires a C++20 compiler (GCC 11+, Clang 14+) and CMake 3.20+.
No external runtime dependencies for the M1 kernel.

## Roadmap (18-month MVP)

| Phase | Window | Status | Deliverable |
|-------|--------|--------|-------------|
| **M1** | months 1–3 | ✅ shipped v0.1.0 | Argus IR, atmosphere/opacity/RT scaffolding |
| **M2** | months 3–6 | ✅ shipped v0.3.0 | Hydrostatic geometry, Hui-Armstrong-Wray Voigt (~1e-6), LineListOpacity, dual-number autograd, HITRAN .par parser, TIPS partition functions, real-data tests, finite-diff autograd validation. CUDA residency (M2.5) and WASP-39b benchmark vs. petitRADTRANS (M3 wedge) outstanding. |
| **M3** | months 6–9 | ✅ shipped v0.5.3 | MCMC + emcee ensemble + HMC (autograd-gradient) + Retrieval API + R̂/ESS + CSV I/O + 3 prior types + posterior-predictive + nn (Linear/Activation/Sequential) + AffineCoupling (Real NVP) + NormalizingFlow with save/load + reverse-mode autograd (Wengert tape) + Adam/SGD optimizers + end-to-end NN training (MLP fits sin(2x)) + end-to-end flow training (scale+shift recovers analytic optimum). Remaining for M3.5: ConditionalNF for true amortized SBI, WASP-39b benchmark vs. POSEIDON/CHIMERA. |
| **M4** | months 9–12 | ✅ shipped v0.7.5 | Strong-lensing pass: SIS + SIE + NFW + ExternalShear + CompoundLens + numerical image solver + lensing potential + Fermat time delays. Two substrate-claim retrievals: SIS recovers (θ_E, lens, source) to 3σ via `EnsembleSampler`; **SIE capstone** recovers all 7 SIE+source params (including q and φ) to 0.5σ via source-plane chi². |
| **M5** | months 12–18 | ✅ shipped v0.7.6 | Radio-interferometry pass: PointSource + GaussianSource + visibility predictor + UV-coverage primitive (snapshot + Earth-rotation track). Two substrate-claim retrievals: 4-param Gaussian-source recovery to 3σ via single-chain MH; **2-component capstone** recovers 8 params (compact core + extended jet) over a 5-hour Earth-rotation track to 3σ via `EnsembleSampler`. |

### What v0.3.5 (M2-complete) adds over v0.2.0
- **Verified Hui-Armstrong-Wray Voigt** (~1e-6 accuracy) replacing pseudo-Voigt (1%).
  Tested against analytic Gaussian/Lorentzian limits, numerical convolution,
  Lorentz asymptotic wing, dual-number derivatives vs central finite differences,
  and an independent closed-form reference (`std::erfc`, 1e-15 accurate).
- **`argus::Hitran`** — fixed-width 160-char HITRAN .par parser with file I/O.
- **`argus::Partition`** — TIPS-anchored Q(T) for H2O, CO2, CH4, CO, NH3.
- **Full HITRAN intensity scaling** — `S(T) = S(296) · Q(296)/Q(T) · exp(-c2 E"·Δ(1/T)) · induced_emission`.
- **Self-broadening + pressure shift**: γ_total = γ_air·(1-VMR) + γ_self·VMR; nu_eff = nu0 + δ_air·P_atm.
- **Bundled real HITRAN test fixtures** — 16 H2O lines (2.7 μm + 1.4 μm bands) + 10 CO2 lines (4.3 μm band).
- **`argus::guillot()`** — Guillot 2010 analytic hot-Jupiter T-P profile.
- **`argus::CloudDeckOpacity`** — canonical gray cloud-deck retrieval model.
- **9 new hard tests** since v0.2.0 (4→20 total) — see CHANGELOG.
- **Performance baseline**: ~1.7 ms / forward call, ~9 ns / Voigt evaluation.

### Test suite
49 tests · all pass under `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`:

**Physics layer (M2):**

| test | what it asserts |
|---|---|
| `test_atmosphere` | layer ordering, isothermal builder validation |
| `test_radiative_transfer` | spectrum shape, grey-opacity flat profile |
| `test_ir` | content-address determinism, graph topology |
| `test_geometry` | chord path = 2r at b=0, hydrostatic z(P) to 1e-9 |
| `test_voigt` | analytic limits, symmetry, area-normalisation, numerical convolution, asymptotic wing, autograd vs FD |
| `test_voigt_reference` | independent closed-form check via `std::erfc` (1e-15 ref) |
| `test_line_list` | line-shape ordering, sum rules |
| `test_self_broadening` | γ_self contribution + δ_air pressure shift |
| `test_dual` | arithmetic, exp/log/sqrt, chain rule |
| `test_partition` | TIPS anchors at 296/1000/2000 K, monotonicity, throw on bad input |
| `test_hitran` | round-trip parse, CR/LF tolerance, malformed rejection, filter |
| `test_real_h2o` | end-to-end with real HITRAN: sum rule, T-dependence, VMR monotonicity, saturation |
| `test_multi_molecule` | H2O + CO2 in one atmosphere, additivity, band positions |
| `test_finite_diff` | autograd ∂/∂σ, ∂/∂γ, ∂/∂x, ∂/∂S vs central FD to 1e-5 |
| `test_stress` | 15 edge cases: malformed, single-layer, extreme T/P, empty, garbage |
| `test_file_io` | HITRAN .par round-trip via `/tmp` file |
| `test_reproducibility` | bit-exact forward-model output across calls |
| `test_benchmark` | wall-time baseline (~1.7 ms / forward call, ~9 ns / Voigt) |
| `test_guillot` | Guillot 2010 T-P: skin temperature, monotonicity, hot/cold transit |
| `test_clouds` | gray cloud deck: zero above P_cloud, opaque below, transit-depth ordering |
| `test_rayleigh` | λ⁻⁴ scaling exact, T,P-independence, end-to-end Rayleigh slope |
| `test_hot_jupiter` | comprehensive: Guillot + H2O+CO2+CH4 + cloud + Rayleigh — all bands detected |
| `test_atmosphere_retrieval_full` | **M2 capstone**: 4-param MH retrieval through the full multi-physics forward (Guillot + H₂O + CO₂ + CH₄ + cloud + Rayleigh); recovers (T_irr, log₁₀VMR_H₂O, log₁₀VMR_CO₂, log₁₀P_cloud) within 3σ in 1.8 s |

**Inference layer (M3):**

| test | what it asserts |
|---|---|
| `test_mcmc` | Metropolis-Hastings recovers analytic 1D and 2D Gaussians; determinism |
| `test_ensemble` | Goodman-Weare 2010 stretch move; affine-invariance on highly-correlated 2D |
| `test_diagnostics` | Gelman-Rubin R̂ (converged ≈ 1.0, non-converged > 1.5) + ESS via Geyer's cutoff |
| `test_chain_io` | CSV chain round-trip via hex-float — bit-exact equality |
| `test_priors` | Uniform / Gaussian / LogUniform prior shapes recovered via MCMC |
| `test_posterior` | Per-parameter median/mean/stddev/16-84% quantiles |
| `test_retrieval` | End-to-end retrieval: recover injected T + log10(VMR) within 3σ |
| `test_posterior_predictive` | 16-84% spectrum bands cover ≥ 85% of observations within 2σ_noise |
| `test_nn` | Linear/Activation/Sequential + AffineCoupling forward/inverse + analytic log-det vs FD |
| `test_normalizing_flow` | Stacked Real NVP: bit-exact roundtrip, log_density vs change-of-variables, sample reproducibility, save/load round-trip |
| `test_hmc` | Hamiltonian MC with leapfrog + autograd gradient (forward-mode); recovers analytic Gaussian |
| `test_ad` | Reverse-mode autograd (Wengert tape); cross-validates against forward-mode Dual<T> on 5-D function |
| `test_optimizer` | Adam + SGD: minimise (x-3)²; train (w, b) on noisy linear regression to within 2% |
| `test_nn_training` | Train a 1-16-1 MLP to fit sin(2x) on 64 points end-to-end via Adam |
| `test_flow_training` | Train scale+shift normalizing flow on N(2.5, 0.7²) — recovers (s, t) to closed-form optimum |
| `test_conditional_flow` | **M3.5**: ConditionalAffineCoupling + ConditionalNormalizingFlow — forward/inverse round-trip, log_density consistency, conditional dependence (different y → different density), determinism. The architecture amortized SBI is built on (NPE, SNPE, DINGO). |
| `test_hmc_atmosphere` | **M3 differentiable-physics proof**: Hamiltonian MC end-to-end through a 3-line H₂O-like isothermal atmospheric forward — Voigt + temperature-scaled widths + chi² — with leapfrog gradients via forward-mode autograd (Dual<T>). Recovers (T_K, log₁₀ VMR) within 3σ at >50% acceptance. |

**Strong-lensing layer (M4):**

| test | what it asserts |
|---|---|
| `test_lensing` | SIS / SIE / NFW / ExternalShear / CompoundLens — 29 test groups: deflection / potential closed forms, lens-equation closure, 4-image cusp config via numerical solver, NFW Wright-Brainerd h(1) limit, find_images on shear-broken configurations |
| `test_lensing_retrieval` | **Substrate proof (SIS)**: same `Retrieval` API recovers (θ_E, lens_x, lens_y, source_x, source_y) via `EnsembleSampler` to 3σ; posterior-predictive 5–95% band brackets every observation; bit-exact determinism |
| `test_sie_retrieval` | **Substrate proof (SIE capstone)**: 7 free params (θ_E, q, φ, lens, source) recovered to **< 0.5σ** from 4 observed quad-image positions via source-plane chi²; 32×4000 ensemble sampler; bit-exact determinism |
| `test_time_delays` | SIS potential closed form; ∇ψ_SIE = α via FD; SIE q=1 potential reduces to SIS; Fermat ∇τ=0 at lens-equation roots; SIS on-axis Δτ = 2θ_E·β; off-axis numerical agreement; translation invariance; SIE 4-image cusp pairwise τ-diffs |

**Radio-interferometry layer (M5):**

| test | what it asserts |
|---|---|
| `test_interferometry` | 15 test groups: PointSource at origin → V = F; off-origin phase = -2π(u·l + v·m); visibility additivity; translation theorem; conjugate symmetry; Gaussian decay exp(-2π² σ² r²); σ=0 reduces to point bit-exactly; UV-coverage snapshot; **Earth-rotation track ellipse identity** u²+(v/sin δ)² = (B_E/λ)² |
| `test_interferometry_retrieval` | **Substrate proof (single Gaussian)**: same `Retrieval` API recovers (l, m, F, σ) on a 7-antenna VLA-like array (21 baselines, 42 components) to 3σ; posterior-predictive coverage ≥ 85%; bit-exact determinism |
| `test_multi_component_retrieval` | **Substrate proof (2-component capstone)**: 8 free params (compact core + extended jet) recovered to 3σ from 210 visibility components sampled across a 5-hour Earth-rotation track on a 7-antenna VLA-like array; label-switching broken by non-overlapping size priors |
| `test_lensing_benchmark` | Wall-time baseline: SIE deflection ≈ 40 ns; NFW ≈ 26 ns; CompoundLens(SIS+γ) ≈ 6 ns; find_images(SIE 60²) ≈ 1.2 ms |
| `test_interferometry_benchmark` | Wall-time baseline: predict_visibility ≈ 13 ns; predict_visibilities(2-Gauss × 10 530 UV) ≈ 12 ns/eval; uv_coverage_track(27 ant × 30 HA) ≈ 22 µs |

## Design principles

1. **C++ kernel, Python edges.** Production hot path is C++20/CUDA. Python bindings (pybind11) live above the kernel for retrieval scripts, agent reasoning, plotting, UI.
2. **Differentiable end-to-end.** Every forward model on the kernel side is autograd-aware. SBI, normalizing-flow posteriors, and gradient-based MAP all hang off the same tape.
3. **Pluggable opacity sources.** `OpacityKernel` is an interface; HITRAN, HITEMP, ExoMol, neural emulators all implement it without leaking into the IR.
4. **Content-addressed runs.** Every retrieval is a frozen graph + frozen opacity + frozen sampler. Bit-reproducible across hardware and years.
5. **Open under Apache 2.0.** A kernel with a closed core has no ecosystem.

## License

Apache 2.0. See `LICENSE`.

---

*Future Project Brief — 2026.* The accompanying research deck
(`docs/Astronomy-Compute-Crisis.pdf`) covers the top 10 computational problems
in 2026 astronomy, the AI/kernel solutions, and the reasoning behind picking
exoplanet atmospheric retrieval as the Argus wedge.
