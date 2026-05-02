# Argus

**A differentiable atmospheric retrieval kernel.**

> "LLVM for the sky — starting with a single planet."

Argus is a C++20/CUDA core for differentiable line-by-line radiative transfer of
exoplanet atmospheres, with an autograd-aware intermediate representation and a
built-in simulation-based inference engine. The wedge into a substrate that
eventually compiles photons into physics across exoplanets, lensing, and
interferometry.

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
│   └── ir.hpp                       Argus IR (typed physics graph + content addressing)
├── src/                             implementations
├── tests/                           assert-based hard tests (30 tests)
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
| **M3** | months 6–9 | 🟡 substantially shipped (v0.4.9) | MCMC + emcee ensemble + Retrieval API + R̂/ESS diagnostics + CSV chain I/O + Gaussian/LogUniform/Uniform priors + posterior-predictive checks + neural-net primitives (Linear/Activation/Sequential) + AffineCoupling (Real NVP) + NormalizingFlow with sample / log_density / save / load. Outstanding: NF training loop (needs reverse-mode autograd), WASP-39b benchmark vs. POSEIDON/CHIMERA. |
| **M4** | months 9–12 | ⏳ planned | Lensing pass — proves the IR generalizes |
| **M5** | months 12–18 | ⏳ planned | Interferometric imaging pass — three sub-fields, one kernel, substrate claim |

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
32 tests · all pass under `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`:

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
