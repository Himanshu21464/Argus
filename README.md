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
│   ├── radiative_transfer.hpp       TransmissionModel (forward pass)
│   └── ir.hpp                       Argus IR (typed physics graph + content addressing)
├── src/                             implementations
├── tests/                           assert-based smoke tests
├── examples/
│   └── 01_transmission_spectrum.cpp first end-to-end demo
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

| Phase | Window | Deliverable |
|-------|--------|-------------|
| **M1** | months 1–3 | Argus IR, atmosphere/opacity/RT scaffolding (this repo's starting point) |
| **M2** | months 3–6 | GPU-resident HITRAN/HITEMP/ExoMol opacity, Voigt CUDA kernel, differentiable transmission solver, WASP-39b benchmark vs. petitRADTRANS |
| **M3** | months 6–9 | Amortized SBI (normalizing flows), 10× speed vs. POSEIDON/CHIMERA on public JWST spectra without binning |
| **M4** | months 9–12 | Lensing pass — proves the IR generalizes |
| **M5** | months 12–18 | Interferometric imaging pass — three sub-fields, one kernel, substrate claim |

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
