# Claude project notes — Argus

## What this repo is

Argus is a **differentiable astrophysical inference kernel** — the same
C++20 substrate is used for atmospheric retrieval (M2/M3, shipped),
strong-gravitational lensing (M4, shipped), and radio interferometry
(M5, shipped). Read `Argus.tex` (the one-page brief) and
`docs/Astronomy-Compute-Crisis.pdf` (the 23-page research deck) before
suggesting architectural changes.

## Architectural commitments (don't drift)

- **C++20 kernel, Python edges.** The hot path stays C++. Python bindings
  (planned via pybind11) sit above, never inside.
- **Stable IR is sacred.** `include/argus/ir.hpp` is the integration contract
  for every future pass (cosmology, photometry, GW). Adding a new
  `NodeKind` is fine; renaming or repurposing existing kinds breaks
  downstream content addresses.
- **Pluggable opacity sources.** `OpacityKernel` is the substitution point for
  HITRAN, HITEMP, ExoMol, neural emulators. Don't leak source-specific types
  into the IR or the RT layer.
- **Content-addressed runs.** Every retrieval is reproducible bit-for-bit.
  Do not introduce nondeterminism (uninitialized memory, wall-clock
  seeds, hash randomization) in the kernel. The seven retrieval
  substrate-proof tests verify bit-equal-determinism on a fixed seed.
- **One substrate, three physics layers.** Atmospheric / lensing /
  interferometry retrievals all reuse the same `Spectrum` /
  `Retrieval::log_posterior` / `EnsembleSampler` / `MetropolisHastings`
  / `HMC` / `PosteriorSummary` / `posterior_predictive` infrastructure —
  only the forward closure changes. Don't add domain-specific code into
  the inference layer.
- **CMake-only build.** No external runtime deps. CUDA, pybind11 and
  HITEMP/ExoMol loaders are all on the wishlist (see site WANT block);
  the public surface absorbs them without a reshape.

## Working principles

- Prefer editing existing files over creating new ones.
- One short comment when the *why* is non-obvious — never narrate *what*.
- No backwards-compat shims yet — pre-1.0, break the interface freely.
- Tests use plain `assert()`, no Catch2/GoogleTest dependency.
- When introducing a new pass: add the header in `include/argus/`, the impl in
  `src/`, a smoke test in `tests/`, and a one-paragraph design rationale at
  the top of the header.

## Roadmap reference

See README.md "Roadmap" table. **Current state: M1–M5 all shipped (v0.7.18).**

50/50 hard tests pass under
`-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -O2`
on GCC 15.2; 50/50 also pass clean under
`-O1 -fsanitize=address,undefined`.

**Open wishlist** (not yet shipped):
- M2.5: CUDA-resident Voigt opacity kernel
- M3.5: ConditionalNF training loop (architecture shipped v0.7.9);
  WASP-39b multi-molecule (Na+H2O+CO2+CO) full-PRISM fit
  (real-data infra + H2O-only retrieval shipped v0.7.17)
- M4.5: real-data fixtures (HE 0435-1223 quad-lens positions);
  AD-templated lens deflection so HMC can climb lens posteriors
- M5.5: CLEAN-style image reconstruction; self-cal; UV-FITS /
  Measurement-Set ingest
- DOC: architecture diagrams + design rationale beyond the deck
- PY: pybind11 bindings for the full surface
