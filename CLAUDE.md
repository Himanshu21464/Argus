# Claude project notes — Argus

## What this repo is

Argus is the **differentiable atmospheric retrieval kernel** wedge of a larger
"LLVM for astronomy" substrate play. Read `Argus.tex` (the one-page brief) and
`docs/Astronomy-Compute-Crisis.pdf` (the 23-page research deck) before
suggesting architectural changes.

## Architectural commitments (don't drift)

- **C++20 kernel, Python edges.** The hot path stays C++. Python bindings
  (planned via pybind11) sit above, never inside.
- **Stable IR is sacred.** `include/argus/ir.hpp` is the integration contract
  for every future pass (lensing, interferometry, cosmology). Adding a
  new `NodeKind` is fine; renaming or repurposing existing kinds breaks
  downstream content addresses.
- **Pluggable opacity sources.** `OpacityKernel` is the substitution point for
  HITRAN, HITEMP, ExoMol, neural emulators. Don't leak source-specific types
  into the IR or the RT layer.
- **Content-addressed runs.** Every retrieval is reproducible bit-for-bit. Do
  not introduce nondeterminism (uninitialized memory, wall-clock seeds, hash
  randomization) in the kernel.
- **No external runtime deps in M1.** CMake-only build. CUDA arrives in M2;
  HITRAN loaders arrive in M2; pybind11 arrives in M3.

## Working principles

- Prefer editing existing files over creating new ones.
- One short comment when the *why* is non-obvious — never narrate *what*.
- No backwards-compat shims yet — pre-1.0, break the interface freely.
- Tests use plain `assert()`, no Catch2/GoogleTest dependency in M1.
- When introducing a new pass: add the header in `include/argus/`, the impl in
  `src/`, a smoke test in `tests/`, and a one-paragraph design rationale at
  the top of the header.

## Roadmap reference

See README.md "Roadmap" table. Current phase: **M1 (months 1–3)** — IR + scaffold.
Next milestone: HITRAN-backed CUDA Voigt opacity kernel.
