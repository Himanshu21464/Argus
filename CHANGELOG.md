# Changelog

## 0.1.0 — 2026-05-01

Initial scaffolding (M1 starting point).

- Argus IR skeleton: `Node`, `Graph`, content-address (FNV-1a 64-bit stub).
- Atmosphere model: layered T–P–VMR with isothermal builder.
- Opacity interface: `OpacityKernel` virtual base, `GreyOpacity` placeholder.
- Transmission forward model: chord-integrated optical depth (M1 stub).
- C++20 build via CMake 3.20+, zero runtime dependencies.
- Smoke tests for atmosphere, transmission, IR content-addressing.
- Example: `examples/01_transmission_spectrum.cpp` — end-to-end demo.
- `Argus.tex` one-page project brief (matches Vita/Kepler portfolio style).
- `docs/Astronomy-Compute-Crisis.tex` 23-page research deck on the top 10
  computational problems in 2026 astronomy and the kernel-level solutions.
