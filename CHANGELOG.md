# Changelog

## 0.2.0 — 2026-05-01 — M2 starting wedge

First piece of M2: real physics in place of the M1 chord-stub, plus the
forward-mode autograd primitive.

### Added
- **`argus::Geometry`** (`include/argus/geometry.hpp`, `src/geometry.cpp`)
  — hydrostatic ray geometry. `build_geometry()` derives layer altitudes
  from the local scale height; `chord_path_length()` returns the geometric
  path length of a tangent ray through one spherical shell.
- **`argus::voigt`** (`include/argus/voigt.hpp`) — header-only,
  area-normalised Voigt line-shape via the Thompson-Cox-Hastings (1987)
  pseudo-Voigt approximation. Templated so the same code path runs under
  `Dual<double>` for forward-mode autograd. `gaussian()` and `lorentz()`
  fallbacks alongside.
- **`argus::Line`, `argus::LineListOpacity`** (`include/argus/line_list.hpp`,
  `src/line_list.cpp`) — HITRAN-subset line records and an `OpacityKernel`
  that sums Voigt-shaped contributions from a list of lines. Includes
  Doppler width + air-broadened Lorentz HWHM with HITRAN T-dependence
  exponent, and an E_lower-only intensity-temperature scaling (M3 will
  add full partition-function ratios).
- **`argus::Dual<T>`** (`include/argus/dual.hpp`) — forward-mode autograd
  dual numbers. Arithmetic, `exp`, `log`, `sqrt`, `pow(a,T)`, `sin`, `cos`.
  Composes through the templated `voigt()` so derivatives propagate
  through the line-shape evaluator.
- **Atmosphere** gained `bulk_mmw_amu` (default 2.3, H2/He hot-Jupiter) and
  `star_radius_rsun` (default 1.0, solar-type host) for the geometry pass.
- New tests: `test_geometry`, `test_voigt`, `test_line_list`, `test_dual`
  (4 new, 3 existing — all 7 pass under `-Wall -Wextra -Wpedantic
  -Wshadow -Wconversion -Wsign-conversion`).
- `examples/02_voigt_h2o.cpp` — 4-line H2O-like list, prints the
  resulting transmission spectrum (now with absorption peaks at the line
  centres, not flat). Demonstrates `Dual<double>` through the Voigt call:
  `dV/d(gamma_l)` at line centre.

### Changed
- **`TransmissionModel::forward`** now does proper transit-radius
  integration: hydrostatic geometry + chord through concentric shells +
  `R_eff^2 = R_p^2 + 2 * sum_b (1 - exp(-tau(b))) * b * db`. Spectrum
  values are now physical transit depths `(R_eff / R_star)^2` (M1 used a
  flat 1-km path stub).
- Version → 0.2.0.

### Roadmap context
- M2 ⅓ delivered. Outstanding for M2 completion: HITRAN/HITEMP/ExoMol
  loaders, GPU residency for opacity tables, CUDA Voigt kernel, WASP-39b
  benchmark vs. petitRADTRANS.

---

## 0.1.0 — 2026-05-01 — M1 starting point

Initial scaffolding.

- Argus IR skeleton: `Node`, `Graph`, content-address (FNV-1a 64-bit stub).
- Atmosphere model: layered T–P–VMR with isothermal builder.
- Opacity interface: `OpacityKernel` virtual base, `GreyOpacity` placeholder.
- Transmission forward model: chord-integrated optical depth (M1 stub).
- C++20 build via CMake 3.20+, zero runtime dependencies.
- Smoke tests for atmosphere, RT, IR content-addressing.
- Example: `examples/01_transmission_spectrum.cpp` — end-to-end demo.
- `Argus.tex` one-page Future Project Brief (matches Vita/Kepler portfolio).
- `docs/Astronomy-Compute-Crisis.tex` — 23-page research deck on the top
  10 computational problems in 2026 astronomy and the kernel-level
  solutions. Argus is positioned as the wedge into the substrate.
- `site/` — astronomy-themed landing page (later redesigned with native
  astronomy vocabulary in commit 45bb494).
