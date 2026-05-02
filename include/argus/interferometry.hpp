#pragma once

#include <vector>

namespace argus::interferometry {

// Radio interferometry forward model — predicts complex visibilities
// V(u, v) = ∫∫ I(l, m) exp(-2πi (u·l + v·m)) dl dm
// from a small-angle sky brightness distribution I(l, m).
//
// Coordinates:
//   * (u, v) in wavelengths (i.e. baseline length divided by λ)
//   * (l, m) in radians (direction cosines, small-angle limit)
//   * flux  in arbitrary units (Jy if u/v are wavelengths)
//
// M5's substrate claim: the same Argus Retrieval pattern that recovers
// atmospheric T, VMR from JWST spectra and lens parameters from image
// positions also recovers source positions / fluxes / sizes from
// interferometric visibilities. Only the physics layer (this file)
// changes.

// One UV-plane sample (a single baseline at one integration time).
struct UVPoint {
  double u = 0.0;   // East-West component, wavelengths
  double v = 0.0;   // North-South component, wavelengths
};

// Complex visibility V = real + i · imag.
struct Visibility {
  double real = 0.0;
  double imag = 0.0;
};

// A point source — Dirac δ in (l, m). Visibility magnitude is constant
// across the UV plane: V(u,v) = F · exp(-2πi (u·l + v·m)).
struct PointSource {
  double l    = 0.0;   // direction cosine, radians (E)
  double m    = 0.0;   // direction cosine, radians (N)
  double flux = 1.0;   // in flux units consistent with the visibilities
};

// A 2-D circular Gaussian source. The 1-σ half-width is `sigma` in
// radians (FWHM = 2·√(2 ln 2)·σ ≈ 2.355 σ). Visibility:
//   V(u,v) = F · exp(-2π² σ² (u²+v²)) · exp(-2πi (u·l + v·m))
// — the FT of a 2D Gaussian is another Gaussian (in UV) times the
// position phase.
struct GaussianSource {
  double l     = 0.0;
  double m     = 0.0;
  double flux  = 1.0;
  double sigma = 0.0;     // radians; sigma=0 reduces to PointSource
};

// Predict the visibility of one component at one UV point.
Visibility predict_visibility(const PointSource&    src, UVPoint uv);
Visibility predict_visibility(const GaussianSource& src, UVPoint uv);

// Predict visibilities at every UV point — summed over a list of
// components. Two overloads: for point sources only, and for arbitrary
// Gaussian components (which include points as the σ=0 limit).
std::vector<Visibility> predict_visibilities(
    const std::vector<PointSource>& sources,
    const std::vector<UVPoint>&     uv_points);

std::vector<Visibility> predict_visibilities(
    const std::vector<GaussianSource>& sources,
    const std::vector<UVPoint>&        uv_points);

// Synthesise the UV coverage of an interferometric array. Each pair
// of antennas contributes one baseline; here we evaluate the snapshot
// (u, v) of every pair by treating antenna positions as East-North
// metres and dividing by wavelength. For a multi-snapshot observation
// the caller can union the outputs.
//
// Returns N·(N-1)/2 UVPoints — one per antenna pair, in lexicographic
// (i, j) order with i < j.
std::vector<UVPoint> uv_coverage_snapshot(
    const std::vector<double>& antenna_east_m,
    const std::vector<double>& antenna_north_m,
    double                     wavelength_m);

// Earth-rotation UV synthesis (the standard real-array trick): as
// the Earth rotates, projected baselines sweep tracks in the UV
// plane. For a source at (HA, dec) seen from a station at latitude L
// the local ENU baseline (B_E, B_N, 0) is rotated into equatorial
// coords then projected onto the UV plane perpendicular to the source
// direction:
//
//     B_X = -sin(L) · B_N,    B_Y = B_E,    B_Z = +cos(L) · B_N
//     u =  sin(h) · B_X + cos(h) · B_Y
//     v = -sin(δ) cos(h) · B_X + sin(δ) sin(h) · B_Y + cos(δ) · B_Z
//
// (u, v) are then divided by wavelength.
//
// Returns N·(N-1)/2 · N_HA UVPoints in HA-major order: all baselines
// at hour_angles[0], then all at hour_angles[1], ... All antennas are
// assumed to be at altitude 0 (no B_U); for real arrays this is a
// negligible approximation at the array footprint.
std::vector<UVPoint> uv_coverage_track(
    const std::vector<double>& antenna_east_m,
    const std::vector<double>& antenna_north_m,
    double                     latitude_rad,
    const std::vector<double>& hour_angles_rad,
    double                     declination_rad,
    double                     wavelength_m);

}  // namespace argus::interferometry
