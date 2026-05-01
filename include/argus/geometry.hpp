#pragma once

#include <vector>

#include "atmosphere.hpp"

namespace argus {

// Hydrostatic geometry of a 1-D layered atmosphere.
//
// The kernel models the planet's limb as a stack of concentric spherical
// shells. Each `Atmosphere` layer becomes one shell, with the layer centre
// at radius `radius_m[i]`. M2 ships the proper geometric chord through
// these shells; the M1 stub used a flat 1-km path per layer.
struct Geometry {
  // Planet "solid" reference radius (m) — corresponds to the bottom layer's
  // pressure (the highest pressure in the atmosphere).
  double planet_radius_m = 0.0;
  // [n_layers] — radius (m) at the centre of each layer.
  std::vector<double> radius_m;
  // [n_layers] — radial thickness (m) of each layer (used for impact-parameter
  // integration weights in the transit-radius formula).
  std::vector<double> thickness_m;
};

// Build the geometry of `atm` from hydrostatic equilibrium.
// Bottom layer (highest pressure) sits at planet_radius_m.
// Each layer's altitude follows from the local scale height
//     H = k_B * T / (mu * g)
// integrated layer by layer.
Geometry build_geometry(const Atmosphere& atm);

// Path length (m) of a chord at impact parameter `b_m` through one
// spherical shell bounded by [r_inner_m, r_outer_m]. Returns 0 if the chord
// does not enter the shell. Single-pass length (the full chord traverses
// each shell twice — caller multiplies by 2 if needed).
double chord_path_length(double r_inner_m, double r_outer_m, double b_m);

// Useful constants exposed for tests.
inline constexpr double kJupiterRadiusM = 7.1492e7;     // R_J in metres
inline constexpr double kSolarRadiusM   = 6.957e8;      // R_sun in metres
inline constexpr double kBoltzmannSI    = 1.380649e-23; // J / K
inline constexpr double kAmuKg          = 1.66053906660e-27;

}  // namespace argus
