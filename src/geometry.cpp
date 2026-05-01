#include "argus/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace argus {

Geometry build_geometry(const Atmosphere& atm) {
  atm.validate();
  const std::size_t n = atm.num_layers();
  if (n < 2) {
    throw std::invalid_argument(
        "build_geometry: need at least 2 layers");
  }

  Geometry g;
  g.planet_radius_m = atm.planet_radius_rj * kJupiterRadiusM;
  g.radius_m.assign(n, 0.0);
  g.thickness_m.assign(n, 0.0);

  // Convention used throughout the kernel:
  //   layer 0 is the TOP of the atmosphere (lowest pressure)
  //   layer n-1 is the BOTTOM (highest pressure, sits at planet_radius_m)
  // Pressures must be sorted so this invariant holds.
  const bool top_first = atm.pressure_bar.front() < atm.pressure_bar.back();
  if (!top_first) {
    throw std::invalid_argument(
        "build_geometry: atmosphere layers must be ordered with the lowest "
        "pressure first (top of atmosphere = layer 0)");
  }

  // Walk from the bottom (z = 0) upward, accumulating altitude using the
  // local scale height between adjacent layers.
  const double mu_kg = atm.bulk_mmw_amu * kAmuKg;
  const double g_si = atm.planet_gravity_si;

  std::vector<double> z(n, 0.0);
  z[n - 1] = 0.0;
  for (std::size_t i = n - 1; i > 0; --i) {
    // Use the average temperature between layers i and i-1 for the scale
    // height in that interval.
    const double T_avg = 0.5 * (atm.temperature_k[i] + atm.temperature_k[i - 1]);
    const double H_m = kBoltzmannSI * T_avg / (mu_kg * g_si);
    const double dlnP = std::log(atm.pressure_bar[i] / atm.pressure_bar[i - 1]);
    // dlnP > 0 because layer i is below layer i-1 (higher pressure).
    z[i - 1] = z[i] + H_m * dlnP;
  }

  for (std::size_t i = 0; i < n; ++i) {
    g.radius_m[i] = g.planet_radius_m + z[i];
  }

  // Per-layer radial thickness. For interior layers, thickness is the
  // half-distance to each neighbour. For the boundary layers, just use the
  // distance to the single neighbour.
  for (std::size_t i = 0; i < n; ++i) {
    double dr_above = (i > 0)        ? std::fabs(g.radius_m[i] - g.radius_m[i - 1]) : 0.0;
    double dr_below = (i + 1 < n)    ? std::fabs(g.radius_m[i + 1] - g.radius_m[i]) : 0.0;
    if (i == 0)         g.thickness_m[i] = dr_below;
    else if (i + 1 == n) g.thickness_m[i] = dr_above;
    else                 g.thickness_m[i] = 0.5 * (dr_above + dr_below);
  }

  return g;
}

double chord_path_length(double r_inner_m, double r_outer_m, double b_m) {
  if (r_outer_m <= b_m) return 0.0;          // chord is above this shell
  const double outer = std::sqrt(r_outer_m * r_outer_m - b_m * b_m);
  const double inner = (r_inner_m > b_m)
                       ? std::sqrt(r_inner_m * r_inner_m - b_m * b_m)
                       : 0.0;
  return outer - inner;
}

}  // namespace argus
