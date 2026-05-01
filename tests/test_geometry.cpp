#include <cassert>
#include <cmath>

#include "argus/argus.hpp"

int main() {
  using namespace argus;

  // 1. Build an isothermal atmosphere and inspect its hydrostatic geometry.
  Species h2o{"H2O", 18.015};
  Atmosphere atm = isothermal(/*T=*/1200.0,
                              /*P_top=*/1.0e-6,
                              /*P_bot=*/1.0e2,
                              /*n_layers=*/40,
                              h2o,
                              /*VMR=*/1.0e-3);

  Geometry g = build_geometry(atm);

  // Top-of-atmosphere radius is larger than planet radius.
  assert(g.radius_m[0] > g.planet_radius_m);
  // Bottom-of-atmosphere radius equals planet radius (z=0 anchor).
  assert(std::fabs(g.radius_m.back() - g.planet_radius_m) < 1e-6);

  // Radii are monotonically decreasing from top to bottom.
  for (std::size_t i = 1; i < g.radius_m.size(); ++i) {
    assert(g.radius_m[i] < g.radius_m[i - 1]);
  }

  // For an isothermal atmosphere the altitude difference between adjacent
  // layers should match the scale height times the log-pressure increment:
  //   dz = H * dlnP
  // with H = k_B T / (mu g).
  const double mu_kg = atm.bulk_mmw_amu * kAmuKg;
  const double H_m   = kBoltzmannSI * atm.temperature_k[0] /
                       (mu_kg * atm.planet_gravity_si);
  const double dlnP_expected =
      std::log(atm.pressure_bar[1] / atm.pressure_bar[0]);
  const double dz_expected = -H_m * dlnP_expected;  // positive (going up)
  const double dz_actual   = g.radius_m[0] - g.radius_m[1];
  assert(std::fabs(dz_actual - dz_expected) / dz_expected < 1e-3);

  // 2. Chord path length sanity.
  // Chord at b = r_inner is tangent to the inner shell, so the path through
  // the shell is just the outer arc.
  const double ds = chord_path_length(/*r_in=*/1.0, /*r_out=*/2.0, /*b=*/1.0);
  assert(std::fabs(ds - std::sqrt(3.0)) < 1e-12);

  // Chord above the shell -> zero length.
  assert(chord_path_length(1.0, 2.0, 3.0) == 0.0);

  // Chord through the centre (b=0) -> full thickness.
  const double ds0 = chord_path_length(1.0, 2.0, 0.0);
  assert(std::fabs(ds0 - 1.0) < 1e-12);

  return 0;
}
