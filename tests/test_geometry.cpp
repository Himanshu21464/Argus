// Hard tests for the hydrostatic geometry pass.
// Asserts the chord formula matches analytic spheres to machine
// precision and that the layer altitudes match the analytic isothermal
// scale-height profile to 6 significant figures.

#include <cassert>
#include <cmath>
#include <stdexcept>

#include "argus/argus.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

}  // namespace

int main() {
  using namespace argus;

  // 1. Chord through the diameter (b = 0) of a sphere of radius r.
  for (double r : {1.0, 7.0e7, 6.957e8}) {
    assert(close(chord_path_length(0.0, r, 0.0), r, 1.0e-12));
  }

  // 2. Chord at b = r outer (tangent) -> 0 length.
  for (double r : {1.0, 7.0e7}) {
    assert(chord_path_length(0.5 * r, r, r) == 0.0);
    assert(chord_path_length(0.0,    r, r) == 0.0);
  }

  // 3. Chord at b = r_inner (skimming the inner shell) — analytic answer.
  //    Through a shell [r_in, r_out] at b = r_in: path = sqrt(r_out² - r_in²)
  for (double r_in : {1.0, 0.5, 7.0e7}) {
    const double r_out = 2.0 * r_in;
    const double expected = std::sqrt(r_out * r_out - r_in * r_in);
    assert(close(chord_path_length(r_in, r_out, r_in), expected, 1.0e-12));
  }

  // 4. Chord far above shell -> 0.
  assert(chord_path_length(1.0, 2.0, 3.0) == 0.0);

  // 5. Build hydrostatic geometry for an isothermal H2/He atmosphere and
  //    cross-check against the closed-form  z(P) = -H ln(P/P_bot).
  Species h2o{"H2O", 18.015};
  const double T = 1500.0;
  Atmosphere atm = isothermal(T, 1.0e-7, 1.0e2, 80, h2o, 1.0e-3);
  // hot-Jupiter defaults: planet_radius_rj=1, gravity=25 m/s², mmw=2.3
  Geometry g = build_geometry(atm);

  const double mu_kg = atm.bulk_mmw_amu * kAmuKg;
  const double H_m   = kBoltzmannSI * T / (mu_kg * atm.planet_gravity_si);
  const double R_p   = g.planet_radius_m;
  const double P_bot = atm.pressure_bar.back();

  for (std::size_t i = 0; i < g.radius_m.size(); ++i) {
    const double z_actual = g.radius_m[i] - R_p;
    const double z_expect = -H_m * std::log(atm.pressure_bar[i] / P_bot);
    // The integration path uses average T between adjacent layers; for an
    // isothermal atmosphere this is exact, so tolerance is tight.
    assert(close(z_actual, z_expect, 1.0e-9, 1.0e-3));
  }

  // 6. Anchor: bottom layer sits exactly at the planet radius.
  assert(close(g.radius_m.back(), R_p, 1.0e-12));

  // 7. Top of atmosphere is at altitude > 0.
  assert(g.radius_m.front() > R_p);

  // 8. Total atmosphere height ≈ -H * ln(P_top / P_bot).
  const double h_total_expected = -H_m *
      std::log(atm.pressure_bar.front() / P_bot);
  const double h_total_actual = g.radius_m.front() - R_p;
  assert(close(h_total_actual, h_total_expected, 1.0e-9, 1.0e-3));

  // 9. build_geometry rejects an atmosphere with non-monotonic
  //    pressures (interior misordering, not just front-vs-back).
  //    A boundary-only sort check would let this through and produce
  //    negative scale-height contributions in the integration loop.
  {
    Atmosphere bad;
    bad.species = {{"H2O", 18.015}};
    bad.pressure_bar = {1.0e-3, 1.0e-1, 1.0e-2, 1.0};   // out of order
    bad.temperature_k.assign(4, 1500.0);
    bad.mixing_ratios = Tensor({4, 1});
    for (std::size_t i = 0; i < 4; ++i) bad.mixing_ratios.at(i, 0) = 1.0e-3;

    bool threw = false;
    try { (void)build_geometry(bad); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // 10. build_geometry rejects an atmosphere with a non-positive
  //     pressure (would NaN the log() in the scale-height integral).
  {
    Atmosphere bad;
    bad.species = {{"H2O", 18.015}};
    bad.pressure_bar = {0.0, 1.0e-2, 1.0};               // zero at top
    bad.temperature_k.assign(3, 1500.0);
    bad.mixing_ratios = Tensor({3, 1});
    for (std::size_t i = 0; i < 3; ++i) bad.mixing_ratios.at(i, 0) = 1.0e-3;

    bool threw = false;
    try { (void)build_geometry(bad); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  return 0;
}
