#include <cassert>
#include <cmath>

#include "argus/argus.hpp"

int main() {
  using namespace argus;

  Species h2o{"H2O", 18.015};
  Atmosphere a = isothermal(1200.0, 1e-6, 1e2, 40, h2o, 1e-3);
  a.validate();

  assert(a.num_layers() == 40);
  assert(a.num_species() == 1);
  assert(a.species[0].key == "H2O");

  // Pressure should grow monotonically.
  for (std::size_t i = 1; i < a.pressure_bar.size(); ++i) {
    assert(a.pressure_bar[i] > a.pressure_bar[i - 1]);
  }
  // Top and bottom pressures should match what we asked for.
  assert(std::abs(a.pressure_bar.front() - 1e-6) < 1e-12);
  assert(std::abs(a.pressure_bar.back() - 1e2) < 1e-9);

  // All temperatures should be the isothermal value.
  for (double T : a.temperature_k) {
    assert(std::abs(T - 1200.0) < 1e-9);
  }

  return 0;
}
