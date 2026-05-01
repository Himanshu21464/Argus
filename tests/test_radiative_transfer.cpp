#include <cassert>
#include <cmath>

#include "argus/argus.hpp"

int main() {
  using namespace argus;

  Species h2o{"H2O", 18.015};
  Atmosphere atm = isothermal(1200.0, 1e-6, 1e2, 40, h2o, 1e-3);

  TransmissionModel model;
  model.add_opacity(std::make_shared<GreyOpacity>("H2O", 1e-22));

  std::vector<double> wn{2000.0, 2500.0, 3000.0};
  Spectrum s = model.forward(atm, wn);

  // Spectrum has the right shape.
  assert(s.wavenumber_cm.size() == wn.size());
  assert(s.values.size() == wn.size());

  // Grey opacity → spectrum should be flat and in (0, 1).
  for (std::size_t i = 0; i < s.values.size(); ++i) {
    assert(s.values[i] > 0.0);
    assert(s.values[i] < 1.0);
    if (i > 0) {
      assert(std::abs(s.values[i] - s.values[0]) < 1e-9);
    }
  }
  return 0;
}
