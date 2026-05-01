#include <cassert>
#include <cmath>

#include "argus/argus.hpp"

int main() {
  using namespace argus;

  // Build a tiny H2O-like line list with two lines.
  std::vector<Line> lines;
  lines.push_back(Line{
    /*nu0=*/3000.0,  /*S=*/1.0e-21, /*g_air=*/0.05,
    /*g_self=*/0.0,   /*n_air=*/0.5, /*delta=*/0.0, /*E_low=*/0.0});
  lines.push_back(Line{
    /*nu0=*/3500.0,  /*S=*/5.0e-22, /*g_air=*/0.05,
    /*g_self=*/0.0,   /*n_air=*/0.5, /*delta=*/0.0, /*E_low=*/0.0});

  LineListOpacity h2o("H2O", lines, 18.015);
  assert(h2o.species_key() == "H2O");
  assert(h2o.num_lines() == 2);

  // Sample exactly at the line centres and at a midpoint. Line widths are
  // ~0.05 cm^-1 so a coarse grid would land in the wings.
  //   wn[0] = 3000 cm^-1 (centre of stronger line)
  //   wn[1] = 3500 cm^-1 (centre of weaker line)
  //   wn[2] = 3250 cm^-1 (midpoint between lines, far from either)
  std::vector<double> wn{3000.0, 3500.0, 3250.0};
  std::vector<double> T{1200.0};
  std::vector<double> P{1.0};

  Tensor sigma = h2o.cross_section(wn, T, P);
  assert(sigma.shape().size() == 3);
  assert(sigma.shape()[0] == 1);
  assert(sigma.shape()[1] == 1);
  assert(sigma.shape()[2] == wn.size());

  // 1. Cross-section is positive everywhere.
  for (std::size_t i = 0; i < sigma.size(); ++i) {
    assert(sigma[i] >= 0.0);
  }

  // 2. Each line peak is much larger than the wing midpoint.
  assert(sigma[0] > 1e3 * sigma[2]);
  assert(sigma[1] > 1e3 * sigma[2]);

  // 3. Stronger line (3000) is taller than weaker line (3500).
  assert(sigma[0] > sigma[1]);

  return 0;
}
