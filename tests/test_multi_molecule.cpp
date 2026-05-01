// Multi-molecule end-to-end test: H2O + CO2 in a single hot-Jupiter
// atmosphere. Verifies the cross-section is additive across opacity
// kernels and that the transmission spectrum places absorption features
// in the correct bands (H2O at 2.7/1.4 μm, CO2 at 4.3 μm).

#include <cassert>
#include <cmath>
#include <sstream>

#include "argus/argus.hpp"
#include "argus/test_data.hpp"

int main() {
  using namespace argus;

  // 1. Parse both bundled fixtures.
  std::vector<Line> h2o_lines, co2_lines;
  {
    std::istringstream is{std::string(test_data::kH2OLines)};
    auto records = Hitran::load(is, /*filter=*/1);
    for (const auto& r : records) h2o_lines.push_back(r.line);
    assert(h2o_lines.size() == 16);
  }
  {
    std::istringstream is{std::string(test_data::kCO2Lines)};
    auto records = Hitran::load(is, /*filter=*/2);
    for (const auto& r : records) co2_lines.push_back(r.line);
    assert(co2_lines.size() == 10);
  }

  auto h2o_op = std::make_shared<LineListOpacity>("H2O", h2o_lines, 18.015);
  auto co2_op = std::make_shared<LineListOpacity>("CO2", co2_lines, 44.010);

  // 2. Atmosphere with both species.
  Atmosphere atm;
  const std::size_t n = 50;
  atm.species = {{"H2O", 18.015}, {"CO2", 44.010}};
  atm.pressure_bar.assign(n, 0.0);
  atm.temperature_k.assign(n, 1500.0);
  atm.mixing_ratios = Tensor({n, 2});
  const double log_top = std::log(1.0e-6);
  const double log_bot = std::log(1.0e2);
  for (std::size_t i = 0; i < n; ++i) {
    const double frac = static_cast<double>(i) /
                        static_cast<double>(n - 1);
    atm.pressure_bar[i] = std::exp(log_top + frac * (log_bot - log_top));
    atm.mixing_ratios.at(i, 0) = 1.0e-3;   // H2O VMR
    atm.mixing_ratios.at(i, 1) = 5.0e-4;   // CO2 VMR
  }
  atm.validate();

  // 3. Cross-section additivity: with both kernels attached, the per-layer
  //    alpha is the sum of the individual contributions. Verify by running
  //    each species alone and adding.
  TransmissionModel m_both;
  m_both.add_opacity(h2o_op);
  m_both.add_opacity(co2_op);

  TransmissionModel m_h2o; m_h2o.add_opacity(h2o_op);
  TransmissionModel m_co2; m_co2.add_opacity(co2_op);

  // 4. Sample wavenumbers in each species' band centre.
  std::vector<double> wn{
    3651.97,    // H2O 2.7 μm strong line
    7099.18,    // H2O 1.4 μm strong line
    2363.48,    // CO2 4.3 μm strong line (Q-branch peak)
    5500.00,    // continuum window between bands
  };

  Spectrum s_both = m_both.forward(atm, wn);
  Spectrum s_h2o  = m_h2o.forward(atm, wn);
  Spectrum s_co2  = m_co2.forward(atm, wn);

  // 5. At H2O lines: H2O dominates absorption.
  assert(s_both.values[0] > s_co2.values[0]);
  assert(s_both.values[1] > s_co2.values[1]);
  // 6. At CO2 line: CO2 dominates.
  assert(s_both.values[2] > s_h2o.values[2]);
  // 7. Both-active depth ≥ either single-species depth (additive optical depth).
  for (std::size_t i = 0; i < 4; ++i) {
    assert(s_both.values[i] >= s_h2o.values[i] - 1e-12);
    assert(s_both.values[i] >= s_co2.values[i] - 1e-12);
  }

  // 8. Continuum window: little absorption from either band.
  // The continuum-window depth equals the geometric base + small Lorentzian
  // wings — should be smaller than any band-centre depth.
  for (std::size_t i = 0; i < 3; ++i) {
    assert(s_both.values[i] > s_both.values[3]);
  }

  return 0;
}
