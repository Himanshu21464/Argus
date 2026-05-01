// Tests for the gray cloud-deck opacity model.

#include <cassert>
#include <cmath>
#include <stdexcept>

#include "argus/argus.hpp"

int main() {
  using namespace argus;

  // 1. Below cloud-top pressure (low P): zero cross-section.
  // 2. At/above cloud-top pressure: large cross-section.
  {
    CloudDeckOpacity cloud("CLOUD_DECK", /*P_cloud=*/1.0e-2, /*sigma=*/1.0e-18);
    Tensor s = cloud.cross_section({1000.0, 5000.0},
                                   {1500.0},
                                   {1.0e-4, 1.0e-2, 1.0});
    // shape: [1, 3, 2]
    assert(s.shape()[0] == 1 && s.shape()[1] == 3 && s.shape()[2] == 2);

    // At P=1e-4 (above cloud, low pressure): zero.
    assert(s[0 * 3 * 2 + 0 * 2 + 0] == 0.0);
    assert(s[0 * 3 * 2 + 0 * 2 + 1] == 0.0);
    // At P=1e-2 (== P_cloud): opaque.
    assert(s[0 * 3 * 2 + 1 * 2 + 0] == 1.0e-18);
    // At P=1.0 (well below cloud): opaque.
    assert(s[0 * 3 * 2 + 2 * 2 + 0] == 1.0e-18);
  }

  // 3. End-to-end: cloud raises the transit-depth floor by raising the
  //    effective opaque radius. Compare a clear vs cloudy hot-Jupiter
  //    transmission spectrum.
  {
    Species h2o{"H2O", 18.015};
    Species cloud_sp{"CLOUD_DECK", 18.0};

    auto build = [&](double P_cloud_bar) {
      Atmosphere a;
      a.species = {h2o, cloud_sp};
      const std::size_t n = 50;
      a.pressure_bar.resize(n);
      a.temperature_k.assign(n, 1500.0);
      a.mixing_ratios = Tensor({n, 2});
      const double log_top = std::log(1.0e-6);
      const double log_bot = std::log(1.0e2);
      for (std::size_t i = 0; i < n; ++i) {
        const double frac = static_cast<double>(i) /
                            static_cast<double>(n - 1);
        a.pressure_bar[i] = std::exp(log_top + frac * (log_bot - log_top));
        a.mixing_ratios.at(i, 0) = 1.0e-3;     // H2O VMR
        a.mixing_ratios.at(i, 1) = 1.0;         // cloud "VMR" anchor
      }
      a.validate();
      auto h2o_op   = std::make_shared<GreyOpacity>("H2O", 1.0e-22);
      auto cloud_op = std::make_shared<CloudDeckOpacity>(
          "CLOUD_DECK", P_cloud_bar, 1.0e-18);
      TransmissionModel m;
      m.add_opacity(h2o_op);
      m.add_opacity(cloud_op);
      return m.forward(a, std::vector<double>{4000.0}).values[0];
    };

    const double clear   = build(1.0e10);   // cloud "very deep" — no effect
    const double mid     = build(1.0e-2);   // cloud at 10 mbar
    const double high    = build(1.0e-4);   // cloud at 0.1 mbar — high deck

    // Higher cloud deck (lower P_cloud) -> more layers opaque -> deeper
    // transit. Verify monotonic.
    assert(clear < mid);
    assert(mid < high);
  }

  // 4. Bad inputs throw.
  bool threw = false;
  try { CloudDeckOpacity("X", -1.0, 1.0e-18); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);

  threw = false;
  try { CloudDeckOpacity("X", 1.0, -1.0); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);

  return 0;
}
