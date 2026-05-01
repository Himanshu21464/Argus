// Tests for the Rayleigh-scattering opacity kernel.

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

  // 1. Cross-section at 1 μm equals the reference value exactly.
  {
    RayleighOpacity ray("H2", /*sigma_1um=*/8.49e-29);
    // 1 μm = 10000 cm⁻¹
    Tensor s = ray.cross_section({10000.0}, {1500.0}, {1.0});
    assert(close(s[0], 8.49e-29, 1.0e-12));
  }

  // 2. λ⁻⁴ scaling: σ(0.5 μm) = 16 · σ(1 μm); σ(2 μm) = σ(1μm)/16.
  {
    RayleighOpacity ray("H2", 8.49e-29);
    Tensor s = ray.cross_section({20000.0, 10000.0, 5000.0},
                                 {1500.0}, {1.0});
    // 20000 cm⁻¹ = 0.5 μm  -> sigma = sigma_1um * 16
    // 5000 cm⁻¹ = 2 μm     -> sigma = sigma_1um / 16
    assert(close(s[0], 8.49e-29 * 16.0, 1.0e-9));
    assert(close(s[1], 8.49e-29,        1.0e-12));
    assert(close(s[2], 8.49e-29 / 16.0, 1.0e-9));
  }

  // 3. Constant in T and P (Rayleigh has no T,P dependence).
  {
    RayleighOpacity ray("H2", 8.49e-29);
    Tensor a = ray.cross_section({10000.0}, {100.0}, {0.001});
    Tensor b = ray.cross_section({10000.0}, {3000.0}, {100.0});
    assert(close(a[0], b[0], 1.0e-12));
  }

  // 4. Bad input throws.
  bool threw = false;
  try { RayleighOpacity("X", -1.0); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);

  // 5. End-to-end Rayleigh slope: in a clear, line-free atmosphere with
  //    Rayleigh scattering on the bulk gas, transit depth at short λ
  //    must be greater than at long λ.
  {
    Atmosphere atm;
    atm.species = {{"H2", 2.016}};
    const std::size_t n = 60;
    atm.pressure_bar.resize(n);
    atm.temperature_k.assign(n, 1500.0);
    atm.mixing_ratios = Tensor({n, 1});
    const double log_top = std::log(1.0e-6);
    const double log_bot = std::log(1.0e2);
    for (std::size_t i = 0; i < n; ++i) {
      const double frac = static_cast<double>(i) /
                          static_cast<double>(n - 1);
      atm.pressure_bar[i] = std::exp(log_top + frac * (log_bot - log_top));
      atm.mixing_ratios.at(i, 0) = 1.0;          // bulk H2
    }
    atm.validate();

    auto ray = std::make_shared<RayleighOpacity>("H2", 8.49e-29);
    TransmissionModel m;
    m.add_opacity(ray);

    // 0.5 μm vs 2.0 μm
    Spectrum s = m.forward(atm, std::vector<double>{20000.0, 5000.0});
    assert(s.values[0] > s.values[1]);

    // Rayleigh slope: depth_short / depth_long should be substantial,
    // since σ ratio is 256× across this range.
    // Transit-depth difference grows as ln(σ_R) — Rayleigh atmospheres
    // produce a slope of d(R_eff²)/d(ln λ) ≈ -4 H R_p in the optically
    // thin limit. Just verify the ordering is right and the difference
    // is large enough to exceed numerical noise.
    assert(s.values[0] - s.values[1] > 1.0e-8);
  }

  return 0;
}
