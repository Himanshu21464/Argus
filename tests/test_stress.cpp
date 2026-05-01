// Stress + edge-case tests. Production-grade kernels must not crash, hang,
// or silently produce garbage on degenerate inputs. This file pokes the
// boundary conditions hard.

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "argus/argus.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

}  // namespace

int main() {
  using namespace argus;
  Species h2o{"H2O", 18.015};

  // 1. Atmosphere validation: malformed inputs throw, not segfault.
  {
    Atmosphere a;
    a.pressure_bar = {1.0e-3, 1.0};
    a.temperature_k = {1500.0};                // mismatched length
    a.species = {h2o};
    a.mixing_ratios = Tensor({2, 1});
    bool threw = false;
    try { a.validate(); } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // 2. Mixing ratio outside [0,1] is rejected.
  {
    Atmosphere a;
    a.pressure_bar = {1.0e-3, 1.0};
    a.temperature_k = {1500.0, 1500.0};
    a.species = {h2o};
    a.mixing_ratios = Tensor({2, 1}, {0.5, 1.5});       // 1.5 invalid
    bool threw = false;
    try { a.validate(); } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // 3. isothermal builder rejects degenerate inputs.
  {
    bool threw = false;
    try { (void)isothermal(1500.0, 1.0, 0.5, 10, h2o, 1e-3); }   // P_top > P_bot
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    threw = false;
    try { (void)isothermal(1500.0, 1e-6, 1.0, 1, h2o, 1e-3); }   // n_layers=1
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // 4. Geometry rejects 1-layer atmosphere.
  {
    Atmosphere a;
    a.pressure_bar = {1.0};
    a.temperature_k = {1500.0};
    a.species = {h2o};
    a.mixing_ratios = Tensor({1, 1}, {1e-3});
    bool threw = false;
    try { (void)build_geometry(a); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // 5. Geometry rejects mis-ordered (top-first invariant violated).
  {
    Atmosphere a;
    a.pressure_bar = {1.0e2, 1.0e-6};   // bottom first — wrong order
    a.temperature_k = {1500.0, 1500.0};
    a.species = {h2o};
    a.mixing_ratios = Tensor({2, 1}, {1e-3, 1e-3});
    bool threw = false;
    try { (void)build_geometry(a); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // 6. Tensor shape mismatch in constructor throws.
  {
    bool threw = false;
    try { Tensor t({2, 3}, {1.0, 2.0}); }   // 6 elements expected, 2 given
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // 7. TransmissionModel with zero opacity returns the geometric baseline.
  //    Optical depth = 0 everywhere => transit depth = (R_p / R_*)^2.
  {
    Atmosphere atm = isothermal(1500.0, 1.0e-6, 1.0e2, 30, h2o, 1.0e-3);
    TransmissionModel m;
    Spectrum s = m.forward(atm, std::vector<double>{1000.0, 5000.0});
    const double R_p = atm.planet_radius_rj * kJupiterRadiusM;
    const double R_s = atm.star_radius_rsun * kSolarRadiusM;
    const double expected = (R_p / R_s) * (R_p / R_s);
    assert(close(s.values[0], expected, 1.0e-9));
    assert(close(s.values[1], expected, 1.0e-9));
  }

  // 8. add_opacity rejects null pointer.
  {
    TransmissionModel m;
    bool threw = false;
    try { m.add_opacity(std::shared_ptr<OpacityKernel>{}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // 9. LineListOpacity rejects non-positive molar mass.
  {
    bool threw = false;
    try { LineListOpacity("X", {}, 0.0); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    threw = false;
    try { LineListOpacity("X", {}, -1.0); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // 10. Voigt is finite + positive at extreme conditions.
  //     T=100 K (cold), T=3000 K (very hot), tiny gamma, huge gamma.
  {
    const double conditions[][2] = {
      {0.001, 0.001},     // both narrow
      {0.001, 1.0},       // narrow Doppler, wide Lorentz
      {1.0, 0.001},       // wide Doppler, narrow Lorentz (HAW degraded; Gaussian fallback)
      {1.0, 1.0},         // both wide
      {0.05, 100.0},      // huge Lorentz
    };
    for (const auto& c : conditions) {
      const double V = voigt<double>(0.0, c[0], c[1]);
      assert(std::isfinite(V));
      assert(V > 0.0);
    }
  }

  // 11. LineListOpacity at extreme T (50K, 4000K) does not throw and gives
  //     positive cross-sections.
  {
    std::vector<Line> lines{
      Line{3000.0, 1e-21, 0.05, 0.0, 0.5, 0.0, 100.0},
    };
    LineListOpacity lo("H2O", lines, 18.015);
    Tensor cold = lo.cross_section({3000.0}, {50.0}, {1.0});
    Tensor hot  = lo.cross_section({3000.0}, {4000.0}, {1.0});
    assert(std::isfinite(cold[0]) && cold[0] > 0.0);
    assert(std::isfinite(hot[0])  && hot[0]  > 0.0);
  }

  // 12. Single-wavelength forward model.
  {
    Atmosphere atm = isothermal(1500.0, 1.0e-6, 1.0e2, 20, h2o, 1.0e-3);
    auto op = std::make_shared<GreyOpacity>("H2O", 1.0e-22);
    TransmissionModel m;
    m.add_opacity(op);
    Spectrum s = m.forward(atm, std::vector<double>{2000.0});
    assert(s.values.size() == 1);
    assert(s.values[0] > 0.0 && s.values[0] < 1.0);
  }

  // 13. Empty wavenumber grid returns empty spectrum, no crash.
  {
    Atmosphere atm = isothermal(1500.0, 1.0e-6, 1.0e2, 20, h2o, 1.0e-3);
    TransmissionModel m;
    Spectrum s = m.forward(atm, std::vector<double>{});
    assert(s.values.empty());
    assert(s.wavenumber_cm.empty());
  }

  // 14. Empty line list returns zero cross-section.
  {
    LineListOpacity lo("H2O", {}, 18.015);
    Tensor t = lo.cross_section({3000.0, 4000.0}, {1500.0}, {1.0});
    assert(t[0] == 0.0);
    assert(t[1] == 0.0);
  }

  // 15. Hitran::parse_line on completely garbage input returns nullopt.
  {
    assert(!Hitran::parse_line("").has_value());
    assert(!Hitran::parse_line("garbage").has_value());
    assert(!Hitran::parse_line("    ").has_value());
  }

  return 0;
}
