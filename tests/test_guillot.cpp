// Tests for the Guillot 2010 analytic T-P profile builder.

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
  Species h2o{"H2O", 18.015};

  // 1. Build a typical hot-Jupiter T-P profile.
  Atmosphere atm = guillot(/*T_int=*/200.0, /*T_irr=*/1500.0,
                           /*gamma=*/0.5,
                           /*P_top=*/1.0e-6, /*P_bot=*/1.0e2,
                           /*n_layers=*/60,
                           h2o, /*VMR=*/1.0e-3);
  atm.validate();

  assert(atm.num_layers() == 60);
  // 2. All temperatures finite & positive, in plausible range.
  for (double T : atm.temperature_k) {
    assert(std::isfinite(T));
    assert(T > 100.0 && T < 5000.0);
  }

  // 3. The deep atmosphere (high P) should be the hottest layer thanks
  //    to the internal-heat term scaling with τ ∝ P.
  const double T_top  = atm.temperature_k.front();
  const double T_bot  = atm.temperature_k.back();
  assert(T_bot > T_top);

  // 4. T at very low P (τ → 0) is the "skin temperature" of an
  //    irradiated atmosphere:
  //      T_skin⁴ ≈ (3/4) T_int⁴ · (2/3) + (3/4) T_irr⁴ · [2/3 + 1/(γ√3) + (γ/√3 - 1/(γ√3))]
  //             = (3/4) T_int⁴ · (2/3) + (3/4) T_irr⁴ · [2/3 + γ/√3]
  // For T_int=200, T_irr=1500, γ=0.5:
  //   T_skin⁴ = 0.75 · 200⁴ · 0.667 + 0.75 · 1500⁴ · (0.667 + 0.289)
  //           = 800000 + 3.629e12 · ~0.956 = 3.469e12
  //   T_skin ≈ 1287
  {
    const double sqrt3 = std::sqrt(3.0);
    const double T_int4 = std::pow(200.0, 4.0);
    const double T_irr4 = std::pow(1500.0, 4.0);
    const double a = 2.0 / 3.0;
    const double b = 1.0 / (0.5 * sqrt3);
    const double c = (0.5 / sqrt3) - b;
    const double T4_skin = 0.75 * T_int4 * a
                         + 0.75 * T_irr4 * (a + b + c);  // exp(0) = 1
    const double T_skin = std::pow(T4_skin, 0.25);
    assert(close(T_top, T_skin, 1.0e-3));
  }

  // 5. Profile is monotonically non-decreasing top -> bottom (the
  //    Guillot profile may show a stratosphere inversion; for γ < 1
  //    this is monotone).
  for (std::size_t i = 1; i < atm.temperature_k.size(); ++i) {
    assert(atm.temperature_k[i] >= atm.temperature_k[i - 1] - 1.0);
  }

  // 6. Geometry on Guillot atmosphere works (non-isothermal hydrostatic).
  Geometry g = build_geometry(atm);
  assert(g.radius_m.size() == atm.num_layers());
  assert(g.radius_m.front() > g.radius_m.back());

  // 7. Forward model on Guillot atmosphere produces sensible spectrum.
  auto opacity = std::make_shared<GreyOpacity>("H2O", 1.0e-22);
  TransmissionModel m;
  m.add_opacity(opacity);
  Spectrum s = m.forward(atm, std::vector<double>{4000.0});
  assert(s.values[0] > 0.0 && s.values[0] < 1.0);

  // 8. Higher T_irr -> hotter atmosphere -> larger scale height -> deeper transit.
  Atmosphere a_cold = guillot(200.0, 1000.0, 0.5, 1.0e-6, 1.0e2, 60, h2o, 1.0e-3);
  Atmosphere a_hot  = guillot(200.0, 2500.0, 0.5, 1.0e-6, 1.0e2, 60, h2o, 1.0e-3);
  Spectrum sc = m.forward(a_cold, std::vector<double>{4000.0});
  Spectrum sh = m.forward(a_hot,  std::vector<double>{4000.0});
  assert(sh.values[0] > sc.values[0]);

  // 9. Bad inputs throw.
  bool threw = false;
  try { (void)guillot(200.0, 1500.0, /*gamma=*/-0.5, 1e-6, 1e2, 60, h2o, 1e-3); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);

  threw = false;
  try { (void)guillot(200.0, 1500.0, 0.5, /*P_top=*/1.0, /*P_bot=*/0.5, 60, h2o, 1e-3); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);

  return 0;
}
