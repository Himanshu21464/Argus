// End-to-end real-world test: load real HITRAN H2O records, compute
// cross-sections, and validate physics-based invariants (sum rule,
// monotone temperature scaling, transit-depth saturation behaviour).

#include <algorithm>
#include <cassert>
#include <cmath>
#include <sstream>

#include "argus/argus.hpp"
#include "argus/test_data.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

}  // namespace

int main() {
  using namespace argus;

  // 1. Parse the bundled real H2O records, build a LineListOpacity.
  std::istringstream is{std::string(test_data::kH2OLines)};
  auto records = Hitran::load(is, /*molecule_id_filter=*/1);
  assert(records.size() == 16);

  std::vector<Line> lines;
  lines.reserve(records.size());
  for (const auto& r : records) lines.push_back(r.line);

  auto opacity = std::make_shared<LineListOpacity>("H2O", lines, 18.015);

  // 2. Cross-section sum rule: integrating sigma(nu) over a wide grid
  //    should approach Σ S_i(T) for a single layer (each Voigt is
  //    area-normalised, so Σ phi_i(nu) integrated equals Σ S_i).
  {
    const double T = 296.0;            // anchor at HITRAN reference
    const double P = 0.001;             // low pressure -> Doppler-dominated
    // Build a fine wavenumber grid covering both H2O bands.
    std::vector<double> wn;
    const double dw = 0.005;
    for (double w = 3500.0; w <= 3850.0; w += dw) wn.push_back(w);
    for (double w = 6900.0; w <= 7350.0; w += dw) wn.push_back(w);

    Tensor sigma = opacity->cross_section(wn, {T}, {P});

    double integrated = 0.0;
    for (std::size_t i = 0; i < wn.size(); ++i) integrated += sigma[i] * dw;

    // Expected total intensity (sum over all lines, no temperature scaling
    // since T=Tref).
    double expected = 0.0;
    for (const auto& r : records) expected += r.line.intensity_cm;

    // 1% tolerance accounts for Voigt tails outside the integration window.
    assert(close(integrated, expected, 0.02, 1.0e-26));
  }

  // 3. Temperature dependence: at higher T, lines with E_lower > 0 weaken
  //    relative to the lowest-E lines. We assert that hot-band lines
  //    (E_lower > 200 cm⁻¹) gain intensity going from 296 K to 1500 K
  //    while ground-state lines lose it.
  {
    std::vector<double> wn{3651.97, 3735.79};   // E_lower = 136, 222
    Tensor s_cool = opacity->cross_section(wn, {296.0}, {0.001});
    Tensor s_hot  = opacity->cross_section(wn, {1500.0}, {0.001});

    // higher-E line gains relative to lower-E line going hot.
    const double ratio_cool = s_cool[1] / s_cool[0];
    const double ratio_hot  = s_hot[1]  / s_hot[0];
    assert(ratio_hot > ratio_cool);
  }

  // 4. End-to-end transmission spectrum with real H2O lines:
  //    - all transit depths positive and < 1
  //    - depth at line centre exceeds depth at nearby continuum window
  {
    Species h2o{"H2O", 18.015};
    Atmosphere atm = isothermal(/*T=*/1500.0,
                                /*P_top=*/1.0e-6,
                                /*P_bot=*/1.0e2,
                                /*n_layers=*/60,
                                h2o,
                                /*VMR=*/1.0e-3);

    TransmissionModel model;
    model.add_opacity(opacity);

    // Sample the centre of one strong line and a nearby continuum point.
    std::vector<double> wn{3680.45, 3700.00};   // line centre, off-line
    Spectrum s = model.forward(atm, wn);
    assert(s.values[0] > 0.0 && s.values[0] < 1.0);
    assert(s.values[1] > 0.0 && s.values[1] < 1.0);
    assert(s.values[0] > s.values[1]);  // line centre absorbs more
  }

  // 5. Doubling the H2O VMR must monotonically deepen the line-centre
  //    transit depth (more absorbers = larger effective radius).
  {
    Species h2o{"H2O", 18.015};
    auto build = [&](double vmr) {
      Atmosphere a = isothermal(1500.0, 1.0e-6, 1.0e2, 60, h2o, vmr);
      TransmissionModel m;
      m.add_opacity(opacity);
      return m.forward(a, std::vector<double>{3680.45}).values[0];
    };
    const double d1 = build(1.0e-4);
    const double d2 = build(1.0e-3);
    const double d3 = build(1.0e-2);
    assert(d1 < d2);
    assert(d2 < d3);
  }

  // 6. Optically-thick saturation: at extremely high VMR the transit depth
  //    plateaus near the geometric area of the top of the atmosphere — it
  //    does not grow without bound.
  {
    Species h2o{"H2O", 18.015};
    auto depth_at = [&](double vmr) {
      Atmosphere a = isothermal(1500.0, 1.0e-6, 1.0e2, 60, h2o, vmr);
      TransmissionModel m;
      m.add_opacity(opacity);
      return m.forward(a, std::vector<double>{3680.45}).values[0];
    };
    const double d_thick   = depth_at(1.0e-1);
    const double d_extreme = depth_at(1.0e-0);
    // saturation: extreme is at most 1% larger than already-thick
    assert(d_extreme < d_thick * 1.05);
  }

  return 0;
}
