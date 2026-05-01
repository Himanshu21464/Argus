// Comprehensive end-to-end test: a realistic hot Jupiter retrieval setup
// with EVERYTHING active — Guillot T-P profile, H2O + CO2 + CH4 line
// opacity (real HITRAN), gray cloud deck, Rayleigh scattering on the
// bulk H2/He.
//
// Asserts the qualitative features a real WASP-39b-class spectrum
// shows: Rayleigh slope at short λ, distinct molecular bands sticking
// up above continuum, cloud floor capping the deepest features.

#include <cassert>
#include <cmath>
#include <sstream>
#include <vector>

#include "argus/argus.hpp"
#include "argus/test_data.hpp"

int main() {
  using namespace argus;

  // ─── 1. Load all three molecule line lists from real HITRAN data. ──
  auto load_species = [](std::string_view raw, int filter) {
    std::istringstream is{std::string(raw)};
    auto records = Hitran::load(is, filter);
    std::vector<Line> lines;
    lines.reserve(records.size());
    for (const auto& r : records) lines.push_back(r.line);
    return lines;
  };

  auto h2o_lines = load_species(test_data::kH2OLines, /*mol=*/1);
  auto co2_lines = load_species(test_data::kCO2Lines, /*mol=*/2);
  auto ch4_lines = load_species(test_data::kCH4Lines, /*mol=*/6);
  assert(h2o_lines.size() == 16);
  assert(co2_lines.size() == 10);
  assert(ch4_lines.size() == 8);

  auto h2o_op = std::make_shared<LineListOpacity>("H2O", h2o_lines, 18.015);
  auto co2_op = std::make_shared<LineListOpacity>("CO2", co2_lines, 44.010);
  auto ch4_op = std::make_shared<LineListOpacity>("CH4", ch4_lines, 16.043);
  auto cld_op = std::make_shared<CloudDeckOpacity>(
      "CLOUD_DECK", /*P_cloud=*/1.0e-2, /*sigma=*/1.0e-18);
  auto ray_op = std::make_shared<RayleighOpacity>("H2_BULK", 8.49e-29);

  // ─── 2. Build a Guillot 2010 hot-Jupiter atmosphere with 5 species. ──
  // Species: H2 (bulk + Rayleigh), H2O, CO2, CH4, CLOUD_DECK
  Atmosphere atm;
  const std::size_t n = 80;
  atm.species = {
    {"H2_BULK",    2.016},
    {"H2O",        18.015},
    {"CO2",        44.010},
    {"CH4",        16.043},
    {"CLOUD_DECK", 1.0},
  };
  atm.pressure_bar.resize(n);
  atm.mixing_ratios = Tensor({n, 5});

  const double log_top = std::log(1.0e-6);
  const double log_bot = std::log(1.0e2);
  for (std::size_t i = 0; i < n; ++i) {
    const double frac = static_cast<double>(i) /
                        static_cast<double>(n - 1);
    atm.pressure_bar[i] = std::exp(log_top + frac * (log_bot - log_top));
    atm.mixing_ratios.at(i, 0) = 1.0;        // H2 bulk
    atm.mixing_ratios.at(i, 1) = 1.0e-3;      // H2O VMR
    atm.mixing_ratios.at(i, 2) = 5.0e-4;      // CO2 VMR
    atm.mixing_ratios.at(i, 3) = 5.0e-5;      // CH4 VMR
    atm.mixing_ratios.at(i, 4) = 1.0;        // cloud anchor (kernel handles cutoff by P)
  }

  atm.temperature_k = guillot_profile(atm.pressure_bar,
                                      /*T_int=*/200.0,
                                      /*T_irr=*/1500.0,
                                      /*gamma=*/0.5);
  atm.validate();

  // ─── 3. Forward model with all kernels stacked. ──────────────────────
  TransmissionModel model;
  model.add_opacity(ray_op);
  model.add_opacity(h2o_op);
  model.add_opacity(co2_op);
  model.add_opacity(ch4_op);
  model.add_opacity(cld_op);

  // JWST-PRISM-like grid (0.5–5 μm).
  std::vector<double> wn;
  for (double w = 2000.0; w <= 20000.0; w += 100.0) wn.push_back(w);
  Spectrum s = model.forward(atm, wn);

  // ─── 4. ASSERTIONS — physics-based ───────────────────────────────────

  // 4a. Spectrum is well-formed: positive, finite, < 1 transit depth.
  for (double v : s.values) {
    assert(std::isfinite(v));
    assert(v > 0.0 && v < 1.0);
  }

  // 4b. Rayleigh slope: depth at 0.5 μm (20000 cm⁻¹) should be
  //     significantly greater than depth at the continuum window
  //     between bands at ~1.7 μm (~5882 cm⁻¹).
  auto idx = [&](double w_target) {
    std::size_t best = 0;
    double best_dx = 1.0e30;
    for (std::size_t i = 0; i < wn.size(); ++i) {
      const double dx = std::fabs(wn[i] - w_target);
      if (dx < best_dx) { best_dx = dx; best = i; }
    }
    return best;
  };
  const std::size_t i_blue = idx(20000.0);   // 0.5 μm
  const std::size_t i_cont = idx(5800.0);    // ~1.7 μm continuum window
  assert(s.values[i_blue] > s.values[i_cont]);

  // 4c. H2O 1.4 μm band (~7100 cm⁻¹) absorbs more than the continuum
  //     either side (5800 cm⁻¹ and 9000 cm⁻¹).
  const std::size_t i_h2o14 = idx(7100.0);
  const std::size_t i_cont_blue = idx(9000.0);
  assert(s.values[i_h2o14] > s.values[i_cont]);
  assert(s.values[i_h2o14] > s.values[i_cont_blue]);

  // 4d. CO2 4.3 μm band (~2363 cm⁻¹) absorbs more than the nearby
  //     continuum at ~2200 cm⁻¹.
  const std::size_t i_co2 = idx(2363.0);
  const std::size_t i_co2_cont = idx(2200.0);
  assert(s.values[i_co2] > s.values[i_co2_cont]);

  // 4e. CH4 3.3 μm band (~3000 cm⁻¹) absorbs more than the continuum.
  const std::size_t i_ch4 = idx(3000.0);
  const std::size_t i_ch4_cont = idx(2700.0);
  assert(s.values[i_ch4] > s.values[i_ch4_cont]);

  // 4f. Cloud floor: the deepest transit depth in the spectrum should
  //     not exceed a "cloudy worst case" upper bound — clouds cap the
  //     dynamic range. Without a cloud the spectrum can have arbitrarily
  //     deep features; with the deck at P_cloud=10 mbar most absorption
  //     comes from above the cloud and is bounded.
  double max_depth = 0.0;
  for (double v : s.values) max_depth = std::max(max_depth, v);
  // For a hot Jupiter on sun: max depth ~ 1.5%; with our setup expect
  // < 5% (loose bound).
  assert(max_depth < 0.05);

  // 4g. Reproducibility — running twice gives byte-equal output.
  Spectrum s2 = model.forward(atm, wn);
  assert(s.values.size() == s2.values.size());
  for (std::size_t i = 0; i < s.values.size(); ++i) {
    assert(s.values[i] == s2.values[i]);
  }

  return 0;
}
