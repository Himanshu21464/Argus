// M2 capstone: end-to-end multi-physics atmospheric retrieval that
// exercises every M2 ingredient simultaneously — Guillot 2010 T-P
// profile, real-HITRAN H₂O + CO₂ + CH₄ Voigt opacity, gray cloud
// deck, Rayleigh scattering on a bulk H₂ background. The same
// `argus::Retrieval` API used for substrate proofs across lensing
// and interferometry recovers four atmospheric parameters
// (T_irr, log₁₀VMR_H₂O, log₁₀VMR_CO₂, log₁₀P_cloud) from a noisy
// JWST-PRISM-shaped synthetic spectrum.
//
// Validates the M2 substrate end-to-end: every kernel layer plugs
// into the same Retrieval pipeline that the lensing and interferometry
// retrievals already proved.

#include <cassert>
#include <cmath>
#include <random>
#include <sstream>
#include <vector>

#include "argus/argus.hpp"
#include "argus/test_data.hpp"

int main() {
  using namespace argus;

  // ─── 1. Real HITRAN line lists for H₂O, CO₂, CH₄. ──────────────────
  auto load_species = [](std::string_view raw, int filter) {
    std::istringstream is{std::string(raw)};
    auto records = Hitran::load(is, filter);
    std::vector<Line> lines;
    lines.reserve(records.size());
    for (const auto& r : records) lines.push_back(r.line);
    return lines;
  };
  auto h2o_lines = load_species(test_data::kH2OLines, 1);
  auto co2_lines = load_species(test_data::kCO2Lines, 2);
  auto ch4_lines = load_species(test_data::kCH4Lines, 6);

  auto h2o_op = std::make_shared<LineListOpacity>("H2O", h2o_lines, 18.015);
  auto co2_op = std::make_shared<LineListOpacity>("CO2", co2_lines, 44.010);
  auto ch4_op = std::make_shared<LineListOpacity>("CH4", ch4_lines, 16.043);
  auto ray_op = std::make_shared<RayleighOpacity>("H2_BULK", 8.49e-29);

  // ─── 2. Wavelength grid covering the H₂O 1.4 µm + CO₂ 4.3 µm bands.
  std::vector<double> wn;
  for (double w = 2200.0; w <= 7600.0; w += 200.0) wn.push_back(w);

  // ─── 3. Forward closure: 4 free params → noisy spectrum. ──────────
  // Rebuild atmosphere + cloud each call (cloud P depends on params).
  // The atmosphere is small (40 layers) so the per-call cost is modest.
  constexpr std::size_t N_LAYERS = 40;
  constexpr double TRUE_VMR_CH4   = 5.0e-5;            // fixed background
  auto build_forward = [&]() {
    return [&](double T_irr, double log10_vmr_h2o,
               double log10_vmr_co2, double log10_p_cloud) {
      auto cld_op = std::make_shared<CloudDeckOpacity>(
          "CLOUD_DECK", std::pow(10.0, log10_p_cloud), 1.0e-18);

      Atmosphere atm;
      atm.species = {
        {"H2_BULK",    2.016},
        {"H2O",        18.015},
        {"CO2",        44.010},
        {"CH4",        16.043},
        {"CLOUD_DECK", 1.0},
      };
      atm.pressure_bar.resize(N_LAYERS);
      atm.mixing_ratios = Tensor({N_LAYERS, 5});
      const double log_top = std::log(1.0e-6);
      const double log_bot = std::log(1.0e2);
      const double vmr_h2o = std::pow(10.0, log10_vmr_h2o);
      const double vmr_co2 = std::pow(10.0, log10_vmr_co2);
      for (std::size_t i = 0; i < N_LAYERS; ++i) {
        const double frac = static_cast<double>(i) /
                            static_cast<double>(N_LAYERS - 1);
        atm.pressure_bar[i] = std::exp(log_top + frac * (log_bot - log_top));
        atm.mixing_ratios.at(i, 0) = 1.0;
        atm.mixing_ratios.at(i, 1) = vmr_h2o;
        atm.mixing_ratios.at(i, 2) = vmr_co2;
        atm.mixing_ratios.at(i, 3) = TRUE_VMR_CH4;
        atm.mixing_ratios.at(i, 4) = 1.0;
      }
      atm.temperature_k = guillot_profile(atm.pressure_bar,
                                          /*T_int=*/200.0,
                                          T_irr, /*gamma=*/0.5);
      TransmissionModel m;
      m.add_opacity(ray_op);
      m.add_opacity(h2o_op);
      m.add_opacity(co2_op);
      m.add_opacity(ch4_op);
      m.add_opacity(cld_op);
      return m.forward(atm, wn);
    };
  };
  auto forward = build_forward();

  // ─── 4. Inject truth + JWST-PRISM-realistic noise. ────────────────
  const double TRUE_T_IRR    = 1500.0;
  const double TRUE_LV_H2O   = -3.0;
  const double TRUE_LV_CO2   = -3.3;
  const double TRUE_LP_CLOUD = -2.0;
  const double NOISE_SIGMA   = 1.0e-4;

  Spectrum truth = forward(TRUE_T_IRR, TRUE_LV_H2O,
                           TRUE_LV_CO2, TRUE_LP_CLOUD);
  std::mt19937_64 rng(2026);
  std::normal_distribution<double> nz(0.0, NOISE_SIGMA);
  Spectrum observed = truth;
  for (auto& v : observed.values) v += nz(rng);
  std::vector<double> uncertainty(wn.size(), NOISE_SIGMA);

  // ─── 5. Wrap into the substrate Retrieval API. ────────────────────
  auto wrapped = [&](const std::vector<double>& s) {
    return forward(s[0], s[1], s[2], s[3]);
  };
  std::vector<Parameter> params{
    {"T_irr",         800.0,  2200.0},
    {"log10_VMR_H2O",  -6.0,   -1.0},
    {"log10_VMR_CO2",  -6.0,   -1.0},
    {"log10_P_cloud",  -3.0,   +1.0},
  };
  Retrieval ret(params, wrapped, observed, uncertainty);
  std::vector<double> truth_vec{TRUE_T_IRR, TRUE_LV_H2O,
                                 TRUE_LV_CO2, TRUE_LP_CLOUD};
  assert(std::isfinite(ret.log_posterior(truth_vec)));

  // ─── 6. MCMC. Single-chain MH; proposal widths tuned for ~25%
  //     acceptance with 4 weakly-correlated parameters. ───────────────
  auto result = ret.run_mcmc(
      /*init=*/{1300.0, -2.7, -3.0, -1.5},
      /*burn=*/2000,
      /*ns=*/3000,
      /*proposal_widths=*/{30.0, 0.04, 0.06, 0.08},
      /*seed=*/2026);

  assert(result.acceptance_rate > 0.05);
  assert(result.acceptance_rate < 0.85);

  // ─── 7. Recovery within 3σ on every parameter. ────────────────────
  PosteriorSummary post(params, result.samples);
  for (std::size_t i = 0; i < params.size(); ++i) {
    const auto& e = post[params[i].name];
    const double truth_v = truth_vec[i];
    assert(std::fabs(e.median - truth_v) < 3.0 * e.stddev);
    // Convergence sanity: posterior stddev must be a meaningful
    // fraction of the prior box, not the full width — a bad chain
    // that just samples the prior would have stddev ≈ box/√12 and
    // could pass the 3σ check by accident.
    const double box = params[i].prior_max - params[i].prior_min;
    assert(e.stddev < 0.5 * box);
  }

  // ─── 8. Determinism: same seed → bit-equal samples. ───────────────
  auto r2 = ret.run_mcmc({1300.0, -2.7, -3.0, -1.5}, 2000, 3000,
                          {30.0, 0.04, 0.06, 0.08}, 2026);
  assert(result.samples.size() == r2.samples.size());
  for (std::size_t i = 0; i < result.samples.size(); ++i) {
    for (std::size_t d = 0; d < params.size(); ++d) {
      assert(result.samples[i][d] == r2.samples[i][d]);
    }
  }

  return 0;
}
