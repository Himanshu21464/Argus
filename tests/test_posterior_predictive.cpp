// Posterior-predictive check tests.
//
// Generate a synthetic spectrum at known truth, run a retrieval, then
// compute the posterior-predictive 16-84% band. Assert that the
// observed spectrum lies inside the band at the vast majority of
// wavelengths (the standard "are we converged?" check).

#include <cassert>
#include <cmath>
#include <random>
#include <sstream>

#include "argus/argus.hpp"
#include "argus/test_data.hpp"

int main() {
  using namespace argus;

  // ─── Setup: real H2O retrieval over the 1.4 μm band. ────────────────
  std::istringstream is{std::string(test_data::kH2OLines)};
  auto records = Hitran::load(is, /*filter=*/1);
  std::vector<Line> lines;
  for (const auto& r : records) lines.push_back(r.line);
  auto opacity = std::make_shared<LineListOpacity>("H2O", lines, 18.015);

  std::vector<double> wn;
  for (int i = 0; i < 22; ++i) wn.push_back(6900.0 + i * 17.0);
  Species h2o{"H2O", 18.015};

  auto forward = [&](double T_k, double log10_vmr) {
    Atmosphere atm = isothermal(T_k, 1.0e-6, 1.0e2, 40, h2o,
                                std::pow(10.0, log10_vmr));
    TransmissionModel m;
    m.add_opacity(opacity);
    return m.forward(atm, wn);
  };

  const double TRUE_T  = 1500.0;
  const double TRUE_LV = -3.0;
  Spectrum truth = forward(TRUE_T, TRUE_LV);

  std::mt19937_64 rng(2026);
  const double noise_sigma = 5.0e-5;
  std::normal_distribution<double> noise(0.0, noise_sigma);
  Spectrum observed = truth;
  for (auto& v : observed.values) v += noise(rng);
  std::vector<double> uncertainty(wn.size(), noise_sigma);

  std::vector<Parameter> params{
    {"T_K",       800.0,  2200.0},
    {"log10_VMR", -6.0,   -1.0},
  };
  Retrieval ret(params,
                [&](const std::vector<double>& s) {
                  return forward(s[0], s[1]);
                },
                observed, uncertainty);

  auto result = ret.run_mcmc(/*init=*/{1300.0, -2.7},
                             /*burn=*/3000,
                             /*ns=*/5000,
                             /*proposal_widths=*/{20.0, 0.04},
                             /*seed=*/99);

  // ─── 1. Default quantiles {0.16, 0.50, 0.84}. ───────────────────────
  {
    auto pp = ret.posterior_predictive(result.samples, /*thin=*/10);
    assert(pp.wavenumber_cm.size() == wn.size());
    assert(pp.quantiles.size() == 3);
    assert(pp.bands.size() == wn.size());
    for (const auto& band : pp.bands) {
      assert(band.size() == 3);
      // Order: q16 ≤ q50 ≤ q84
      assert(band[0] <= band[1]);
      assert(band[1] <= band[2]);
    }

    // ─── 2. Posterior-predictive coverage: observation should fall
    //     within 2σ_noise of the median band at >85% of wavelengths
    //     for a converged retrieval (the model band is much narrower
    //     than the noise, so the observed-vs-model residuals are
    //     dominated by Gaussian noise: 95% within ±2σ).
    int inside = 0;
    for (std::size_t i = 0; i < wn.size(); ++i) {
      if (observed.values[i] >= pp.bands[i][0] - 2.0 * noise_sigma &&
          observed.values[i] <= pp.bands[i][2] + 2.0 * noise_sigma) {
        ++inside;
      }
    }
    const double frac = static_cast<double>(inside) /
                        static_cast<double>(wn.size());
    assert(frac > 0.85);

    // ─── 3. Median should be close to truth at every wavelength. ──────
    for (std::size_t i = 0; i < wn.size(); ++i) {
      const double rel_err = std::fabs(pp.bands[i][1] - truth.values[i]) /
                             truth.values[i];
      // The median posterior-predictive should be within 1% of truth
      // for a converged retrieval.
      assert(rel_err < 0.01);
    }
  }

  // ─── 4. Custom quantiles. ───────────────────────────────────────────
  {
    auto pp = ret.posterior_predictive(result.samples, /*thin=*/20,
                                       /*quantiles=*/{0.025, 0.5, 0.975});
    assert(pp.quantiles.size() == 3);
    for (const auto& band : pp.bands) {
      assert(band[0] <= band[1]);
      assert(band[1] <= band[2]);
    }
  }

  // ─── 5. Thin = 1 (no thinning) consistent with thin = 5 to <1%. ─────
  {
    auto pp1 = ret.posterior_predictive(result.samples, 1);
    auto pp5 = ret.posterior_predictive(result.samples, 5);
    for (std::size_t i = 0; i < wn.size(); ++i) {
      const double rel = std::fabs(pp1.bands[i][1] - pp5.bands[i][1]) /
                         pp1.bands[i][1];
      assert(rel < 0.005);
    }
  }

  // ─── 6. Bad inputs throw. ───────────────────────────────────────────
  {
    bool threw = false;
    try { (void)ret.posterior_predictive({}, 1); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    threw = false;
    try {
      (void)ret.posterior_predictive(result.samples, 1, {0.5, 1.5});
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  return 0;
}
