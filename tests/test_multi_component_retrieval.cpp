// M5 capstone: 2-component (core + jet) interferometric retrieval over
// an Earth-rotation UV synthesis track.
//
// Real radio sources (AGN, quasars, starburst galaxies) are typically
// multi-component on the relevant angular scales — a compact core
// plus extended jet / dust emission. Argus's substrate claim now
// covers this end-to-end: the same `Retrieval` API recovers all 8
// parameters (l, m, F, σ for each of two Gaussian components) from
// noisy synthetic visibilities sampled across a 5-hour HA track on a
// 7-antenna VLA-like array (105 baselines × 2 = 210 observations).
//
// Label-switching degeneracy is broken via non-overlapping size priors
// (σ_core ∈ [0, 2e-5], σ_jet ∈ [2e-5, 8e-5]). Walkers initialised
// near a deliberately wrong guess; recovery within 3σ on every
// parameter via EnsembleSampler.

#include <cassert>
#include <cmath>
#include <random>
#include <vector>

#include "argus/argus.hpp"

int main() {
  using namespace argus;
  using namespace argus::interferometry;

  // ─── 1. UV coverage from a 5-hour track. ────────────────────────────
  std::vector<double> east { 0.0, -250.0,  250.0,    0.0, -500.0,  500.0,    0.0};
  std::vector<double> north{ 0.0, -433.0, -433.0,  500.0, -866.0, -866.0, 1000.0};
  std::vector<double> ha;
  for (int k = -2; k <= 2; ++k) ha.push_back(k * 0.5);     // 5 HAs
  const double dec = 30.0 * M_PI / 180.0;
  const double lat = 34.0 * M_PI / 180.0;
  const double wavelength = 0.21;                          // HI 21 cm
  const std::vector<UVPoint> uv =
      uv_coverage_track(east, north, lat, ha, dec, wavelength);
  assert(uv.size() == 21 * 5);                             // 21 baselines × 5 HAs

  // ─── 2. 2-component forward closure. ────────────────────────────────
  // 8 free params: (l_core, m_core, F_core, σ_core,
  //                 l_jet,  m_jet,  F_jet,  σ_jet).
  Retrieval::Forward fwd =
      [&uv](const std::vector<double>& s) -> Spectrum {
    Spectrum out;
    out.wavenumber_cm.reserve(2 * uv.size());
    out.values.reserve(2 * uv.size());
    if (s[2] <= 0.0 || s[3] < 0.0 || s[6] <= 0.0 || s[7] < 0.0) {
      out.wavenumber_cm.assign(2 * uv.size(), 0.0);
      out.values.assign(2 * uv.size(), 1.0e6);
      return out;
    }
    GaussianSource core{s[0], s[1], s[2], s[3]};
    GaussianSource jet {s[4], s[5], s[6], s[7]};
    auto vs = predict_visibilities(
        std::vector<GaussianSource>{core, jet}, uv);
    for (std::size_t i = 0; i < vs.size(); ++i) {
      out.wavenumber_cm.push_back(static_cast<double>(2 * i));
      out.wavenumber_cm.push_back(static_cast<double>(2 * i + 1));
      out.values.push_back(vs[i].real);
      out.values.push_back(vs[i].imag);
    }
    return out;
  };

  // ─── 3. Inject truth + Gaussian noise. ─────────────────────────────
  std::vector<double> truth{0.0, 0.0, 1.0, 1.0e-5,
                             1.5e-5, 5.0e-6, 0.3, 3.0e-5};
  Spectrum clean = fwd(truth);

  std::mt19937_64 rng(2026);
  const double NOISE_SIGMA = 0.02;
  std::normal_distribution<double> nz(0.0, NOISE_SIGMA);
  Spectrum observed = clean;
  for (auto& v : observed.values) v += nz(rng);
  std::vector<double> uncertainty(observed.values.size(), NOISE_SIGMA);

  // ─── 4. Retrieval setup with non-overlapping size priors to break
  //     the label-switching degeneracy. ────────────────────────────
  std::vector<Parameter> params{
    {"l_core",    -3.0e-5, 3.0e-5},
    {"m_core",    -3.0e-5, 3.0e-5},
    {"F_core",     0.1,    3.0},
    {"sigma_core", 0.0,    2.0e-5},
    {"l_jet",     -3.0e-5, 3.0e-5},
    {"m_jet",     -3.0e-5, 3.0e-5},
    {"F_jet",      0.05,   1.5},
    {"sigma_jet",  2.0e-5, 8.0e-5},
  };
  Retrieval ret(params, fwd, observed, uncertainty);
  assert(std::isfinite(ret.log_posterior(truth)));

  // ─── 5. Initialise 32 walkers near a deliberately wrong guess. ─────
  auto logp = [&](const std::vector<double>& v) {
    return ret.log_posterior(v);
  };
  const std::vector<double> init{1.0e-6, -1.0e-6, 0.8, 1.5e-5,
                                  1.0e-5,  1.0e-5, 0.4, 2.5e-5};
  std::mt19937_64 wr(2026);
  std::normal_distribution<double> wj_pos (0.0, 1.0e-6);
  std::normal_distribution<double> wj_flux(0.0, 0.02);
  std::normal_distribution<double> wj_sig (0.0, 5.0e-7);
  std::vector<std::vector<double>> walkers;
  walkers.reserve(32);
  for (std::size_t w = 0; w < 32; ++w) {
    std::vector<double> walker(8);
    for (int attempt = 0; attempt < 100; ++attempt) {
      walker[0] = init[0] + wj_pos (wr);
      walker[1] = init[1] + wj_pos (wr);
      walker[2] = init[2] + wj_flux(wr);
      walker[3] = init[3] + wj_sig (wr);
      walker[4] = init[4] + wj_pos (wr);
      walker[5] = init[5] + wj_pos (wr);
      walker[6] = init[6] + wj_flux(wr);
      walker[7] = init[7] + wj_sig (wr);
      if (std::isfinite(logp(walker))) break;
    }
    assert(std::isfinite(logp(walker)));
    walkers.push_back(walker);
  }

  EnsembleSampler sampler(logp, walkers, /*stretch_a=*/2.0, /*seed=*/2026);
  sampler.burn_in(2000);
  auto result = sampler.sample(2000);

  assert(sampler.acceptance_rate() > 0.05);
  assert(sampler.acceptance_rate() < 0.95);

  // ─── 6. Recovery within 3σ on every parameter (8 of them). ─────────
  PosteriorSummary post(params, result.samples);
  for (std::size_t i = 0; i < params.size(); ++i) {
    const auto& e = post[params[i].name];
    const double truth_v = truth[i];
    assert(std::fabs(e.median - truth_v) < 3.0 * e.stddev);
    // Convergence sanity: chain mustn't merely sample the prior.
    const double box = params[i].prior_max - params[i].prior_min;
    assert(e.stddev < 0.5 * box);
  }

  // ─── 7. Determinism: same seed → bit-equal samples. ────────────────
  EnsembleSampler s2(logp, walkers, 2.0, 2026);
  s2.burn_in(2000);
  auto r2 = s2.sample(2000);
  assert(result.samples.size() == r2.samples.size());
  for (std::size_t i = 0; i < result.samples.size(); ++i) {
    for (std::size_t d = 0; d < params.size(); ++d) {
      assert(result.samples[i][d] == r2.samples[i][d]);
    }
  }

  return 0;
}
