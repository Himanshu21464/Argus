// M4 capstone: full SIE substrate retrieval. The same `argus::Retrieval`
// pipeline used for atmospheric and SIS retrievals now recovers all 7
// parameters of an elliptical-galaxy lens (θ_E, q, φ, lens_x, lens_y,
// source_x, source_y) from 4 observed quad-image positions.
//
// Approach: source-plane chi². For each observed image θ_i, the
// back-projection β_back_i = θ_i - α(θ_i; lens_params) gives an
// estimate of the source position. If the lens model is correct, all
// 4 back-projections coincide; the chi² measures their spread relative
// to the free `source` parameter. Avoids `find_images` inside the
// MCMC inner loop (the deflection alone is ~1 µs per call).
//
// Verification: posterior median within 3σ of every parameter; all 7
// params (including q and φ which are SIE-specific) constrained.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <random>
#include <vector>

#include "argus/argus.hpp"

int main() {
  using namespace argus;
  using namespace argus::lensing;

  // ─── 1. Inject SIE truth, solve for the 4 quad images. ─────────────
  const double TRUE_THETA_E = 1.0;
  const double TRUE_Q       = 0.7;
  const double TRUE_PHI     = 0.3;
  const double TRUE_LX      = 0.05;
  const double TRUE_LY      = -0.02;
  const double TRUE_SX      = 0.04;
  const double TRUE_SY      = 0.03;
  const double IMG_NOISE    = 0.01;     // 1-σ image-position noise (arcsec)

  SIE truth_lens(TRUE_THETA_E, TRUE_Q, TRUE_PHI, {TRUE_LX, TRUE_LY});
  Vec2 truth_source{TRUE_SX, TRUE_SY};
  auto imgs = find_images(truth_lens, truth_source, /*radius=*/2.0,
                          /*grid=*/120);
  assert(imgs.size() == 4);             // source inside tangential caustic

  // ─── 2. Add image-plane Gaussian noise. ────────────────────────────
  std::mt19937_64 rng(2026);
  std::normal_distribution<double> nz(0.0, IMG_NOISE);
  std::vector<Vec2> observed;
  observed.reserve(4);
  for (const auto& im : imgs) {
    observed.push_back({im.theta.x + nz(rng), im.theta.y + nz(rng)});
  }

  // ─── 3. Source-plane forward closure. ──────────────────────────────
  // 7 free params: θ_E, q, φ, lens_x, lens_y, source_x, source_y.
  // For each observed image θ_i, back-project β_back_i = θ_i - α(θ_i)
  // and pack (β_back_i.x − source_x, β_back_i.y − source_y) into
  // the spectrum. With observation = 8 zeros, χ² becomes Σ |β_back_i
  // − source|²/σ_β².
  Retrieval::Forward fwd =
      [&observed](const std::vector<double>& s) -> Spectrum {
    Spectrum out;
    out.wavenumber_cm.reserve(8);
    out.values.reserve(8);
    if (!(s[0] > 0.0 && s[1] > 0.0 && s[1] <= 1.0)) {
      out.wavenumber_cm.assign(8, 0.0);
      out.values.assign(8, 1.0e6);
      return out;
    }
    SIE lens(s[0], s[1], s[2], {s[3], s[4]});
    const double sx = s[5], sy = s[6];
    for (std::size_t i = 0; i < observed.size(); ++i) {
      Vec2 a = lens.deflection(observed[i]);
      const double bx = observed[i].x - a.x;
      const double by = observed[i].y - a.y;
      out.wavenumber_cm.push_back(static_cast<double>(2 * i));
      out.wavenumber_cm.push_back(static_cast<double>(2 * i + 1));
      out.values.push_back(bx - sx);
      out.values.push_back(by - sy);
    }
    return out;
  };

  // Source-plane noise σ_β ~ σ_θ / √⟨μ⟩; use 0.02 arcsec uniformly.
  Spectrum obs_zero;
  obs_zero.wavenumber_cm.assign(8, 0.0);
  obs_zero.values.assign(8, 0.0);
  std::vector<double> sigma_beta(8, 0.02);

  std::vector<Parameter> params{
    {"theta_E",  0.5,  2.0},
    {"q",        0.3,  1.0},
    {"phi",     -1.5,  1.5},
    {"lens_x",  -0.3,  0.3},
    {"lens_y",  -0.3,  0.3},
    {"source_x",-0.3,  0.3},
    {"source_y",-0.3,  0.3},
  };
  Retrieval ret(params, fwd, obs_zero, sigma_beta);
  std::vector<double> truth{TRUE_THETA_E, TRUE_Q, TRUE_PHI,
                             TRUE_LX, TRUE_LY, TRUE_SX, TRUE_SY};
  assert(std::isfinite(ret.log_posterior(truth)));

  // ─── 4. EnsembleSampler — 32 walkers around a wrong initial guess. ─
  auto logp = [&](const std::vector<double>& v) {
    return ret.log_posterior(v);
  };
  const std::vector<double> init{0.95, 0.8, 0.0, 0.0, 0.0, 0.05, 0.05};
  std::mt19937_64 wr(2026);
  std::normal_distribution<double> wj(0.0, 0.01);
  std::vector<std::vector<double>> walkers;
  walkers.reserve(32);
  for (std::size_t w = 0; w < 32; ++w) {
    std::vector<double> walker(7);
    for (int attempt = 0; attempt < 50; ++attempt) {
      for (std::size_t i = 0; i < 7; ++i) walker[i] = init[i] + wj(wr);
      if (std::isfinite(logp(walker))) break;
    }
    assert(std::isfinite(logp(walker)));
    walkers.push_back(walker);
  }

  EnsembleSampler sampler(logp, walkers, 2.0, /*seed=*/2026);
  sampler.burn_in(2000);
  auto result = sampler.sample(2000);

  assert(sampler.acceptance_rate() > 0.05);
  assert(sampler.acceptance_rate() < 0.95);

  // ─── 5. Recovery within 3σ on every SIE parameter. ─────────────────
  PosteriorSummary post(params, result.samples);
  for (std::size_t i = 0; i < params.size(); ++i) {
    const auto& e = post[params[i].name];
    const double truth_v = truth[i];
    assert(std::fabs(e.median - truth_v) < 3.0 * e.stddev);
  }

  // ─── 6. Determinism: same seed → bit-equal samples. ────────────────
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
