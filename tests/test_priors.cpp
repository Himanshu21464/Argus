// Tests for the extended prior types: Uniform / Gaussian / LogUniform.
//
// We use a degenerate retrieval (constant model = observation) so the
// likelihood is constant and the posterior IS the prior. Then we
// recover the prior via MCMC and check the moments match analytic.

#include <cassert>
#include <cmath>
#include <random>
#include <vector>

#include "argus/argus.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

double mean_col(const std::vector<std::vector<double>>& s, std::size_t d) {
  double sum = 0.0;
  for (const auto& v : s) sum += v[d];
  return sum / static_cast<double>(s.size());
}

double std_col(const std::vector<std::vector<double>>& s, std::size_t d) {
  const double m = mean_col(s, d);
  double sq = 0.0;
  for (const auto& v : s) { const double r = v[d] - m; sq += r * r; }
  return std::sqrt(sq / static_cast<double>(s.size()));
}

}  // namespace

int main() {
  using namespace argus;

  // Degenerate "observation" — model returns the same spectrum as the
  // observation regardless of state, so likelihood is constant and the
  // posterior reduces to the prior.
  Spectrum obs;
  obs.wavenumber_cm = {1.0, 2.0};
  obs.values        = {0.5, 0.5};
  std::vector<double> sigma{1.0, 1.0};
  auto constant_forward = [obs](const std::vector<double>&) { return obs; };

  // ─── 1. Uniform prior on [-2, 4] ────────────────────────────────────
  // Posterior = Uniform(-2, 4): mean = 1.0, std = sqrt((4 - -2)^2 / 12) = √3
  {
    std::vector<Parameter> params{
      {"x", -2.0, 4.0, PriorType::Uniform, 0.0, 1.0},
    };
    Retrieval ret(params, constant_forward, obs, sigma);
    auto r = ret.run_mcmc({0.0}, /*burn=*/2000, /*ns=*/30000,
                          /*proposal_widths=*/{1.5}, /*seed=*/7);
    assert(close(mean_col(r.samples, 0), 1.0,            0.05, 0.10));
    assert(close(std_col(r.samples, 0),  std::sqrt(3.0), 0.10));
  }

  // ─── 2. Gaussian prior, mean=0, stddev=1, hard-clipped to [-3, 3].
  // The clip is at ±3σ where the truncated Gaussian still ≈ standard,
  // so mean ≈ 0, std ≈ 1.
  {
    std::vector<Parameter> params{
      {"x", -3.0, 3.0, PriorType::Gaussian, /*mean=*/0.0, /*stddev=*/1.0},
    };
    Retrieval ret(params, constant_forward, obs, sigma);
    auto r = ret.run_mcmc({0.0}, /*burn=*/2000, /*ns=*/30000,
                          /*proposal_widths=*/{0.7}, /*seed=*/13);
    // MCMC noise on N(0,1) mean from 30k samples: ±~σ/√Neff ≈ 0.05-0.10.
    // Loosen the bound to absorb seed-dependent variation.
    assert(std::fabs(mean_col(r.samples, 0) - 0.0) < 0.10);
    assert(close(std_col(r.samples, 0), 1.0, 0.10));
  }

  // ─── 3. Gaussian prior centred away from zero.
  // Prior N(5, 0.5), clipped to [3, 7] (well outside ±3σ).
  // Recovered mean ≈ 5, std ≈ 0.5.
  {
    std::vector<Parameter> params{
      {"x", 3.0, 7.0, PriorType::Gaussian, /*mean=*/5.0, /*stddev=*/0.5},
    };
    Retrieval ret(params, constant_forward, obs, sigma);
    auto r = ret.run_mcmc({5.0}, /*burn=*/2000, /*ns=*/30000,
                          /*proposal_widths=*/{0.4}, /*seed=*/27);
    assert(close(mean_col(r.samples, 0), 5.0, 0.02));
    assert(close(std_col(r.samples, 0),  0.5, 0.10));
  }

  // ─── 4. LogUniform prior on [1, 100].
  // p(x) ∝ 1/x. Mean = 99/ln(100) ≈ 21.5
  //               std = sqrt(E[x^2] - mean^2) where E[x^2] = (100^2 - 1)/(2 ln 100) ≈ 1085
  //               std ≈ sqrt(1085 - 21.5^2) = sqrt(1085 - 462) ≈ √623 ≈ 25.0
  {
    std::vector<Parameter> params{
      {"x", 1.0, 100.0, PriorType::LogUniform, 0.0, 1.0},
    };
    Retrieval ret(params, constant_forward, obs, sigma);
    auto r = ret.run_mcmc({10.0}, /*burn=*/3000, /*ns=*/30000,
                          /*proposal_widths=*/{8.0}, /*seed=*/41);
    const double m = mean_col(r.samples, 0);
    const double s = std_col(r.samples, 0);
    // Loose tolerance because LogUniform with wide range mixes slowly.
    assert(std::fabs(m - 21.5) / 21.5 < 0.20);
    assert(std::fabs(s - 25.0) / 25.0 < 0.25);
  }

  // ─── 5. Bounds clip applies to all prior types.
  // Set initial state outside the box -> log_posterior returns -inf.
  {
    std::vector<Parameter> params{
      {"x", -1.0, 1.0, PriorType::Gaussian, 0.0, 100.0},   // huge stddev
    };
    Retrieval ret(params, constant_forward, obs, sigma);
    assert(std::isinf(ret.log_posterior({1.5})));
    assert(std::isinf(ret.log_posterior({-1.5})));
    assert(std::isfinite(ret.log_posterior({0.5})));
  }

  // ─── 6. LogUniform with x ≤ 0 returns -inf even inside the bounds.
  // (Requires construction with prior_min > 0; we test with bound = 0.)
  {
    std::vector<Parameter> params{
      {"x", 0.0, 10.0, PriorType::LogUniform, 0.0, 1.0},
    };
    Retrieval ret(params, constant_forward, obs, sigma);
    assert(std::isinf(ret.log_posterior({0.0})));
    assert(std::isfinite(ret.log_posterior({1.0})));
  }

  return 0;
}
