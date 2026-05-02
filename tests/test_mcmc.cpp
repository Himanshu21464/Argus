// Hard tests for the Metropolis-Hastings sampler.
// Recovers known-truth parameters from analytic log-posteriors.

#include <cassert>
#include <cmath>
#include <stdexcept>

#include "argus/argus.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

// Mean of a sample column.
double mean(const std::vector<std::vector<double>>& s, std::size_t d) {
  double sum = 0.0;
  for (const auto& v : s) sum += v[d];
  return sum / static_cast<double>(s.size());
}

// Standard deviation of a sample column.
double stddev(const std::vector<std::vector<double>>& s, std::size_t d) {
  const double m = mean(s, d);
  double sq = 0.0;
  for (const auto& v : s) { const double diff = v[d] - m; sq += diff * diff; }
  return std::sqrt(sq / static_cast<double>(s.size()));
}

}  // namespace

int main() {
  using namespace argus;

  // ─── 1. 1D Gaussian: log-p(x) = -0.5 * (x - mu)^2 / sigma^2 ──────────
  // Recover mean and standard deviation to <5%.
  {
    const double TRUE_MU    = 3.0;
    const double TRUE_SIGMA = 1.5;
    auto logp = [&](const std::vector<double>& x) {
      const double r = (x[0] - TRUE_MU) / TRUE_SIGMA;
      return -0.5 * r * r;
    };
    MetropolisHastings sampler(logp, {1.5}, /*seed=*/12345);
    std::vector<double> state{0.0};
    sampler.burn_in(state, /*n_steps=*/2000);
    auto out = sampler.sample(state, /*n_samples=*/20000);

    assert(close(mean(out.samples, 0),    TRUE_MU,    0.05));
    assert(close(stddev(out.samples, 0),  TRUE_SIGMA, 0.05));
    // Acceptance rate should be reasonable.
    const double a = sampler.acceptance_rate();
    assert(a > 0.20 && a < 0.80);
  }

  // ─── 2. 2D uncorrelated Gaussian. ────────────────────────────────────
  {
    const double mu1 = -2.0, mu2 = 5.0;
    const double s1 = 0.5,   s2 = 2.0;
    auto logp = [&](const std::vector<double>& x) {
      const double r1 = (x[0] - mu1) / s1;
      const double r2 = (x[1] - mu2) / s2;
      return -0.5 * (r1 * r1 + r2 * r2);
    };
    MetropolisHastings sampler(logp, {0.5, 2.0}, /*seed=*/42);
    std::vector<double> state{0.0, 0.0};
    sampler.burn_in(state, 3000);
    auto out = sampler.sample(state, 30000);

    assert(close(mean(out.samples, 0),   mu1, 0.05));
    assert(close(mean(out.samples, 1),   mu2, 0.05));
    assert(close(stddev(out.samples, 0), s1,  0.10));
    assert(close(stddev(out.samples, 1), s2,  0.10));
  }

  // ─── 3. Determinism: same seed -> bit-equal samples. ─────────────────
  {
    auto logp = [](const std::vector<double>& x) { return -0.5 * x[0] * x[0]; };
    MetropolisHastings a(logp, {1.0}, 7);
    MetropolisHastings b(logp, {1.0}, 7);
    std::vector<double> sa{0.5}, sb{0.5};
    auto ra = a.sample(sa, 200);
    auto rb = b.sample(sb, 200);
    for (std::size_t i = 0; i < ra.samples.size(); ++i) {
      assert(ra.samples[i][0] == rb.samples[i][0]);
    }
  }

  // ─── 4. Bad inputs throw. ────────────────────────────────────────────
  {
    auto logp = [](const std::vector<double>& x) { return -x[0] * x[0]; };
    bool threw = false;
    try { MetropolisHastings(nullptr, {1.0}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { MetropolisHastings(logp, {}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { MetropolisHastings(logp, {0.0}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { MetropolisHastings(logp, {-1.0}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // ─── 5. Initial state with -inf log-posterior throws on burn_in. ────
  {
    auto logp = [](const std::vector<double>& x) {
      return (x[0] > 0.0)
        ? -0.5 * x[0] * x[0]
        : -std::numeric_limits<double>::infinity();
    };
    MetropolisHastings sampler(logp, {0.5}, 1);
    std::vector<double> state{-1.0};
    bool threw = false;
    try { sampler.burn_in(state, 100); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  return 0;
}
