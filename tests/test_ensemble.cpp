// Tests for the affine-invariant ensemble sampler (Goodman-Weare 2010).
// The killer test: sample from a *highly correlated* Gaussian where
// single-chain Metropolis-Hastings struggles, and verify the ensemble
// recovers the correct mean and (rotated) covariance.

#include <cassert>
#include <cmath>
#include <random>
#include <stdexcept>
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

  // ─── 1. 2D uncorrelated Gaussian ────────────────────────────────────
  // Recover means and stds with a small ensemble.
  {
    const double mu0 = -1.0, mu1 = 4.0;
    const double s0 = 0.5,   s1 = 1.5;
    auto logp = [&](const std::vector<double>& x) {
      const double r0 = (x[0] - mu0) / s0;
      const double r1 = (x[1] - mu1) / s1;
      return -0.5 * (r0 * r0 + r1 * r1);
    };
    // 16 walkers initialised in a small ball around the truth.
    std::mt19937_64 rng(123);
    std::normal_distribution<double> jitter(0.0, 0.1);
    std::vector<std::vector<double>> walkers(16, std::vector<double>(2));
    for (auto& w : walkers) {
      w[0] = mu0 + jitter(rng);
      w[1] = mu1 + jitter(rng);
    }

    EnsembleSampler s(logp, walkers, /*stretch_a=*/2.0, /*seed=*/7);
    s.burn_in(500);
    auto out = s.sample(2000);

    assert(close(mean_col(out.samples, 0), mu0, 0.05, 0.05));
    assert(close(mean_col(out.samples, 1), mu1, 0.05));
    assert(close(std_col(out.samples,  0),  s0, 0.10));
    assert(close(std_col(out.samples,  1),  s1, 0.10));

    const double a = s.acceptance_rate();
    assert(a > 0.10 && a < 0.90);
  }

  // ─── 2. *Highly correlated* 2D Gaussian — the affine-invariance test
  //     The covariance has condition number 100; single-chain MH would
  //     mix slowly along the principal axis. The ensemble sampler
  //     should still recover the marginal stds.
  //
  //     Σ in basis (x, y):
  //         var(x) = 1.0
  //         var(y) = 1.0
  //         cov(x,y) = 0.99
  //     Inverse:
  //         A = (1/det) * [[1, -0.99], [-0.99, 1]]
  //         det = 1 - 0.99^2 = 0.0199
  //
  //     log_p = -0.5 * x^T A x
  {
    const double det = 1.0 - 0.99 * 0.99;
    auto logp = [&](const std::vector<double>& v) {
      const double x = v[0], y = v[1];
      // Quadratic form: (x^2 - 2*0.99*x*y + y^2) / det
      const double q = (x * x - 2.0 * 0.99 * x * y + y * y) / det;
      return -0.5 * q;
    };
    std::mt19937_64 rng(2026);
    std::normal_distribution<double> jitter(0.0, 0.5);
    std::vector<std::vector<double>> walkers(20, std::vector<double>(2));
    for (auto& w : walkers) {
      w[0] = jitter(rng);
      w[1] = jitter(rng);
    }
    EnsembleSampler s(logp, walkers, 2.0, /*seed=*/9999);
    s.burn_in(2000);
    auto out = s.sample(5000);

    // Marginal means should be ≈ 0
    assert(std::fabs(mean_col(out.samples, 0)) < 0.10);
    assert(std::fabs(mean_col(out.samples, 1)) < 0.10);
    // Marginal stds should be ≈ 1 (diagonal of Σ)
    assert(close(std_col(out.samples, 0), 1.0, 0.20));
    assert(close(std_col(out.samples, 1), 1.0, 0.20));
  }

  // ─── 3. Determinism: same seed -> bit-equal samples. ────────────────
  {
    auto logp = [](const std::vector<double>& v) {
      return -0.5 * (v[0] * v[0] + v[1] * v[1]);
    };
    std::vector<std::vector<double>> w(8, std::vector<double>{0.3, 0.4});
    EnsembleSampler a(logp, w, 2.0, 42);
    EnsembleSampler b(logp, w, 2.0, 42);
    auto ra = a.sample(50);
    auto rb = b.sample(50);
    assert(ra.samples.size() == rb.samples.size());
    for (std::size_t i = 0; i < ra.samples.size(); ++i) {
      assert(ra.samples[i][0] == rb.samples[i][0]);
      assert(ra.samples[i][1] == rb.samples[i][1]);
    }
  }

  // ─── 4. Bad inputs throw. ───────────────────────────────────────────
  {
    auto logp = [](const std::vector<double>& v) { return -v[0] * v[0]; };
    bool threw = false;
    // null callable
    try { EnsembleSampler(nullptr, {{0.0}, {1.0}, {2.0}, {3.0}}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    // too few walkers
    threw = false;
    try { EnsembleSampler(logp, {{0.0}, {1.0}}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    // odd number of walkers
    threw = false;
    try {
      EnsembleSampler(logp, {{0.0}, {1.0}, {2.0}, {3.0}, {4.0}});
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    // dimensionality mismatch
    threw = false;
    try {
      EnsembleSampler(logp, {{0.0}, {1.0}, {2.0, 3.0}, {4.0}});
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    // stretch_a <= 1 invalid
    threw = false;
    try {
      EnsembleSampler(logp, {{0.0}, {1.0}, {2.0}, {3.0}}, /*a=*/0.5);
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    // walker with -inf logp at init
    auto neg_inf_logp = [](const std::vector<double>& v) {
      return (v[0] >= 0.0)
        ? -v[0] * v[0]
        : -std::numeric_limits<double>::infinity();
    };
    threw = false;
    try {
      EnsembleSampler(neg_inf_logp,
                      {{1.0}, {2.0}, {-1.0}, {3.0}});
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // ─── 5. Result shape matches API contract ───────────────────────────
  {
    auto logp = [](const std::vector<double>& v) {
      return -0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    };
    std::vector<std::vector<double>> w(10, std::vector<double>{0.0, 0.0, 0.0});
    EnsembleSampler s(logp, w, 2.0, 1);
    auto out = s.sample(7);
    assert(out.n_steps   == 7);
    assert(out.n_walkers == 10);
    assert(out.n_dim     == 3);
    assert(out.samples.size() == 7 * 10);
    assert(out.log_posteriors.size() == 7 * 10);
    for (const auto& s_row : out.samples) {
      assert(s_row.size() == 3);
    }
  }

  return 0;
}
