// Tests for the Hamiltonian Monte Carlo sampler.
//
// HMC uses the forward-mode autograd path through Dual<T> to compute
// gradients, so the user's log-posterior must be templated. We test
// against analytic Gaussians where the gradient is closed-form.

#include <cassert>
#include <cmath>
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

  // ─── 1. argus::grad on a 1-D quadratic. ─────────────────────────────
  {
    auto f = [](const std::vector<Dual<double>>& x) {
      return x[0] * x[0] + Dual<double>(3.0) * x[0] - Dual<double>(7.0);
    };
    // ∇f(x) = 2x + 3
    for (double x : {-2.0, 0.0, 1.5, 3.7}) {
      auto g = grad(f, std::vector<double>{x});
      assert(close(g[0], 2.0 * x + 3.0, 1.0e-12));
    }
  }

  // ─── 2. argus::grad on a 2-D quadratic. ─────────────────────────────
  {
    auto f = [](const std::vector<Dual<double>>& x) {
      return x[0] * x[0] + Dual<double>(2.0) * x[1] * x[1]
             + Dual<double>(3.0) * x[0] * x[1];
    };
    // ∂f/∂x0 = 2x0 + 3x1
    // ∂f/∂x1 = 4x1 + 3x0
    auto g = grad(f, std::vector<double>{1.0, -2.0});
    assert(close(g[0], 2.0 * 1.0 + 3.0 * (-2.0), 1.0e-12));
    assert(close(g[1], 4.0 * (-2.0) + 3.0 * 1.0, 1.0e-12));
  }

  // ─── 3. HMC on a 2-D Gaussian: recover mean and stddev. ─────────────
  {
    const double mu0 = -1.0, mu1 = 4.0;
    const double s0 = 0.5,   s1 = 1.5;
    auto logp = [=](const std::vector<Dual<double>>& x) {
      Dual<double> r0 = (x[0] - Dual<double>(mu0)) / Dual<double>(s0);
      Dual<double> r1 = (x[1] - Dual<double>(mu1)) / Dual<double>(s1);
      return Dual<double>(-0.5) * (r0 * r0 + r1 * r1);
    };

    HMC sampler(logp, /*step_size=*/0.2, /*n_leapfrog=*/15, /*seed=*/2026);
    std::vector<double> state{0.0, 0.0};
    auto out = sampler.sample(state, 5000);
    // Drop first 1000 as burn-in.
    std::vector<std::vector<double>> post(out.samples.begin() + 1000,
                                          out.samples.end());
    assert(std::fabs(mean_col(post, 0) - mu0) < 0.10);
    assert(std::fabs(mean_col(post, 1) - mu1) < 0.20);
    assert(close(std_col(post, 0), s0, 0.20));
    assert(close(std_col(post, 1), s1, 0.20));
    // HMC with small step gives near-100% acceptance; larger step
    // gives ~70-80%. Either is acceptable as long as it isn't 0.
    const double a = sampler.acceptance_rate();
    assert(a > 0.30);
  }

  // ─── 4. Determinism: same seed -> identical samples. ─────────────────
  {
    auto logp = [](const std::vector<Dual<double>>& x) {
      return Dual<double>(-0.5) * (x[0] * x[0] + x[1] * x[1]);
    };
    HMC a(logp, 0.3, 10, 42);
    HMC b(logp, 0.3, 10, 42);
    std::vector<double> sa{0.5, -0.5};
    std::vector<double> sb{0.5, -0.5};
    auto ra = a.sample(sa, 100);
    auto rb = b.sample(sb, 100);
    for (std::size_t i = 0; i < ra.samples.size(); ++i) {
      assert(ra.samples[i][0] == rb.samples[i][0]);
      assert(ra.samples[i][1] == rb.samples[i][1]);
    }
  }

  // ─── 5. Bad inputs throw. ───────────────────────────────────────────
  {
    auto logp = [](const std::vector<Dual<double>>& x) {
      return Dual<double>(-0.5) * x[0] * x[0];
    };
    bool threw = false;
    try { HMC h(logp, /*step_size=*/-0.1, 10); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { HMC h(logp, 0.1, /*n_leapfrog=*/0); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    auto bad_logp = [](const std::vector<Dual<double>>& x) {
      Dual<double> result;
      result.v = (x[0].v >= 0.0)
        ? -0.5 * x[0].v * x[0].v
        : -std::numeric_limits<double>::infinity();
      result.d = 0.0;
      return result;
    };
    HMC h2(bad_logp, 0.1, 10);
    std::vector<double> bad{-1.0};
    threw = false;
    try { h2.sample(bad, 5); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  return 0;
}
