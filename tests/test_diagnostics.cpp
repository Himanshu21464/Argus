// Tests for the Gelman-Rubin R̂ + effective-sample-size diagnostics.

#include <cassert>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

#include "argus/argus.hpp"

int main() {
  using namespace argus;

  // ─── 1. Converged chains: 4 chains all sampling from the same
  //     standard Gaussian. R̂ should be very close to 1. ──────────────
  {
    std::vector<std::vector<std::vector<double>>> chains(4);
    std::mt19937_64 rng(42);
    std::normal_distribution<double> nd(0.0, 1.0);
    for (auto& c : chains) {
      for (int i = 0; i < 2000; ++i) c.push_back({nd(rng)});
    }
    auto diag = compute_diagnostics(chains);
    assert(diag.r_hat.size() == 1);
    // R̂ very close to 1 for genuinely independent chains
    assert(diag.r_hat[0] < 1.05);
    assert(diag.r_hat[0] > 0.95);
    // ESS for IID samples ≈ N * M (independent samples)
    assert(diag.ess[0] > 0.5 * 2000.0 * 4.0);
  }

  // ─── 2. Non-converged chains: each chain stuck near a different
  //     mean. R̂ should be substantially > 1. ─────────────────────────
  {
    std::vector<std::vector<std::vector<double>>> chains(3);
    std::mt19937_64 rng(7);
    std::normal_distribution<double> tight(0.0, 0.1);
    const double offsets[] = {-2.0, 0.0, +3.0};
    for (std::size_t c = 0; c < 3; ++c) {
      for (int i = 0; i < 1000; ++i) {
        chains[c].push_back({offsets[c] + tight(rng)});
      }
    }
    auto diag = compute_diagnostics(chains);
    assert(diag.r_hat[0] > 1.5);   // far from 1 — clear non-convergence
  }

  // ─── 3. Highly autocorrelated samples: AR(1) with ρ = 0.9 means
  //     ESS ≈ N * (1-ρ)/(1+ρ) ≈ N / 19. Verify ESS << raw count. ─────
  {
    std::vector<std::vector<std::vector<double>>> chains(2);
    std::mt19937_64 rng(13);
    std::normal_distribution<double> innov(0.0, std::sqrt(1.0 - 0.9 * 0.9));
    for (auto& c : chains) {
      double x = innov(rng);
      for (int i = 0; i < 5000; ++i) {
        x = 0.9 * x + innov(rng);
        c.push_back({x});
      }
    }
    auto diag = compute_diagnostics(chains);
    // R̂ should still be ≈ 1 (chains explore same distribution)
    assert(diag.r_hat[0] < 1.10);
    // ESS should be much less than the raw sample count.
    const double raw = 5000.0 * 2.0;
    assert(diag.ess[0] < raw / 5.0);
    assert(diag.ess[0] > 0.0);
  }

  // ─── 4. Multi-dimensional case ──────────────────────────────────────
  {
    std::vector<std::vector<std::vector<double>>> chains(4);
    std::mt19937_64 rng(99);
    std::normal_distribution<double> n0(0.0, 1.0);
    std::normal_distribution<double> n1(5.0, 2.0);
    for (auto& c : chains) {
      for (int i = 0; i < 1500; ++i) {
        c.push_back({n0(rng), n1(rng)});
      }
    }
    auto diag = compute_diagnostics(chains);
    assert(diag.n_dim == 2);
    assert(diag.r_hat.size() == 2);
    assert(diag.r_hat[0] < 1.05 && diag.r_hat[1] < 1.05);
    // ESS for both should be substantial.
    assert(diag.ess[0] > 1500.0);
    assert(diag.ess[1] > 1500.0);
  }

  // ─── 5. Bad inputs throw ────────────────────────────────────────────
  {
    std::vector<std::vector<std::vector<double>>> empty;
    bool threw = false;
    try { compute_diagnostics(empty); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    threw = false;
    try {
      // < 4 samples per chain
      std::vector<std::vector<std::vector<double>>> too_short{
        {{1.0}, {2.0}, {3.0}}};
      compute_diagnostics(too_short);
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    threw = false;
    try {
      // unequal chain lengths
      std::vector<std::vector<std::vector<double>>> uneven{
        {{1.0}, {2.0}, {3.0}, {4.0}, {5.0}},
        {{1.0}, {2.0}, {3.0}, {4.0}}};
      compute_diagnostics(uneven);
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // ─── 6. Single-chain edge case: R̂ defined as 1.0 ───────────────────
  {
    std::vector<std::vector<std::vector<double>>> chains(1);
    std::mt19937_64 rng(1);
    std::normal_distribution<double> nd(0.0, 1.0);
    for (int i = 0; i < 200; ++i) chains[0].push_back({nd(rng)});
    auto diag = compute_diagnostics(chains);
    assert(diag.r_hat[0] == 1.0);
  }

  // ─── 7. Diagnostics from an EnsembleSampler::Result ─────────────────
  {
    auto logp = [](const std::vector<double>& v) {
      return -0.5 * (v[0] * v[0] + v[1] * v[1]);
    };
    std::mt19937_64 rng(2026);
    std::normal_distribution<double> jitter(0.0, 0.5);
    std::vector<std::vector<double>> walkers(20, std::vector<double>(2));
    for (auto& w : walkers) {
      w[0] = jitter(rng);
      w[1] = jitter(rng);
    }
    EnsembleSampler s(logp, walkers, 2.0, /*seed=*/1234);
    s.burn_in(500);
    auto result = s.sample(1000);
    auto diag = compute_diagnostics(result);
    assert(diag.n_chains == 20);
    assert(diag.n_steps  == 1000);
    assert(diag.n_dim    == 2);
    assert(diag.r_hat[0] < 1.10);
    assert(diag.r_hat[1] < 1.10);
    assert(diag.ess[0] > 0.0);
    assert(diag.ess[1] > 0.0);
  }

  // ─── Corrupt EnsembleSampler::Result (samples.size() inconsistent
  //     with n_steps · n_walkers) — earlier code would OOB inside the
  //     reshape loop; now it throws cleanly. ─────────────────────────
  {
    EnsembleSampler::Result r;
    r.n_walkers = 4;
    r.n_steps   = 10;
    r.n_dim     = 2;
    // Forge inconsistent samples: 10 entries instead of 4 · 10 = 40.
    for (std::size_t i = 0; i < 10; ++i) r.samples.push_back({0.0, 0.0});
    bool threw = false;
    try { (void)compute_diagnostics(r); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  return 0;
}
