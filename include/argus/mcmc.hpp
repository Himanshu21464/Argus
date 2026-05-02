#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dual.hpp"

namespace argus {

// Single-chain Metropolis-Hastings sampler with isotropic Gaussian
// proposal in the bounded parameter space.
//
// The log-posterior is supplied as a std::function — the caller wires
// in priors and likelihood. The sampler is deterministic given a seed,
// so MCMC retrieval results are bit-reproducible across runs (tested
// in tests/test_mcmc.cpp).
//
// M3.5 will add ensemble samplers (emcee-style) and Hamiltonian MC
// (NUTS); M3 ships the Gaussian-proposal baseline.
class MetropolisHastings {
 public:
  using LogPosterior = std::function<double(const std::vector<double>&)>;

  // proposal_widths[i] is the standard deviation of the Gaussian step
  // along parameter i. Caller is responsible for tuning these so the
  // acceptance rate falls in a sensible range (~20-50% for a single
  // chain, lower for many parameters).
  MetropolisHastings(LogPosterior log_posterior,
                     std::vector<double> proposal_widths,
                     std::uint64_t seed = 0);

  // Run `n_steps` of MH starting from `state` (modified in place to the
  // final state). Burn-in samples are NOT returned; use this to settle
  // the chain before calling sample().
  void burn_in(std::vector<double>& state, std::size_t n_steps);

  // Run `n_samples` of MH starting from `state` (modified in place).
  // Returns the samples as a [n_samples, n_dim] flat row-major vector
  // alongside the log-posterior at each sample.
  struct Result {
    std::vector<std::vector<double>> samples;   // [n_samples][n_dim]
    std::vector<double> log_posteriors;          // [n_samples]
  };
  Result sample(std::vector<double>& state, std::size_t n_samples);

  // Diagnostics: fraction of proposals accepted across all completed
  // calls (burn_in + sample).
  double acceptance_rate() const noexcept;

  std::size_t accepted() const noexcept { return n_accepted_; }
  std::size_t proposed() const noexcept { return n_proposed_; }

 private:
  bool step(std::vector<double>& state, double& current_logp);

  LogPosterior log_p_;
  std::vector<double> proposal_widths_;
  std::mt19937_64 rng_;
  std::size_t n_accepted_ = 0;
  std::size_t n_proposed_ = 0;
};

// Compute ∇f(x) for a templated scalar function f via forward-mode
// autograd. f must be callable as f(std::vector<Dual<double>>) and
// return Dual<double>; we evaluate it D times (one per coordinate)
// with the seed derivative on that coordinate set to 1.
//
// Cost: D evaluations of f. Suitable for retrieval problems with
// O(10) parameters. M3.5 will add reverse-mode autograd for larger
// parameter spaces (one evaluation, one backward pass).
template <typename F>
inline std::vector<double> grad(F&& f, const std::vector<double>& x) {
  const std::size_t D = x.size();
  std::vector<double> g(D);
  for (std::size_t i = 0; i < D; ++i) {
    std::vector<Dual<double>> dx(D);
    for (std::size_t j = 0; j < D; ++j) {
      dx[j] = Dual<double>{x[j], (i == j) ? 1.0 : 0.0};
    }
    Dual<double> y = f(dx);
    g[i] = y.d;
  }
  return g;
}

// Hamiltonian Monte Carlo with a leapfrog integrator. Uses forward-
// mode autograd to compute ∇log p, so it works on any user-supplied
// templated log-posterior callable.
//
// Compared to MetropolisHastings, HMC mixes much better on curved /
// correlated posteriors at the cost of D + 1 forward evaluations
// per leapfrog step (D for the gradient + 1 for the Metropolis
// correction).
//
// The user supplies a templated log-posterior:
//   template <typename T> T log_p(const std::vector<T>& x)
template <typename LogP>
class HMC {
 public:
  HMC(LogP log_p,
      double step_size,
      std::size_t n_leapfrog,
      std::uint64_t seed = 0)
      : log_p_(std::move(log_p)),
        step_size_(step_size),
        n_leapfrog_(n_leapfrog),
        rng_(seed) {
    if (!(step_size_ > 0.0)) {
      throw std::invalid_argument("HMC: step_size must be positive");
    }
    if (n_leapfrog_ == 0) {
      throw std::invalid_argument("HMC: n_leapfrog must be >= 1");
    }
  }

  struct Result {
    std::vector<std::vector<double>> samples;
    std::vector<double> log_posteriors;
  };

  Result sample(std::vector<double>& state, std::size_t n_samples) {
    auto logp_double = [&](const std::vector<double>& x) {
      // Wrap x in Dual with d=0 to evaluate as plain double via the
      // templated path.
      std::vector<Dual<double>> dx(x.size());
      for (std::size_t i = 0; i < x.size(); ++i) dx[i] = Dual<double>{x[i], 0.0};
      return log_p_(dx).v;
    };

    auto grad_logp = [&](const std::vector<double>& x) {
      return grad([&](const std::vector<Dual<double>>& xd) {
        return log_p_(xd);
      }, x);
    };

    Result out;
    out.samples.reserve(n_samples);
    out.log_posteriors.reserve(n_samples);

    std::normal_distribution<double> nd(0.0, 1.0);
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    double current_lp = logp_double(state);
    if (!std::isfinite(current_lp)) {
      throw std::invalid_argument(
          "HMC::sample: initial state has non-finite log-posterior");
    }

    for (std::size_t s = 0; s < n_samples; ++s) {
      // 1. Sample momentum p ~ N(0, I).
      std::vector<double> p(state.size());
      for (auto& pi : p) pi = nd(rng_);

      // Initial Hamiltonian: H0 = -log_p(q) + 0.5 * |p|^2
      double K0 = 0.0;
      for (double pi : p) K0 += pi * pi;
      K0 *= 0.5;
      const double H0 = -current_lp + K0;

      // 2. Leapfrog: alternate q and p updates.
      std::vector<double> q = state;
      std::vector<double> g = grad_logp(q);

      // Half-step on p
      for (std::size_t i = 0; i < p.size(); ++i) {
        p[i] += 0.5 * step_size_ * g[i];
      }
      for (std::size_t L = 0; L < n_leapfrog_; ++L) {
        // Full-step on q
        for (std::size_t i = 0; i < q.size(); ++i) {
          q[i] += step_size_ * p[i];
        }
        g = grad_logp(q);
        // Full-step on p, except last is half-step
        if (L + 1 < n_leapfrog_) {
          for (std::size_t i = 0; i < p.size(); ++i) {
            p[i] += step_size_ * g[i];
          }
        } else {
          for (std::size_t i = 0; i < p.size(); ++i) {
            p[i] += 0.5 * step_size_ * g[i];
          }
        }
      }
      // Negate momentum to make proposal symmetric (cosmetic; H is
      // unchanged because K depends only on |p|^2).
      for (auto& pi : p) pi = -pi;

      // 3. Metropolis correction.
      const double new_lp = logp_double(q);
      double K1 = 0.0;
      for (double pi : p) K1 += pi * pi;
      K1 *= 0.5;
      const double H1 = -new_lp + K1;

      ++n_proposed_;
      if (std::isfinite(new_lp) &&
          std::log(u01(rng_)) < (H0 - H1)) {
        state = std::move(q);
        current_lp = new_lp;
        ++n_accepted_;
      }

      out.samples.push_back(state);
      out.log_posteriors.push_back(current_lp);
    }
    return out;
  }

  double acceptance_rate() const noexcept {
    if (n_proposed_ == 0) return 0.0;
    return static_cast<double>(n_accepted_) /
           static_cast<double>(n_proposed_);
  }

 private:
  LogP log_p_;
  double step_size_;
  std::size_t n_leapfrog_;
  std::mt19937_64 rng_;
  std::size_t n_accepted_ = 0;
  std::size_t n_proposed_ = 0;
};

// Affine-invariant ensemble sampler — the Goodman & Weare (2010)
// "stretch move" algorithm. Industry standard in astronomy via emcee
// (Foreman-Mackey et al. 2013).
//
// Key properties vs single-chain Metropolis-Hastings:
//   * affine-invariant: performance does not depend on the
//     condition number of the posterior covariance — handles tightly
//     correlated parameters as well as uncorrelated ones
//   * minimal tuning: only one hyperparameter `a` (the stretch
//     scale, default 2.0) plus the number of walkers
//   * embarrassingly parallel across the half-ensemble update
//
// API mirrors MetropolisHastings: provide a log-posterior callable
// and an initial walker ensemble; sampler returns samples shaped
// [n_steps, n_walkers, n_dim].
class EnsembleSampler {
 public:
  using LogPosterior = std::function<double(const std::vector<double>&)>;

  // walkers0 must contain `n_walkers` initial states, all of which must
  // satisfy log_posterior(walker) > -inf. n_walkers must be >= 2*n_dim
  // and even.
  EnsembleSampler(LogPosterior log_posterior,
                  std::vector<std::vector<double>> walkers0,
                  double stretch_a = 2.0,
                  std::uint64_t seed = 0);

  // Run `n_steps` of the stretch move, discarding samples (burn-in).
  // Final walker positions are kept for the next call.
  void burn_in(std::size_t n_steps);

  // Run `n_steps` and record every walker position at every step.
  // Returns flat samples [n_steps * n_walkers][n_dim] and
  // log_posteriors [n_steps * n_walkers], in step-major order
  // (walker varies fastest).
  struct Result {
    std::vector<std::vector<double>> samples;
    std::vector<double> log_posteriors;
    std::size_t n_steps;
    std::size_t n_walkers;
    std::size_t n_dim;
  };
  Result sample(std::size_t n_steps);

  double acceptance_rate() const noexcept;
  std::size_t n_walkers() const noexcept { return walkers_.size(); }
  std::size_t n_dim()     const noexcept {
    return walkers_.empty() ? 0 : walkers_.front().size();
  }

 private:
  // Half-step: update walkers in [start, end) by drawing complement
  // partners from the OTHER half [other_start, other_end). Returns
  // number of proposals and acceptances added during this half-step.
  void half_step(std::size_t start, std::size_t end,
                 std::size_t other_start, std::size_t other_end);

  LogPosterior log_p_;
  std::vector<std::vector<double>> walkers_;       // [n_walkers][n_dim]
  std::vector<double> walker_logp_;                // [n_walkers]
  double a_;                                        // stretch scale
  std::mt19937_64 rng_;
  std::size_t n_accepted_ = 0;
  std::size_t n_proposed_ = 0;
};

}  // namespace argus
