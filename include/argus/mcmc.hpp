#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>
#include <vector>

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
