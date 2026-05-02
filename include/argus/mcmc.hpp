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

}  // namespace argus
