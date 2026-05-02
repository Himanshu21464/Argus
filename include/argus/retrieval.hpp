#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "mcmc.hpp"
#include "radiative_transfer.hpp"

namespace argus {

// Prior distribution shape. Hard-clipped to [prior_min, prior_max] in
// every type — the prior is zero outside the box even for Gaussian.
enum class PriorType {
  Uniform,        // log_prior = 0 inside box
  Gaussian,       // log_prior = -0.5 * ((x - mean) / stddev)^2 inside box
  LogUniform,     // log_prior = -ln(x) inside box (x > 0)
};

// Free parameter in a retrieval. Defaults to a uniform prior on
// [prior_min, prior_max] for backward compatibility.
struct Parameter {
  std::string name;
  double prior_min;
  double prior_max;
  // Prior shape inside the box. Defaults to Uniform.
  PriorType prior_type = PriorType::Uniform;
  // Used only when prior_type == Gaussian.
  double prior_mean   = 0.0;
  double prior_stddev = 1.0;
};

// Per-parameter posterior summary: median + symmetric 1-sigma credible
// interval (16th and 84th percentiles).
struct PosteriorEntry {
  std::string name;
  double median;
  double q16;       // 16th percentile (lower 1-sigma bound)
  double q84;       // 84th percentile (upper 1-sigma bound)
  double mean;
  double stddev;
};

class PosteriorSummary {
 public:
  PosteriorSummary(const std::vector<Parameter>& params,
                   const std::vector<std::vector<double>>& samples);

  const std::vector<PosteriorEntry>& entries() const noexcept {
    return entries_;
  }

  // Convenience lookup by parameter name. Throws std::out_of_range if
  // the name is not in the parameter list.
  const PosteriorEntry& operator[](const std::string& name) const;

 private:
  std::vector<PosteriorEntry> entries_;
};

// Retrieval — wraps a forward model + uniform priors + chi-squared
// likelihood + MCMC sampler.
//
// `forward(state)` returns a Spectrum given the parameter vector.
// The state vector layout matches the order of `params` passed in.
class Retrieval {
 public:
  using Forward = std::function<Spectrum(const std::vector<double>&)>;

  Retrieval(std::vector<Parameter> params,
            Forward forward,
            Spectrum observation,
            std::vector<double> uncertainty);

  // Build the full log-posterior: log_prior + log_likelihood. Uniform
  // priors return -inf outside [prior_min, prior_max].
  double log_posterior(const std::vector<double>& state) const;

  // Run an MCMC retrieval. Proposal widths default to (max - min)/40
  // per parameter — a reasonable starting heuristic. Returns the
  // (n_samples, n_dim) sample matrix and the underlying sampler so
  // callers can read its acceptance rate.
  struct Result {
    std::vector<std::vector<double>> samples;
    std::vector<double> log_posteriors;
    double acceptance_rate;
  };

  Result run_mcmc(std::vector<double> init_state,
                  std::size_t burn_in_steps,
                  std::size_t n_samples,
                  std::vector<double> proposal_widths = {},
                  std::uint64_t seed = 0) const;

  const std::vector<Parameter>& parameters() const noexcept { return params_; }

 private:
  std::vector<Parameter> params_;
  Forward forward_;
  Spectrum observation_;
  std::vector<double> uncertainty_;     // 1-sigma per wavelength
};

}  // namespace argus
