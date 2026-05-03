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
  double prior_min = 0.0;       // Retrieval ctor validates max > min;
  double prior_max = 0.0;       // defaults make a default-constructed
                                // Parameter detect-ably "empty".
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
  double median = 0.0;
  double q16    = 0.0;       // 16th percentile (lower 1-sigma bound)
  double q84    = 0.0;       // 84th percentile (upper 1-sigma bound)
  double mean   = 0.0;
  double stddev = 0.0;
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
    double acceptance_rate = 0.0;
  };

  Result run_mcmc(std::vector<double> init_state,
                  std::size_t burn_in_steps,
                  std::size_t n_samples,
                  std::vector<double> proposal_widths = {},
                  std::uint64_t seed = 0) const;

  const std::vector<Parameter>& parameters() const noexcept { return params_; }

  // Run the forward model at each (thinned) posterior sample to build
  // a posterior-predictive spectrum: per-wavelength quantiles (16th,
  // 50th, 84th percentile by default) of the model values across the
  // sample set. The standard post-MCMC sanity check — the observed
  // spectrum should fall inside the 16-84% band at most wavelengths
  // if the retrieval converged well.
  //
  //   thin: keep one sample every `thin` (e.g. thin=10 means 10x
  //         fewer forward calls). Pass thin=1 for no thinning.
  //
  // Returns a [n_wn][n_quantile] matrix in the order requested by
  // `quantiles` (default {0.16, 0.50, 0.84}).
  struct PosteriorPredictive {
    std::vector<double> wavenumber_cm;
    std::vector<double> quantiles;                 // requested percentiles
    std::vector<std::vector<double>> bands;        // [n_wn][n_quant]
  };
  PosteriorPredictive posterior_predictive(
      const std::vector<std::vector<double>>& samples,
      std::size_t thin = 1,
      std::vector<double> quantiles = {0.16, 0.50, 0.84}) const;

 private:
  std::vector<Parameter> params_;
  Forward forward_;
  Spectrum observation_;
  std::vector<double> uncertainty_;     // 1-sigma per wavelength
};

}  // namespace argus
