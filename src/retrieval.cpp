#include "argus/retrieval.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace argus {

// ─── PosteriorSummary ─────────────────────────────────────────────────

namespace {

double percentile(std::vector<double> v, double q) {
  std::sort(v.begin(), v.end());
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  const double idx = q * static_cast<double>(v.size() - 1);
  const std::size_t lo = static_cast<std::size_t>(std::floor(idx));
  const std::size_t hi = std::min(lo + 1, v.size() - 1);
  const double frac = idx - static_cast<double>(lo);
  return v[lo] * (1.0 - frac) + v[hi] * frac;
}

}  // namespace

PosteriorSummary::PosteriorSummary(
    const std::vector<Parameter>& params,
    const std::vector<std::vector<double>>& samples) {
  if (samples.empty()) {
    throw std::invalid_argument("PosteriorSummary: empty sample set");
  }
  const std::size_t n_dim = params.size();
  if (samples.front().size() != n_dim) {
    throw std::invalid_argument(
        "PosteriorSummary: sample dimension does not match parameter count");
  }

  entries_.reserve(n_dim);
  for (std::size_t d = 0; d < n_dim; ++d) {
    std::vector<double> col;
    col.reserve(samples.size());
    for (const auto& s : samples) col.push_back(s[d]);

    double sum = 0.0;
    for (double x : col) sum += x;
    const double mean = sum / static_cast<double>(col.size());
    double sq = 0.0;
    for (double x : col) { const double d2 = x - mean; sq += d2 * d2; }
    const double stddev = std::sqrt(sq / static_cast<double>(col.size()));

    PosteriorEntry e;
    e.name   = params[d].name;
    e.median = percentile(col, 0.50);
    e.q16    = percentile(col, 0.16);
    e.q84    = percentile(col, 0.84);
    e.mean   = mean;
    e.stddev = stddev;
    entries_.push_back(e);
  }
}

const PosteriorEntry& PosteriorSummary::operator[](
    const std::string& name) const {
  for (const auto& e : entries_) {
    if (e.name == name) return e;
  }
  throw std::out_of_range("PosteriorSummary: no parameter named '" + name + "'");
}

// ─── Retrieval ────────────────────────────────────────────────────────

Retrieval::Retrieval(std::vector<Parameter> params,
                     Forward forward,
                     Spectrum observation,
                     std::vector<double> uncertainty)
    : params_(std::move(params)),
      forward_(std::move(forward)),
      observation_(std::move(observation)),
      uncertainty_(std::move(uncertainty)) {
  if (params_.empty()) {
    throw std::invalid_argument("Retrieval: parameter list must be non-empty");
  }
  if (!forward_) {
    throw std::invalid_argument("Retrieval: forward function is null");
  }
  if (uncertainty_.size() != observation_.values.size()) {
    throw std::invalid_argument(
        "Retrieval: uncertainty.size() must match observation.values.size()");
  }
  for (const auto& p : params_) {
    if (!(p.prior_max > p.prior_min)) {
      throw std::invalid_argument(
          "Retrieval: prior_max must exceed prior_min for parameter '" +
          p.name + "'");
    }
  }
  for (double u : uncertainty_) {
    if (!(u > 0.0)) {
      throw std::invalid_argument(
          "Retrieval: all uncertainties must be positive");
    }
  }
}

double Retrieval::log_posterior(const std::vector<double>& state) const {
  if (state.size() != params_.size()) {
    throw std::invalid_argument(
        "Retrieval::log_posterior: state.size() does not match parameter count");
  }
  // Bounds + prior contribution per parameter.
  double log_prior = 0.0;
  for (std::size_t i = 0; i < params_.size(); ++i) {
    const double x = state[i];
    const Parameter& p = params_[i];
    // All prior types hard-clip to [prior_min, prior_max].
    if (x < p.prior_min || x > p.prior_max) {
      return -std::numeric_limits<double>::infinity();
    }
    switch (p.prior_type) {
      case PriorType::Uniform:
        // log_prior = -ln(prior_max - prior_min); constant offset
        // omitted (does not affect MCMC acceptance ratios).
        break;
      case PriorType::Gaussian: {
        const double r = (x - p.prior_mean) / p.prior_stddev;
        log_prior += -0.5 * r * r;
        break;
      }
      case PriorType::LogUniform:
        // p(x) ∝ 1/x, so log_prior += -ln(x). Requires x > 0.
        if (!(x > 0.0)) {
          return -std::numeric_limits<double>::infinity();
        }
        log_prior += -std::log(x);
        break;
    }
  }
  // Chi-squared log-likelihood:
  //     ln L = -0.5 * sum_i ((obs_i - model_i) / sigma_i)^2  + const
  Spectrum model = forward_(state);
  if (model.values.size() != observation_.values.size()) {
    return -std::numeric_limits<double>::infinity();
  }
  double chi2 = 0.0;
  for (std::size_t i = 0; i < observation_.values.size(); ++i) {
    const double r = (observation_.values[i] - model.values[i]) /
                     uncertainty_[i];
    chi2 += r * r;
  }
  return log_prior - 0.5 * chi2;
}

Retrieval::PosteriorPredictive Retrieval::posterior_predictive(
    const std::vector<std::vector<double>>& samples,
    std::size_t thin,
    std::vector<double> quantiles) const {
  if (samples.empty()) {
    throw std::invalid_argument(
        "Retrieval::posterior_predictive: empty sample set");
  }
  if (thin == 0) thin = 1;
  for (double q : quantiles) {
    if (!(q >= 0.0 && q <= 1.0)) {
      throw std::invalid_argument(
          "Retrieval::posterior_predictive: quantile must be in [0, 1]");
    }
  }

  // 1. Run forward model at each thinned sample, collect per-wavelength
  //    columns of model values.
  const std::size_t n_keep = (samples.size() + thin - 1) / thin;
  std::vector<std::vector<double>> per_wn;     // [n_wn][n_keep]
  std::vector<double> wn;
  for (std::size_t s = 0; s < samples.size(); s += thin) {
    Spectrum m = forward_(samples[s]);
    if (per_wn.empty()) {
      per_wn.assign(m.values.size(), std::vector<double>{});
      for (auto& col : per_wn) col.reserve(n_keep);
      wn = m.wavenumber_cm;
    } else if (m.values.size() != per_wn.size()) {
      throw std::runtime_error(
          "Retrieval::posterior_predictive: forward model returned "
          "spectrum of inconsistent length across samples");
    }
    for (std::size_t i = 0; i < m.values.size(); ++i) {
      per_wn[i].push_back(m.values[i]);
    }
  }

  // 2. For each wavelength, sort the column then read out quantiles.
  PosteriorPredictive out;
  out.wavenumber_cm = std::move(wn);
  out.quantiles     = quantiles;
  out.bands.assign(per_wn.size(), std::vector<double>(quantiles.size()));
  for (std::size_t i = 0; i < per_wn.size(); ++i) {
    std::vector<double>& col = per_wn[i];
    std::sort(col.begin(), col.end());
    for (std::size_t q = 0; q < quantiles.size(); ++q) {
      const double idx = quantiles[q] *
                          static_cast<double>(col.size() - 1);
      const std::size_t lo =
          static_cast<std::size_t>(std::floor(idx));
      const std::size_t hi = std::min(lo + 1, col.size() - 1);
      const double frac = idx - static_cast<double>(lo);
      out.bands[i][q] = col[lo] * (1.0 - frac) + col[hi] * frac;
    }
  }
  return out;
}

Retrieval::Result Retrieval::run_mcmc(
    std::vector<double> init_state,
    std::size_t burn_in_steps,
    std::size_t n_samples,
    std::vector<double> proposal_widths,
    std::uint64_t seed) const {
  if (init_state.size() != params_.size()) {
    throw std::invalid_argument(
        "Retrieval::run_mcmc: init_state.size() does not match parameter count");
  }
  if (proposal_widths.empty()) {
    proposal_widths.resize(params_.size());
    for (std::size_t i = 0; i < params_.size(); ++i) {
      proposal_widths[i] = (params_[i].prior_max - params_[i].prior_min) / 40.0;
    }
  } else if (proposal_widths.size() != params_.size()) {
    throw std::invalid_argument(
        "Retrieval::run_mcmc: proposal_widths.size() does not match parameter count");
  }

  // Bind the log-posterior — capturing `this` is fine for the lifetime of
  // the call.
  auto logp = [this](const std::vector<double>& s) {
    return this->log_posterior(s);
  };
  MetropolisHastings sampler(logp, proposal_widths, seed);

  if (burn_in_steps > 0) sampler.burn_in(init_state, burn_in_steps);
  auto mh = sampler.sample(init_state, n_samples);

  Result out;
  out.samples         = std::move(mh.samples);
  out.log_posteriors  = std::move(mh.log_posteriors);
  out.acceptance_rate = sampler.acceptance_rate();
  return out;
}

Retrieval::AdaptiveResult Retrieval::run_mcmc_adaptive(
    std::vector<double> init_state,
    std::size_t burn_in_steps,
    std::size_t n_samples,
    std::vector<double> initial_widths,
    double target_accept,
    std::size_t adapt_interval,
    std::uint64_t seed) const {
  if (init_state.size() != params_.size()) {
    throw std::invalid_argument(
        "Retrieval::run_mcmc_adaptive: init_state.size() does not match parameter count");
  }
  if (initial_widths.empty()) {
    initial_widths.resize(params_.size());
    for (std::size_t i = 0; i < params_.size(); ++i) {
      initial_widths[i] = (params_[i].prior_max - params_[i].prior_min) / 40.0;
    }
  } else if (initial_widths.size() != params_.size()) {
    throw std::invalid_argument(
        "Retrieval::run_mcmc_adaptive: initial_widths.size() does not match parameter count");
  }

  auto logp = [this](const std::vector<double>& s) {
    return this->log_posterior(s);
  };
  MetropolisHastings sampler(logp, initial_widths, seed);

  if (burn_in_steps > 0) {
    sampler.burn_in_adaptive(init_state, burn_in_steps,
                             target_accept, adapt_interval);
  }
  // Reset the persistent counters so the reported acceptance rate
  // reflects the FROZEN-width sampling phase only — that's the rate
  // that matters for downstream mixing diagnostics.
  const auto burn_acc  = sampler.accepted();
  const auto burn_prop = sampler.proposed();
  auto mh = sampler.sample(init_state, n_samples);

  AdaptiveResult out;
  out.samples         = std::move(mh.samples);
  out.log_posteriors  = std::move(mh.log_posteriors);
  // Acceptance during the sample phase only.
  const auto sample_acc  = sampler.accepted()  - burn_acc;
  const auto sample_prop = sampler.proposed() - burn_prop;
  out.acceptance_rate = sample_prop > 0
      ? static_cast<double>(sample_acc) / static_cast<double>(sample_prop)
      : 0.0;
  out.tuned_widths    = sampler.proposal_widths();
  return out;
}

}  // namespace argus
