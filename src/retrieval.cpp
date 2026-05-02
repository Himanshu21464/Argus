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
  // Uniform prior: 0 inside the box, -inf outside.
  for (std::size_t i = 0; i < params_.size(); ++i) {
    if (state[i] < params_[i].prior_min || state[i] > params_[i].prior_max) {
      return -std::numeric_limits<double>::infinity();
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
  return -0.5 * chi2;
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

}  // namespace argus
