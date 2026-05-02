#include "argus/mcmc.hpp"

#include <cmath>
#include <stdexcept>

namespace argus {

MetropolisHastings::MetropolisHastings(LogPosterior log_posterior,
                                       std::vector<double> proposal_widths,
                                       std::uint64_t seed)
    : log_p_(std::move(log_posterior)),
      proposal_widths_(std::move(proposal_widths)),
      rng_(seed) {
  if (!log_p_) {
    throw std::invalid_argument("MetropolisHastings: log_posterior is null");
  }
  if (proposal_widths_.empty()) {
    throw std::invalid_argument(
        "MetropolisHastings: proposal_widths must be non-empty");
  }
  for (double w : proposal_widths_) {
    if (!(w > 0.0)) {
      throw std::invalid_argument(
          "MetropolisHastings: proposal_widths must all be positive");
    }
  }
}

bool MetropolisHastings::step(std::vector<double>& state,
                              double& current_logp) {
  const std::size_t n = state.size();
  std::vector<double> proposal(n);
  for (std::size_t i = 0; i < n; ++i) {
    std::normal_distribution<double> dist(0.0, proposal_widths_[i]);
    proposal[i] = state[i] + dist(rng_);
  }
  const double proposed_logp = log_p_(proposal);
  ++n_proposed_;
  if (!std::isfinite(proposed_logp)) {
    // Out-of-bounds prior or numerical failure — reject silently.
    return false;
  }
  // Detailed-balance acceptance ratio for symmetric Gaussian proposal:
  //     alpha = min(1, exp(proposed - current))
  const double log_alpha = proposed_logp - current_logp;
  std::uniform_real_distribution<double> u(0.0, 1.0);
  const double r = u(rng_);
  if (std::log(r) < log_alpha) {
    state = std::move(proposal);
    current_logp = proposed_logp;
    ++n_accepted_;
    return true;
  }
  return false;
}

void MetropolisHastings::burn_in(std::vector<double>& state,
                                 std::size_t n_steps) {
  if (state.size() != proposal_widths_.size()) {
    throw std::invalid_argument(
        "MetropolisHastings::burn_in: state.size() != proposal_widths.size()");
  }
  double current_logp = log_p_(state);
  if (!std::isfinite(current_logp)) {
    throw std::invalid_argument(
        "MetropolisHastings::burn_in: initial state has non-finite log-posterior");
  }
  for (std::size_t i = 0; i < n_steps; ++i) {
    (void)step(state, current_logp);
  }
}

MetropolisHastings::Result MetropolisHastings::sample(
    std::vector<double>& state, std::size_t n_samples) {
  if (state.size() != proposal_widths_.size()) {
    throw std::invalid_argument(
        "MetropolisHastings::sample: state.size() != proposal_widths.size()");
  }
  double current_logp = log_p_(state);
  if (!std::isfinite(current_logp)) {
    throw std::invalid_argument(
        "MetropolisHastings::sample: initial state has non-finite log-posterior");
  }

  Result out;
  out.samples.reserve(n_samples);
  out.log_posteriors.reserve(n_samples);
  for (std::size_t i = 0; i < n_samples; ++i) {
    (void)step(state, current_logp);
    out.samples.push_back(state);
    out.log_posteriors.push_back(current_logp);
  }
  return out;
}

double MetropolisHastings::acceptance_rate() const noexcept {
  if (n_proposed_ == 0) return 0.0;
  return static_cast<double>(n_accepted_) /
         static_cast<double>(n_proposed_);
}

}  // namespace argus
