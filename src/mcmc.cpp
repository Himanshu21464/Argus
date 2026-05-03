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

// ─── EnsembleSampler — Goodman & Weare 2010 stretch move ──────────────

EnsembleSampler::EnsembleSampler(LogPosterior log_posterior,
                                 std::vector<std::vector<double>> walkers0,
                                 double stretch_a,
                                 std::uint64_t seed)
    : log_p_(std::move(log_posterior)),
      walkers_(std::move(walkers0)),
      a_(stretch_a),
      rng_(seed) {
  if (!log_p_) {
    throw std::invalid_argument("EnsembleSampler: log_posterior is null");
  }
  if (walkers_.size() < 4) {
    throw std::invalid_argument(
        "EnsembleSampler: need at least 4 walkers");
  }
  if (walkers_.size() % 2 != 0) {
    throw std::invalid_argument("EnsembleSampler: n_walkers must be even");
  }
  const std::size_t d = walkers_.front().size();
  if (d == 0) {
    throw std::invalid_argument("EnsembleSampler: walkers must be non-empty");
  }
  for (const auto& w : walkers_) {
    if (w.size() != d) {
      throw std::invalid_argument(
          "EnsembleSampler: all walkers must have the same dimensionality");
    }
  }
  if (!(a_ > 1.0)) {
    throw std::invalid_argument("EnsembleSampler: stretch_a must be > 1");
  }

  walker_logp_.resize(walkers_.size());
  for (std::size_t i = 0; i < walkers_.size(); ++i) {
    walker_logp_[i] = log_p_(walkers_[i]);
    if (!std::isfinite(walker_logp_[i])) {
      throw std::invalid_argument(
          "EnsembleSampler: walker " + std::to_string(i) +
          " has non-finite log_posterior at init");
    }
  }
}

void EnsembleSampler::half_step(std::size_t start, std::size_t end,
                                 std::size_t other_start,
                                 std::size_t other_end) {
  // Goodman-Weare 2010 stretch move. Sample stretch factor `s ∈ [1/a, a]`
  // from g(s) ∝ 1/√s; in the inverse-CDF parameterisation we sample
  // u ∈ [0, 1] and compute u_to_sqrt_s, then the actual stretch is
  // `s = (u_to_sqrt_s)²`. This naming is verbose but unambiguous —
  // earlier code used `z` for the sqrt and `zz` for the stretch,
  // which collided with the convention in the reference paper.
  const std::size_t n_dim = walkers_.front().size();
  const double sqrt_a     = std::sqrt(a_);
  const double inv_sqrt_a = 1.0 / sqrt_a;
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  const std::size_t n_other = other_end - other_start;
  std::uniform_int_distribution<std::size_t> partner_dist(0, n_other - 1);

  for (std::size_t i = start; i < end; ++i) {
    const double u            = u01(rng_);
    const double sqrt_stretch = (sqrt_a - inv_sqrt_a) * u + inv_sqrt_a;
    const double stretch      = sqrt_stretch * sqrt_stretch;

    const std::size_t j = other_start + partner_dist(rng_);
    std::vector<double> proposal(n_dim);
    for (std::size_t d = 0; d < n_dim; ++d) {
      proposal[d] = walkers_[j][d] +
                    stretch * (walkers_[i][d] - walkers_[j][d]);
    }

    const double new_logp = log_p_(proposal);
    ++n_proposed_;
    if (!std::isfinite(new_logp)) continue;

    // Acceptance: log α = (n_dim - 1) · log(stretch) + new - old.
    const double log_alpha =
        (static_cast<double>(n_dim) - 1.0) * std::log(stretch)
        + new_logp - walker_logp_[i];
    const double r = u01(rng_);
    if (std::log(r) < log_alpha) {
      walkers_[i] = std::move(proposal);
      walker_logp_[i] = new_logp;
      ++n_accepted_;
    }
  }
}

void EnsembleSampler::burn_in(std::size_t n_steps) {
  const std::size_t n = walkers_.size();
  const std::size_t half = n / 2;
  for (std::size_t step = 0; step < n_steps; ++step) {
    half_step(0,    half, half, n);
    half_step(half, n,    0,    half);
  }
}

EnsembleSampler::Result EnsembleSampler::sample(std::size_t n_steps) {
  const std::size_t n = walkers_.size();
  const std::size_t half = n / 2;
  const std::size_t n_dim = walkers_.front().size();

  Result out;
  out.n_steps   = n_steps;
  out.n_walkers = n;
  out.n_dim     = n_dim;
  out.samples.reserve(n_steps * n);
  out.log_posteriors.reserve(n_steps * n);

  for (std::size_t step = 0; step < n_steps; ++step) {
    half_step(0,    half, half, n);
    half_step(half, n,    0,    half);
    for (std::size_t i = 0; i < n; ++i) {
      out.samples.push_back(walkers_[i]);
      out.log_posteriors.push_back(walker_logp_[i]);
    }
  }
  return out;
}

double EnsembleSampler::acceptance_rate() const noexcept {
  if (n_proposed_ == 0) return 0.0;
  return static_cast<double>(n_accepted_) /
         static_cast<double>(n_proposed_);
}

}  // namespace argus
