#include "argus/diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace argus {

namespace {

// Sample mean of a 1-D series.
double series_mean(const std::vector<double>& v) {
  double s = 0.0;
  for (double x : v) s += x;
  return s / static_cast<double>(v.size());
}

// Sample variance with N-1 denominator (Bessel-corrected).
double series_var(const std::vector<double>& v, double mean) {
  double s = 0.0;
  for (double x : v) { const double d = x - mean; s += d * d; }
  return s / static_cast<double>(v.size() - 1);
}

// Lag-t autocorrelation of a 1-D series, normalised so ρ(0) = 1.
double autocorr(const std::vector<double>& v, double mean, double var,
                std::size_t t) {
  const std::size_t n = v.size();
  if (t >= n) return 0.0;
  double s = 0.0;
  for (std::size_t i = 0; i < n - t; ++i) {
    s += (v[i] - mean) * (v[i + t] - mean);
  }
  return (s / static_cast<double>(n - t)) / var;
}

}  // namespace

ChainDiagnostics compute_diagnostics(
    const std::vector<std::vector<std::vector<double>>>& chains) {
  if (chains.empty()) {
    throw std::invalid_argument(
        "compute_diagnostics: need at least 1 chain");
  }
  const std::size_t M = chains.size();
  const std::size_t N = chains.front().size();
  if (N < 4) {
    throw std::invalid_argument(
        "compute_diagnostics: each chain must have >= 4 samples");
  }
  const std::size_t D = chains.front().front().size();
  for (const auto& c : chains) {
    if (c.size() != N) {
      throw std::invalid_argument(
          "compute_diagnostics: all chains must have equal length");
    }
    for (const auto& s : c) {
      if (s.size() != D) {
        throw std::invalid_argument(
            "compute_diagnostics: all samples must have equal dimension");
      }
    }
  }

  ChainDiagnostics out;
  out.n_chains = M;
  out.n_steps  = N;
  out.n_dim    = D;
  out.r_hat.assign(D, 0.0);
  out.ess.assign(D, 0.0);

  for (std::size_t d = 0; d < D; ++d) {
    // Extract per-chain 1-D series.
    std::vector<std::vector<double>> series(M);
    for (std::size_t m = 0; m < M; ++m) {
      series[m].reserve(N);
      for (std::size_t n = 0; n < N; ++n) {
        series[m].push_back(chains[m][n][d]);
      }
    }
    // Per-chain mean, variance.
    std::vector<double> chain_mean(M), chain_var(M);
    for (std::size_t m = 0; m < M; ++m) {
      chain_mean[m] = series_mean(series[m]);
      chain_var[m]  = series_var(series[m], chain_mean[m]);
    }
    // R-hat (Gelman & Rubin 1992).
    double grand = 0.0;
    for (double cm : chain_mean) grand += cm;
    grand /= static_cast<double>(M);

    if (M >= 2) {
      double B = 0.0;
      for (double cm : chain_mean) {
        const double diff = cm - grand;
        B += diff * diff;
      }
      B *= static_cast<double>(N) / static_cast<double>(M - 1);
      double W = 0.0;
      for (double cv : chain_var) W += cv;
      W /= static_cast<double>(M);
      if (W <= 0.0) {
        out.r_hat[d] = std::numeric_limits<double>::quiet_NaN();
      } else {
        const double V_hat =
            (static_cast<double>(N - 1) / static_cast<double>(N)) * W
            + (1.0 / static_cast<double>(N)) * B;
        out.r_hat[d] = std::sqrt(V_hat / W);
      }
    } else {
      // Single chain → R̂ undefined; report 1.0.
      out.r_hat[d] = 1.0;
    }

    // ESS: pool all chains into one long series for autocorrelation.
    std::vector<double> pooled;
    pooled.reserve(M * N);
    for (const auto& s : series) {
      pooled.insert(pooled.end(), s.begin(), s.end());
    }
    const double pmean = series_mean(pooled);
    const double pvar  = series_var(pooled, pmean);
    if (pvar <= 0.0) {
      out.ess[d] = 0.0;
      continue;
    }
    // Sum autocorrelations until they first drop below 0.05 (initial
    // positive-sequence cutoff).
    double tau = 1.0;          // 1 + 2 * Σ ρ_t with ρ_0 = 1
    const std::size_t max_lag = std::min<std::size_t>(pooled.size() / 4,
                                                      1000);
    for (std::size_t t = 1; t < max_lag; ++t) {
      const double rho = autocorr(pooled, pmean, pvar, t);
      if (rho < 0.05) break;
      tau += 2.0 * rho;
    }
    out.ess[d] = static_cast<double>(pooled.size()) / tau;
  }
  return out;
}

ChainDiagnostics compute_diagnostics(
    const EnsembleSampler::Result& result) {
  // Reshape ensemble result [n_steps * n_walkers] into chains[walker][step].
  const std::size_t W = result.n_walkers;
  const std::size_t S = result.n_steps;
  const std::size_t D = result.n_dim;
  std::vector<std::vector<std::vector<double>>> chains(W);
  for (std::size_t w = 0; w < W; ++w) {
    chains[w].reserve(S);
    for (std::size_t s = 0; s < S; ++s) {
      chains[w].push_back(result.samples[s * W + w]);
    }
  }
  (void)D;
  return compute_diagnostics(chains);
}

}  // namespace argus
