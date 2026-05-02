#pragma once

#include <cstddef>
#include <vector>

#include "mcmc.hpp"

namespace argus {

// Per-parameter MCMC convergence diagnostics. Returned by
// compute_diagnostics() for a multi-chain run.
struct ChainDiagnostics {
  std::vector<double> r_hat;       // Gelman-Rubin R̂ per parameter
  std::vector<double> ess;         // effective sample size per parameter
  std::size_t n_chains = 0;
  std::size_t n_steps  = 0;
  std::size_t n_dim    = 0;
};

// Gelman-Rubin R̂ + effective sample size from M chains, each of
// length N samples in D dimensions.
//
// Input layout: chains[chain_idx][step_idx][dim_idx]
//
// R̂ formula (Gelman & Rubin 1992):
//   B = (N / (M - 1)) · Σ_m (θ̄_m - θ̄_·)²
//   W = (1 / M) · Σ_m s_m²
//   V̂ = ((N - 1)/N) · W + (1/N) · B
//   R̂ = sqrt(V̂ / W)
// R̂ → 1 as chains converge to the same distribution.
//
// ESS computed from the autocorrelation of the pooled samples:
//   ESS = N · M / (1 + 2 · Σ_t ρ_t)
// where ρ_t is the lag-t autocorrelation, summed until it first
// drops below 0.05 ("initial monotone sequence" cutoff, Geyer 1992).
ChainDiagnostics compute_diagnostics(
    const std::vector<std::vector<std::vector<double>>>& chains);

// Convenience: compute diagnostics from an EnsembleSampler::Result by
// treating each walker as an independent chain.
ChainDiagnostics compute_diagnostics(const EnsembleSampler::Result& result);

}  // namespace argus
