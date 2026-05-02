#pragma once

#include <string>
#include <vector>

#include "retrieval.hpp"

namespace argus {

// CSV reader/writer for MCMC chains. The format:
//
//   # Argus chain v0.4.x
//   # n_samples=12345 n_dim=3
//   T_K,log10_VMR,gamma_cloud,log_posterior
//   1500.123,-3.21,0.5,-12.4
//   ...
//
// One header line with the version, one comment line with shape,
// then a CSV header row of parameter names + 'log_posterior', then
// the rows. Doubles are written with hex-float precision so the
// round-trip is bit-exact.
//
// This is the offline-analysis hand-off: any `numpy.loadtxt` or
// `pandas.read_csv` call will load the chain for corner-plotting
// or further analysis.
struct LoadedChain {
  std::vector<std::string> param_names;          // [n_dim]
  std::vector<std::vector<double>> samples;       // [n_samples][n_dim]
  std::vector<double> log_posteriors;             // [n_samples]
};

class ChainIO {
 public:
  // Write samples + log-posteriors to `path`. Throws on file open failure.
  static void save_csv(const std::string& path,
                       const std::vector<Parameter>& params,
                       const std::vector<std::vector<double>>& samples,
                       const std::vector<double>& log_posteriors);

  // Convenience: write a Retrieval::Result directly.
  static void save_csv(const std::string& path,
                       const std::vector<Parameter>& params,
                       const Retrieval::Result& result);

  // Read a chain from a CSV file written by save_csv. Returns the
  // parameter names + samples + log-posteriors. Throws on parse error
  // or file open failure.
  static LoadedChain load_csv(const std::string& path);
};

}  // namespace argus
