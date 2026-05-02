// Argus example 04 — full M3 retrieval pipeline.
//
// 1. Load real HITRAN H2O lines.
// 2. Generate a synthetic JWST-like transmission spectrum at known
//    "truth" T and H2O VMR.
// 3. Add Gaussian noise.
// 4. Run an emcee-style ensemble MCMC retrieval.
// 5. Report acceptance rate + Gelman-Rubin R̂ + effective sample size.
// 6. Save the chain to /tmp/argus_example04.csv (load with numpy/pandas).
// 7. Posterior-predictive: report the 16/50/84% spectrum bands and how
//    many observed points fall inside.
//
// This is the kernel-ready M3 retrieval pipeline that POSEIDON,
// CHIMERA, TauREx etc. all do — but in 100 lines of Argus C++.

#include <cstdio>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

#include "argus/argus.hpp"
#include "argus/test_data.hpp"

int main() {
  using namespace argus;

  std::cout << "Argus " << version_string() << " — example 04 (retrieval)\n";
  std::cout << std::string(60, '-') << "\n";

  // 1. Load real H2O HITRAN records.
  std::istringstream is{std::string(test_data::kH2OLines)};
  auto records = Hitran::load(is, /*filter=*/1);
  std::vector<Line> lines;
  for (const auto& r : records) lines.push_back(r.line);
  auto opacity = std::make_shared<LineListOpacity>("H2O", lines, 18.015);

  // 2. JWST-like wavenumber grid through the H2O 1.4 μm band.
  std::vector<double> wn;
  for (int i = 0; i < 22; ++i) wn.push_back(6900.0 + i * 17.0);

  Species h2o{"H2O", 18.015};
  auto forward = [&](double T_k, double log10_vmr) {
    Atmosphere atm = isothermal(T_k, 1.0e-6, 1.0e2, 40, h2o,
                                std::pow(10.0, log10_vmr));
    TransmissionModel m;
    m.add_opacity(opacity);
    return m.forward(atm, wn);
  };

  // 3. Inject truth and add noise.
  const double TRUE_T  = 1500.0;
  const double TRUE_LV = -3.0;
  Spectrum truth = forward(TRUE_T, TRUE_LV);

  std::mt19937_64 rng(2026);
  const double noise_sigma = 5.0e-5;
  std::normal_distribution<double> noise(0.0, noise_sigma);
  Spectrum observed = truth;
  for (auto& v : observed.values) v += noise(rng);
  std::vector<double> uncertainty(wn.size(), noise_sigma);

  std::cout << "Truth:   T = " << TRUE_T << " K, log10(VMR) = "
            << TRUE_LV << "  (VMR = "
            << std::pow(10.0, TRUE_LV) << ")\n";
  std::cout << "Noise:   σ = " << noise_sigma << " transit-depth\n";
  std::cout << "Grid:    " << wn.size() << " wavenumbers across the H2O 1.4 μm band\n\n";

  // 4. Set up retrieval.
  std::vector<Parameter> params{
    {"T_K",       800.0,  2200.0},
    {"log10_VMR", -6.0,   -1.0},
  };
  auto wrapped = [&](const std::vector<double>& s) {
    return forward(s[0], s[1]);
  };
  Retrieval ret(params, wrapped, observed, uncertainty);

  // 4a. Run the single-chain MH baseline (faster, used for the chain).
  std::cout << "Running single-chain MH ...\n";
  auto result = ret.run_mcmc(/*init=*/{1300.0, -2.7},
                             /*burn=*/3000,
                             /*ns=*/5000,
                             /*proposal_widths=*/{20.0, 0.04},
                             /*seed=*/99);
  std::cout << "  acceptance      = "
            << std::fixed << std::setprecision(1)
            << result.acceptance_rate * 100.0 << "%\n";

  // 4b. Run the ensemble sampler in parallel for diagnostics.
  std::cout << "Running ensemble sampler (16 walkers) ...\n";
  std::mt19937_64 init_rng(2026);
  std::normal_distribution<double> init_T(1500.0, 50.0);
  std::normal_distribution<double> init_LV(-3.0, 0.2);
  std::vector<std::vector<double>> walkers(16, std::vector<double>(2));
  for (auto& w : walkers) {
    w[0] = init_T(init_rng);
    w[1] = init_LV(init_rng);
  }
  auto logp_fn = [&](const std::vector<double>& s) {
    return ret.log_posterior(s);
  };
  EnsembleSampler ens(logp_fn, walkers, /*stretch_a=*/2.0, /*seed=*/77);
  ens.burn_in(500);
  auto ens_result = ens.sample(500);
  auto diag = compute_diagnostics(ens_result);
  std::cout << "  ensemble acc.   = "
            << ens.acceptance_rate() * 100.0 << "%\n";
  for (std::size_t d = 0; d < diag.r_hat.size(); ++d) {
    std::cout << "  " << std::setw(11) << std::left << params[d].name
              << " R̂ = " << std::fixed << std::setprecision(3)
              << diag.r_hat[d]
              << "   ESS = " << std::setprecision(0) << diag.ess[d] << "\n";
  }
  std::cout << "\n";

  // 5. Posterior summary.
  PosteriorSummary post(params, result.samples);
  std::cout << "Posterior:\n";
  for (const auto& e : post.entries()) {
    std::cout << "  " << std::setw(11) << std::left << e.name
              << " = " << std::fixed << std::setprecision(3)
              << std::right << std::setw(8) << e.median
              << "  +" << std::setw(7) << (e.q84 - e.median)
              << "  -" << std::setw(7) << (e.median - e.q16) << "\n";
  }

  // 6. Save chain to CSV.
  const std::string chain_path = "/tmp/argus_example04.csv";
  ChainIO::save_csv(chain_path, params, result);
  std::cout << "\nChain saved to: " << chain_path
            << " (load with `numpy.loadtxt`, `pandas.read_csv`, etc.)\n";

  // 7. Posterior-predictive: spectrum band coverage.
  auto pp = ret.posterior_predictive(result.samples, /*thin=*/10);
  int inside = 0;
  for (std::size_t i = 0; i < wn.size(); ++i) {
    if (observed.values[i] >= pp.bands[i][0] - 2.0 * noise_sigma &&
        observed.values[i] <= pp.bands[i][2] + 2.0 * noise_sigma) {
      ++inside;
    }
  }
  std::cout << "\nPosterior-predictive 16-84% band:\n  "
            << inside << " / " << wn.size()
            << " observed points fall within (band ± 2σ_noise) — "
            << std::setprecision(0)
            << 100.0 * static_cast<double>(inside) /
               static_cast<double>(wn.size()) << "% coverage\n";

  return 0;
}
