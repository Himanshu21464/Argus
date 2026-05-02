// Argus example 04 — end-to-end retrieval.
//
// 1. Load real HITRAN H2O lines.
// 2. Generate a synthetic JWST-like transmission spectrum at known
//    "truth" T and H2O VMR.
// 3. Add Gaussian noise.
// 4. Run an MCMC retrieval with uniform priors and report the
//    recovered posterior alongside the injected truth.
//
// This is the M3 wedge — the core retrieval workflow that POSEIDON,
// CHIMERA, TauREx etc. all do, but in 100 lines of Argus C++.

#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

#include "argus/argus.hpp"
#include "argus/test_data.hpp"

int main() {
  using namespace argus;

  std::cout << "Argus " << version_string() << " — example 04 (retrieval)\n\n";

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

  std::cout << "Injected truth:  T = " << TRUE_T << " K, "
            << "log10(VMR) = " << TRUE_LV << "  (VMR = "
            << std::pow(10.0, TRUE_LV) << ")\n";
  std::cout << "Noise (1-sigma): " << noise_sigma << " transit-depth\n\n";

  // 4. Set up retrieval and run MCMC.
  std::vector<Parameter> params{
    {"T_K",       800.0,  2200.0},
    {"log10_VMR", -6.0,   -1.0},
  };
  Retrieval ret(params, [&](const std::vector<double>& s) {
                  return forward(s[0], s[1]);
                },
                observed, uncertainty);

  std::cout << "Running MCMC ...\n";
  auto result = ret.run_mcmc(/*init=*/{1300.0, -2.7},
                             /*burn_in=*/3000,
                             /*n_samples=*/5000,
                             /*proposal_widths=*/{20.0, 0.04},
                             /*seed=*/99);
  std::cout << "  acceptance = "
            << std::fixed << std::setprecision(3)
            << result.acceptance_rate * 100.0 << "%\n\n";

  // 5. Summarise the posterior.
  PosteriorSummary post(params, result.samples);
  for (const auto& e : post.entries()) {
    std::cout << "  " << std::setw(11) << std::left << e.name
              << " = " << std::fixed << std::setprecision(3)
              << std::right << std::setw(8) << e.median
              << "  +" << std::setw(7) << (e.q84 - e.median)
              << "  -" << std::setw(7) << (e.median - e.q16)
              << "  (mean " << std::setw(7) << e.mean
              << " ± " << std::setw(6) << e.stddev << ")\n";
  }
  std::cout << std::flush;

  return 0;
}
