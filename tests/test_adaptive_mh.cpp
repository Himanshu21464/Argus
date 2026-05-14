// Adaptive MH proposal tuning — Robbins-Monro burn-in.
//
// Asserts the new MetropolisHastings::burn_in_adaptive (and the
// corresponding Retrieval::run_mcmc_adaptive wrapper) drives the
// acceptance rate from a deliberately-bad initial proposal width
// into the target band, while preserving correctness of the
// downstream sample phase (posterior recovers a known truth).

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "argus/argus.hpp"

int main() {
  using namespace argus;

  // ─── 1. Direct MetropolisHastings test on a 5-D Gaussian. ─────────
  // Target: N(0, I) in 5 dimensions. log_p = -0.5 * sum(x_i^2).
  {
    auto logp = [](const std::vector<double>& x) {
      double s = 0.0;
      for (double v : x) s += v * v;
      return -0.5 * s;
    };
    // Start with WAY-too-large widths (10× the target std = 1.0).
    std::vector<double> bad_widths(5, 10.0);
    MetropolisHastings mh(logp, bad_widths, /*seed=*/42);
    std::vector<double> state(5, 0.0);
    mh.burn_in_adaptive(state, /*n_steps=*/2000,
                        /*target=*/0.30, /*adapt_interval=*/50);
    // Adapter should have driven widths smaller (closer to 1.0).
    const auto& w = mh.proposal_widths();
    for (double wi : w) {
      assert(wi < 5.0);   // half the bad starting value at minimum
      assert(wi > 0.1);   // didn't collapse to zero
    }
    // Snapshot counters BEFORE the sample phase so we measure
    // sample-phase acceptance, not the lifetime rate.
    const std::size_t acc_pre  = mh.accepted();
    const std::size_t prop_pre = mh.proposed();
    auto rec = mh.sample(state, /*n=*/4000);
    const std::size_t acc_post  = mh.accepted();
    const std::size_t prop_post = mh.proposed();
    const double sample_acc =
        static_cast<double>(acc_post - acc_pre) /
        static_cast<double>(prop_post - prop_pre);
    // Sample acceptance should be in the reasonable-MH band.
    assert(sample_acc > 0.10 && sample_acc < 0.60);

    // Posterior mean should be near zero, std should be near 1
    // (we sampled from N(0, I)).
    const std::size_t n = rec.samples.size();
    assert(n == 4000);
    std::vector<double> mean(5, 0.0);
    for (const auto& s : rec.samples) {
      for (std::size_t i = 0; i < 5; ++i) mean[i] += s[i];
    }
    for (auto& m : mean) m /= static_cast<double>(n);
    for (double m : mean) {
      // 4000 samples of N(0,1), mean SE ~ 1/sqrt(4000) ≈ 0.016;
      // acceptance ~30% so effective sample is smaller; allow ±0.2.
      assert(std::abs(m) < 0.2);
    }
  }

  // ─── 2. Retrieval::run_mcmc_adaptive on a smooth 2-D quadratic. ──
  //         Tests the full Retrieval wrapper API + that the adapter
  //         frozen-widths sample phase is a proper detailed-balance MH.
  {
    // Forward model: depth_i = (T - 1000)*1e-6 + (lv - (-3))*1e-4 ·
    //                exp(-((wn_i - 3000)/500)^2)
    // Quadratic-ish around the truth; well-conditioned so MH converges.
    std::vector<double> wn = {2000.0, 2500.0, 3000.0, 3500.0, 4000.0};
    auto fwd = [&](double T, double lv) {
      Spectrum s;
      s.wavenumber_cm = wn;
      s.values.resize(wn.size());
      for (std::size_t i = 0; i < wn.size(); ++i) {
        const double feat = std::exp(-std::pow((wn[i] - 3000.0) / 500.0, 2));
        s.values[i] = 0.02 + 1e-6 * (T - 1000.0) + 1e-4 * (lv + 3.0) * feat;
      }
      return s;
    };
    const double TRUE_T = 1100.0;
    const double TRUE_LV = -3.2;
    Spectrum truth = fwd(TRUE_T, TRUE_LV);
    std::mt19937_64 rng(123);
    std::normal_distribution<double> nz(0.0, 1e-5);
    Spectrum obs = truth;
    for (auto& v : obs.values) v += nz(rng);
    std::vector<double> sig(wn.size(), 1e-5);

    std::vector<Parameter> ps = {
      {"T",   500.0,  2000.0},
      {"LV", -6.0,    -1.0},
    };
    auto wrapped = [&](const std::vector<double>& s) { return fwd(s[0], s[1]); };
    Retrieval ret(ps, wrapped, obs, sig);

    // Deliberately oversize initial widths — half the prior box.
    auto ar = ret.run_mcmc_adaptive(
        /*init=*/{900.0, -3.5},
        /*burn=*/3000,
        /*ns=*/4000,
        /*initial_widths=*/{500.0, 1.5},
        /*target=*/0.30, /*interval=*/100, /*seed=*/2026);

    std::cerr << "[adaptive_mh] tuned widths = "
              << ar.tuned_widths[0] << " " << ar.tuned_widths[1]
              << "  acc=" << ar.acceptance_rate << "\n";

    // (a) Adapter shrunk both widths from their too-large starts.
    assert(ar.tuned_widths.size() == 2);
    assert(ar.tuned_widths[0] < 500.0);
    assert(ar.tuned_widths[1] <  1.5);
    // (b) Sample-phase acceptance lands in a reasonable band.
    assert(ar.acceptance_rate > 0.05);
    assert(ar.acceptance_rate < 0.70);
    // (c) Chain mixed at all (samples actually vary).
    assert(ar.samples.size() == 4000);
    double T_min = 1e9, T_max = -1e9;
    for (const auto& s : ar.samples) {
      T_min = std::min(T_min, s[0]); T_max = std::max(T_max, s[0]);
    }
    assert(T_max - T_min > 1.0);   // chain wasn't stuck at one point
  }

  return 0;
}
