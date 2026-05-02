// Tests for the CSV chain I/O round-trip.
//
// Hex-float encoding makes the round-trip bit-exact, which is the
// strictest reasonable test.

#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <stdexcept>
#include <vector>

#include "argus/argus.hpp"

int main() {
  using namespace argus;

  // ─── 1. Synthesise a tiny chain, write, read, compare bit-exact. ────
  {
    const std::vector<Parameter> params{
      {"T_K",         800.0,   2200.0},
      {"log10_VMR",  -6.0,    -1.0},
      {"gamma_cloud", 0.0,     1.0},
    };
    std::mt19937_64 rng(2026);
    std::normal_distribution<double> n(0.0, 1.0);
    std::vector<std::vector<double>> samples;
    std::vector<double> logp;
    for (int i = 0; i < 1000; ++i) {
      samples.push_back({1500.0 + 50.0 * n(rng),
                         -3.0   + 0.3  * n(rng),
                          0.5   + 0.1  * n(rng)});
      logp.push_back(-12.4 + 1.5 * n(rng));
    }

    const std::string path = "/tmp/argus_chain_test.csv";
    ChainIO::save_csv(path, params, samples, logp);

    LoadedChain loaded = ChainIO::load_csv(path);
    assert(loaded.param_names.size() == 3);
    assert(loaded.param_names[0] == "T_K");
    assert(loaded.param_names[1] == "log10_VMR");
    assert(loaded.param_names[2] == "gamma_cloud");
    assert(loaded.samples.size() == samples.size());
    assert(loaded.log_posteriors.size() == logp.size());

    for (std::size_t i = 0; i < samples.size(); ++i) {
      for (std::size_t d = 0; d < 3; ++d) {
        // Bit-exact via hex-float encoding
        assert(loaded.samples[i][d] == samples[i][d]);
      }
      assert(loaded.log_posteriors[i] == logp[i]);
    }
    std::remove(path.c_str());
  }

  // ─── 2. Round-trip via Retrieval::Result overload. ──────────────────
  {
    const std::vector<Parameter> params{
      {"x", -10.0, 10.0},
      {"y", -10.0, 10.0},
    };
    Retrieval::Result r;
    r.acceptance_rate = 0.42;
    for (int i = 0; i < 50; ++i) {
      r.samples.push_back({static_cast<double>(i) * 0.1,
                           static_cast<double>(i) * 0.05});
      r.log_posteriors.push_back(-static_cast<double>(i));
    }

    const std::string path = "/tmp/argus_chain_test2.csv";
    ChainIO::save_csv(path, params, r);
    LoadedChain loaded = ChainIO::load_csv(path);
    assert(loaded.samples.size() == 50);
    for (std::size_t i = 0; i < 50; ++i) {
      assert(loaded.samples[i][0] == r.samples[i][0]);
      assert(loaded.samples[i][1] == r.samples[i][1]);
      assert(loaded.log_posteriors[i] == r.log_posteriors[i]);
    }
    std::remove(path.c_str());
  }

  // ─── 3. Empty samples is allowed (still writes header). ─────────────
  {
    const std::vector<Parameter> params{{"x", 0.0, 1.0}};
    const std::string path = "/tmp/argus_chain_empty.csv";
    ChainIO::save_csv(path, params, {}, {});
    LoadedChain loaded = ChainIO::load_csv(path);
    assert(loaded.samples.empty());
    assert(loaded.log_posteriors.empty());
    assert(loaded.param_names.size() == 1);
    assert(loaded.param_names[0] == "x");
    std::remove(path.c_str());
  }

  // ─── 4. Mismatched samples / logp throws on save. ───────────────────
  {
    const std::vector<Parameter> params{{"x", 0.0, 1.0}};
    bool threw = false;
    try {
      ChainIO::save_csv("/tmp/argus_should_not_open.csv",
                        params, {{0.5}, {0.7}}, {1.0});
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // ─── 5. Sample dim != params count throws on save. ──────────────────
  {
    const std::vector<Parameter> params{{"x", 0.0, 1.0}};
    bool threw = false;
    try {
      ChainIO::save_csv("/tmp/argus_should_not_open2.csv",
                        params, {{0.5, 0.6}}, {1.0});
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // ─── 6. Missing file throws on load. ────────────────────────────────
  {
    bool threw = false;
    try { (void)ChainIO::load_csv("/tmp/argus_does_not_exist_xyz.csv"); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
  }

  // ─── 7. Malformed CSV throws on load. ───────────────────────────────
  {
    const std::string path = "/tmp/argus_chain_bad.csv";
    {
      std::ofstream out(path);
      out << "# this file has no real header\n"
          << "# just comments\n";
    }
    bool threw = false;
    try { (void)ChainIO::load_csv(path); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    std::remove(path.c_str());

    // Header without log_posterior column
    {
      std::ofstream out(path);
      out << "x,y\n0.5,0.6\n";
    }
    threw = false;
    try { (void)ChainIO::load_csv(path); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    std::remove(path.c_str());

    // Wrong column count in data row
    {
      std::ofstream out(path);
      out << "x,log_posterior\n0.5,1.0,99.0\n";
    }
    threw = false;
    try { (void)ChainIO::load_csv(path); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    std::remove(path.c_str());
  }

  // ─── 8. End-to-end: save a real Retrieval result, load it back,
  //     verify PosteriorSummary matches the original chain. ───────────
  {
    // Tiny synthetic 1-parameter retrieval.
    std::vector<Parameter> params{{"x", -5.0, 5.0}};
    Retrieval::Result r;
    r.acceptance_rate = 0.3;
    std::mt19937_64 rng(7);
    std::normal_distribution<double> nd(1.5, 0.7);
    for (int i = 0; i < 5000; ++i) {
      r.samples.push_back({nd(rng)});
      r.log_posteriors.push_back(-0.5 * std::pow((r.samples.back()[0] - 1.5) / 0.7, 2));
    }
    const std::string path = "/tmp/argus_chain_e2e.csv";
    ChainIO::save_csv(path, params, r);
    LoadedChain loaded = ChainIO::load_csv(path);

    PosteriorSummary post_orig(params, r.samples);
    PosteriorSummary post_load(params, loaded.samples);
    assert(post_orig["x"].mean   == post_load["x"].mean);
    assert(post_orig["x"].stddev == post_load["x"].stddev);
    assert(post_orig["x"].median == post_load["x"].median);
    std::remove(path.c_str());
  }

  return 0;
}
