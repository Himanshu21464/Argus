// Tests for the PosteriorSummary class.

#include <cassert>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

#include "argus/argus.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

}  // namespace

int main() {
  using namespace argus;

  // ─── 1. Single-parameter Gaussian samples. ──────────────────────────
  // Generate 100k Gaussian draws with mean=2.0, stddev=0.5 — verify
  // PosteriorSummary recovers mean/std/median accurately.
  {
    std::mt19937_64 rng(42);
    std::normal_distribution<double> nd(2.0, 0.5);
    std::vector<std::vector<double>> samples;
    samples.reserve(100000);
    for (int i = 0; i < 100000; ++i) samples.push_back({nd(rng)});

    std::vector<Parameter> params{{"x", -10.0, 10.0}};
    PosteriorSummary p(params, samples);
    const auto& e = p["x"];
    assert(close(e.mean,   2.0, 0.01));
    assert(close(e.stddev, 0.5, 0.02));
    assert(close(e.median, 2.0, 0.01));
    // 16-84 percentile width should be ~1 sigma -> ~0.5
    assert(close(e.q84 - e.q16, 1.0, 0.05));
    // 16th percentile is mean - sigma
    assert(close(e.q16, 2.0 - 0.5, 0.02));
    assert(close(e.q84, 2.0 + 0.5, 0.02));
  }

  // ─── 2. Two-parameter Gaussian, accessing entries by name. ──────────
  {
    std::mt19937_64 rng(42);
    std::normal_distribution<double> n1(0.0, 1.0);
    std::normal_distribution<double> n2(5.0, 2.0);
    std::vector<std::vector<double>> samples;
    for (int i = 0; i < 50000; ++i) samples.push_back({n1(rng), n2(rng)});

    std::vector<Parameter> params{{"a", -5.0, 5.0}, {"b", -10.0, 20.0}};
    PosteriorSummary p(params, samples);
    // Use absolute tolerance for mean(a)=0 since relative tol degenerates.
    assert(std::fabs(p["a"].mean - 0.0) < 0.05);
    assert(close(p["a"].stddev, 1.0, 0.02));
    assert(close(p["b"].mean,   5.0, 0.05));
    assert(close(p["b"].stddev, 2.0, 0.02));
  }

  // ─── 3. Lookup by unknown name throws. ──────────────────────────────
  {
    std::vector<std::vector<double>> samples{{1.0}, {2.0}, {3.0}};
    std::vector<Parameter> params{{"x", 0.0, 5.0}};
    PosteriorSummary p(params, samples);
    bool threw = false;
    try { (void)p["nonexistent"]; }
    catch (const std::out_of_range&) { threw = true; }
    assert(threw);
  }

  // ─── 4. Empty samples or shape mismatch throws. ─────────────────────
  {
    std::vector<Parameter> params{{"x", 0.0, 1.0}};
    bool threw = false;
    try { PosteriorSummary(params, {}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    threw = false;
    try {
      PosteriorSummary p2(params, {{1.0, 2.0}});  // 2-d sample, 1-d params
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  return 0;
}
