// Tests for the M3.5 ConditionalNormalizingFlow — the architecture
// that backs amortized simulation-based inference (NPE, SNPE, DINGO).
//
// Verifies the algebra is correct (forward/inverse, log-density,
// determinism) and that the flow's density really depends on the
// conditioning vector. Training (gradient through the flow stack)
// is a separate test once the AD-aware coupling layers land in M4.5.

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
  using namespace argus::nn;

  // ─── 1. Bad-input construction throws. ─────────────────────────────
  {
    bool threw;
    threw = false;
    try { ConditionalAffineCoupling(/*dim=*/1, /*cond_dim=*/2,
                                    /*split=*/1, {16}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { ConditionalAffineCoupling(/*dim=*/4, /*cond_dim=*/2,
                                    /*split=*/0, {16}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { ConditionalAffineCoupling(/*dim=*/4, /*cond_dim=*/2,
                                    /*split=*/4, {16}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { ConditionalNormalizingFlow(/*dim=*/4, /*cond_dim=*/2,
                                     /*n_couplings=*/0, /*split=*/2, {16}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { ConditionalNormalizingFlow(/*dim=*/1, /*cond_dim=*/2,
                                     /*n_couplings=*/4, /*split=*/1, {16}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // ─── 2. Wrong-shape forward / inverse throws. ──────────────────────
  {
    ConditionalAffineCoupling c(4, 3, 2, {16});
    c.init_xavier(7);
    bool threw;
    threw = false;
    try { (void)c.forward(std::vector<double>{1, 2, 3}, {0.0, 0.0, 0.0}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { (void)c.forward({1, 2, 3, 4}, {0.0, 0.0}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // ─── 3. ConditionalAffineCoupling: forward → inverse round-trip. ───
  //     Real NVP coupling is exactly invertible — round-trip must give
  //     bit-equality up to FP noise. ───────────────────────────────────
  {
    ConditionalAffineCoupling c(/*dim=*/4, /*cond_dim=*/3,
                                /*split=*/2, {16, 16});
    c.init_xavier(/*seed=*/42);
    std::mt19937_64 rng(2026);
    std::normal_distribution<double> nd(0.0, 1.0);
    for (int trial = 0; trial < 8; ++trial) {
      std::vector<double> x(4), cond(3);
      for (auto& v : x)    v = nd(rng);
      for (auto& v : cond) v = nd(rng);
      auto fw = c.forward(x, cond);
      auto bw = c.inverse(fw.y, cond);
      for (std::size_t i = 0; i < x.size(); ++i) {
        assert(close(bw.y[i], x[i], 1.0e-12, 1.0e-12));
      }
      // Inverse log-det is the negative of the forward log-det.
      assert(close(bw.log_det_jacobian, -fw.log_det_jacobian, 1.0e-12, 1.0e-12));
    }
  }

  // ─── 4. ConditionalNormalizingFlow: forward → inverse round-trip. ──
  {
    ConditionalNormalizingFlow flow(/*dim=*/4, /*cond_dim=*/3,
                                    /*n_couplings=*/6, /*split=*/2,
                                    {16, 16});
    flow.init_xavier(/*seed=*/123);
    std::mt19937_64 rng(2026);
    std::normal_distribution<double> nd(0.0, 1.0);
    for (int trial = 0; trial < 8; ++trial) {
      std::vector<double> x(4), cond(3);
      for (auto& v : x)    v = nd(rng);
      for (auto& v : cond) v = nd(rng);
      auto fw = flow.forward(x, cond);
      auto bw = flow.inverse(fw.y, cond);
      for (std::size_t i = 0; i < x.size(); ++i) {
        assert(close(bw.y[i], x[i], 1.0e-10, 1.0e-10));
      }
      // Cumulative log-det signs cancel.
      assert(close(bw.log_det_jacobian, -fw.log_det_jacobian, 1.0e-10, 1.0e-10));
    }
  }

  // ─── 5. log_density consistent with the change-of-variables formula.
  //     log p_x(x | y) = log N(z; 0, I) + log|det df/dx|
  //                    = -0.5 |z|² - 0.5 D log(2π) + Σ log_det
  //     Recompute the Gaussian piece from forward.y and check. ────────
  {
    ConditionalNormalizingFlow flow(4, 3, 6, 2, {16, 16});
    flow.init_xavier(11);
    std::mt19937_64 rng(2026);
    std::normal_distribution<double> nd(0.0, 1.0);
    constexpr double kLog2Pi = 1.8378770664093453;
    for (int trial = 0; trial < 5; ++trial) {
      std::vector<double> x(4), cond(3);
      for (auto& v : x)    v = nd(rng);
      for (auto& v : cond) v = nd(rng);
      auto fw = flow.forward(x, cond);
      double zz = 0.0;
      for (double v : fw.y) zz += v * v;
      const double expected =
          -0.5 * zz - 0.5 * 4.0 * kLog2Pi + fw.log_det_jacobian;
      assert(close(flow.log_density(x, cond), expected, 1.0e-12, 1.0e-12));
    }
  }

  // ─── 6. Conditional dependence: changing the conditioning vector
  //     changes both the forward output and the log-density. ─────────
  {
    ConditionalNormalizingFlow flow(4, 3, 6, 2, {16, 16});
    flow.init_xavier(99);
    std::vector<double> x{0.3, -0.1, 0.7, -0.5};
    std::vector<double> cond_a{0.0, 0.0, 0.0};
    std::vector<double> cond_b{1.0, -2.0, 0.5};
    auto fa = flow.forward(x, cond_a);
    auto fb = flow.forward(x, cond_b);
    bool any_diff = false;
    for (std::size_t i = 0; i < x.size(); ++i) {
      if (std::fabs(fa.y[i] - fb.y[i]) > 1.0e-6) { any_diff = true; break; }
    }
    assert(any_diff);
    assert(std::fabs(flow.log_density(x, cond_a) -
                     flow.log_density(x, cond_b)) > 1.0e-6);
  }

  // ─── 7. init_xavier reproducibility: same seed → bit-equal output.
  //     Two flows initialised with the same seed must produce
  //     identical forward outputs on identical (x, cond) inputs. ────
  {
    ConditionalNormalizingFlow flow_a(4, 3, 6, 2, {16, 16});
    ConditionalNormalizingFlow flow_b(4, 3, 6, 2, {16, 16});
    flow_a.init_xavier(2026);
    flow_b.init_xavier(2026);
    std::vector<double> x{0.3, 0.1, -0.7, 0.5};
    std::vector<double> cond{0.2, -0.4, 1.0};
    auto a = flow_a.forward(x, cond);
    auto b = flow_b.forward(x, cond);
    for (std::size_t i = 0; i < x.size(); ++i) {
      assert(a.y[i] == b.y[i]);
    }
    assert(a.log_det_jacobian == b.log_det_jacobian);
  }

  // ─── 8. Sampling: inverse(z, cond) returns a vector of size dim;
  //     determinism: identical RNG state → identical sample. ────────
  {
    ConditionalNormalizingFlow flow(4, 3, 6, 2, {16, 16});
    flow.init_xavier(7);
    std::vector<double> cond{0.0, 0.0, 0.0};
    std::mt19937_64 r1(42), r2(42);
    auto s1 = flow.sample(cond, r1);
    auto s2 = flow.sample(cond, r2);
    assert(s1.size() == 4);
    for (std::size_t i = 0; i < s1.size(); ++i) assert(s1[i] == s2[i]);
  }

  // ─── 9. Different conditioning → different sample (with same RNG).
  //     This is the core amortized-SBI property: the flow learns a
  //     family of distributions parameterised by the observation. ───
  {
    ConditionalNormalizingFlow flow(4, 3, 6, 2, {16, 16});
    flow.init_xavier(7);
    std::vector<double> cond_a{0.0, 0.0, 0.0};
    std::vector<double> cond_b{2.0, -1.0, 0.5};
    std::mt19937_64 r1(42), r2(42);
    auto s_a = flow.sample(cond_a, r1);
    auto s_b = flow.sample(cond_b, r2);
    bool any_diff = false;
    for (std::size_t i = 0; i < s_a.size(); ++i) {
      if (std::fabs(s_a[i] - s_b[i]) > 1.0e-6) { any_diff = true; break; }
    }
    assert(any_diff);
  }

  return 0;
}
