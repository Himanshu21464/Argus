// Tests for the NormalizingFlow primitive.
//
// At Xavier init the flow is ≈ identity (small log_det), so:
//   * forward(x) returns z near x with small log_det
//   * inverse(forward(x)) should be bit-exact x
//   * log_density(x) is finite and equals what change-of-variables predicts

#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "argus/argus.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

}  // namespace

int main() {
  using namespace argus::nn;

  // ─── 1. Inverse-of-forward recovers input bit-exact. ────────────────
  {
    NormalizingFlow flow(/*dim=*/4, /*n_couplings=*/6, /*split=*/2,
                         /*hidden_dims=*/{16, 16}, Activation::Tanh);
    flow.init_xavier(2026);
    std::vector<double> x{0.3, -0.7, 1.1, -0.5};
    auto fw = flow.forward(x);
    auto bw = flow.inverse(fw.y);
    for (std::size_t i = 0; i < x.size(); ++i) {
      assert(std::fabs(bw.y[i] - x[i]) < 1.0e-10);
    }
    // Inverse log-det should equal -forward log-det.
    assert(std::fabs(bw.log_det_jacobian + fw.log_det_jacobian) < 1.0e-10);
  }

  // ─── 2. log_density is finite and matches the change-of-variables. ──
  {
    NormalizingFlow flow(3, 4, 1, {12}, Activation::Tanh);
    flow.init_xavier(7);
    std::vector<double> x{0.2, -0.4, 0.7};
    const double lp = flow.log_density(x);
    assert(std::isfinite(lp));
    // Manually compute via forward + base Gaussian.
    auto fw = flow.forward(x);
    double zz = 0.0;
    for (double v : fw.y) zz += v * v;
    constexpr double kLog2Pi = 1.8378770664093453;
    const double expected =
        -0.5 * zz - 0.5 * 3.0 * kLog2Pi + fw.log_det_jacobian;
    assert(close(lp, expected, 1.0e-12));
  }

  // ─── 3. sample(rng) is reproducible given the same seed and gives a
  //     valid (finite) sample. ───────────────────────────────────────────
  {
    NormalizingFlow flow(3, 4, 1, {12}, Activation::Tanh);
    flow.init_xavier(7);
    std::mt19937_64 a(99), b(99);
    auto sa = flow.sample(a);
    auto sb = flow.sample(b);
    assert(sa.size() == 3);
    assert(sb.size() == 3);
    for (std::size_t i = 0; i < 3; ++i) {
      assert(sa[i] == sb[i]);                    // bit-exact same RNG
      assert(std::isfinite(sa[i]));
    }
  }

  // ─── 4. With one coupling and zeroed conditioner the flow IS identity
  //     (no intermediate swaps). forward(x) = x, log_det = 0,
  //     log_density(x) is the standard Gaussian density. ─────────────────
  {
    NormalizingFlow flow(2, /*n_couplings=*/1, /*split=*/1,
                         {4}, Activation::Tanh);
    // Zero all weights and biases on the single coupling's conditioner.
    auto& cond = flow.coupling(0).conditioner();
    for (std::size_t l = 0; l < cond.n_layers(); ++l) {
      const std::size_t W = cond.layer(l).in_dim() *
                            cond.layer(l).out_dim();
      cond.layer(l).set_weights(std::vector<double>(W, 0.0));
      cond.layer(l).set_bias(
          std::vector<double>(cond.layer(l).out_dim(), 0.0));
    }
    std::vector<double> x{0.4, 0.7};
    auto fw = flow.forward(x);
    for (std::size_t i = 0; i < x.size(); ++i) {
      assert(std::fabs(fw.y[i] - x[i]) < 1.0e-12);
    }
    assert(std::fabs(fw.log_det_jacobian) < 1.0e-12);
    // log_density should match the standard 2-D Gaussian at x.
    double zz = x[0] * x[0] + x[1] * x[1];
    const double expected =
        -0.5 * zz - 0.5 * 2.0 * 1.8378770664093453;
    assert(close(flow.log_density(x), expected, 1.0e-12));
  }

  // ─── 5. Save / load round-trip: weights bit-exact via hex-float. ────
  {
    NormalizingFlow flow_a(4, 3, 2, {16, 16}, Activation::Tanh);
    flow_a.init_xavier(2026);
    const std::string path = "/tmp/argus_flow_test.txt";
    flow_a.save(path);

    NormalizingFlow flow_b(4, 3, 2, {16, 16}, Activation::Tanh);
    // Different seed gives different starting weights; load should overwrite.
    flow_b.init_xavier(99);
    flow_b.load(path);

    // Forward output should agree bit-exact for both flows.
    std::vector<double> x{0.3, -0.4, 0.5, -0.6};
    auto fa = flow_a.forward(x);
    auto fb = flow_b.forward(x);
    for (std::size_t i = 0; i < fa.y.size(); ++i) {
      assert(fa.y[i] == fb.y[i]);
    }
    assert(fa.log_det_jacobian == fb.log_det_jacobian);
    std::remove(path.c_str());
  }

  // ─── 6. Load with shape mismatch throws. ─────────────────────────────
  {
    NormalizingFlow flow_a(4, 2, 2, {16}, Activation::Tanh);
    flow_a.init_xavier(7);
    const std::string path = "/tmp/argus_flow_mismatch.txt";
    flow_a.save(path);

    // Try to load into a flow with different architecture.
    NormalizingFlow flow_b(4, 2, 2, {32}, Activation::Tanh);  // wider hidden
    bool threw = false;
    try { flow_b.load(path); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    std::remove(path.c_str());
  }

  // ─── 7. Missing file throws on load. ─────────────────────────────────
  {
    NormalizingFlow flow(4, 2, 2, {8});
    bool threw = false;
    try { flow.load("/tmp/argus_nonexistent_flow_12345.txt"); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
  }

  // ─── 7b. Truncated file: file has fewer layers than the flow expects.
  //     Earlier code silently left the remaining couplings in their
  //     pre-load state — now caught at end-of-stream. ──────────────
  {
    NormalizingFlow flow_full(4, /*n_couplings=*/2, 2, {8});
    flow_full.init_xavier(123);
    const std::string path = "/tmp/argus_flow_truncated.txt";
    flow_full.save(path);

    // Truncate after the first coupling's first layer's bias block.
    // Read full file, find a safe truncation point, write back.
    std::ifstream src(path);
    std::string line, kept;
    int layer_count = 0;
    while (std::getline(src, line)) {
      if (line.rfind("layer ", 0) == 0) {
        ++layer_count;
        if (layer_count > 1) break;        // stop before second layer
      }
      kept += line + "\n";
    }
    src.close();
    std::ofstream dst(path, std::ios::trunc);
    dst << kept;
    dst.close();

    NormalizingFlow flow_target(4, 2, 2, {8});
    bool threw = false;
    try { flow_target.load(path); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    std::remove(path.c_str());
  }

  // ─── 8. Bad inputs throw. ───────────────────────────────────────────
  {
    bool threw = false;
    try { NormalizingFlow(1, 4, 0, {8}); }      // dim < 2
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { NormalizingFlow(4, 0, 1, {8}); }      // n_couplings = 0
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { NormalizingFlow(4, 1, 0, {8}); }      // split = 0
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { NormalizingFlow(4, 1, 4, {8}); }      // split >= dim
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    NormalizingFlow flow(4, 2, 2, {8});
    threw = false;
    try { flow.forward({1.0, 2.0, 3.0}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { (void)flow.coupling(99); }
    catch (const std::out_of_range&) { threw = true; }
    assert(threw);
  }

  return 0;
}
