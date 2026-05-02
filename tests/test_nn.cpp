// Tests for the neural-net primitives (argus::nn::*).
// Building blocks for normalizing flows / amortized SBI in M3.5.

#include <cassert>
#include <cmath>
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

  // ─── 1. Linear forward with hand-set weights and bias ──────────────
  {
    Linear lin(3, 2);
    // y = W x + b
    //   W = [[1, 2, 3], [-1, 0, 4]]
    //   b = [10, -5]
    //   x = [1, 1, 1]
    //   y = [10 + 6, -5 + 3] = [16, -2]
    lin.set_weights({1.0, 2.0, 3.0, -1.0, 0.0, 4.0});
    lin.set_bias({10.0, -5.0});
    auto y = lin.forward({1.0, 1.0, 1.0});
    assert(y.size() == 2);
    assert(close(y[0], 16.0, 1.0e-12));
    assert(close(y[1], -2.0, 1.0e-12));
  }

  // ─── 2. Activations: piecewise checks ──────────────────────────────
  {
    std::vector<double> v{-2.0, -0.5, 0.0, 0.5, 2.0};

    auto vrelu = v;
    apply_activation(vrelu, Activation::ReLU);
    assert(vrelu[0] == 0.0 && vrelu[1] == 0.0);
    assert(vrelu[2] == 0.0 && vrelu[3] == 0.5 && vrelu[4] == 2.0);

    auto vleaky = v;
    apply_activation(vleaky, Activation::LeakyReLU);
    assert(close(vleaky[0], -0.02, 1.0e-12));
    assert(close(vleaky[1], -0.005, 1.0e-12));
    assert(vleaky[3] == 0.5);

    auto vt = v;
    apply_activation(vt, Activation::Tanh);
    for (std::size_t i = 0; i < v.size(); ++i) {
      assert(close(vt[i], std::tanh(v[i]), 1.0e-12));
    }

    auto vs = v;
    apply_activation(vs, Activation::Sigmoid);
    for (std::size_t i = 0; i < v.size(); ++i) {
      assert(close(vs[i], 1.0 / (1.0 + std::exp(-v[i])), 1.0e-12));
    }
  }

  // ─── 3. Sequential: a 2-2-2 MLP with identity weights gives the
  //     same output as the input after Tanh + final Linear. ─────────────
  {
    Sequential mlp(2, 2, {2}, Activation::Tanh);
    // Set hidden Linear to identity:
    mlp.layer(0).set_weights({1.0, 0.0, 0.0, 1.0});
    mlp.layer(0).set_bias({0.0, 0.0});
    // Set output Linear to identity:
    mlp.layer(1).set_weights({1.0, 0.0, 0.0, 1.0});
    mlp.layer(1).set_bias({0.0, 0.0});
    auto y = mlp.forward({0.5, -0.5});
    // After hidden: [0.5, -0.5] -> [tanh(0.5), tanh(-0.5)]
    // After output: same vector
    assert(close(y[0], std::tanh(0.5),  1.0e-12));
    assert(close(y[1], std::tanh(-0.5), 1.0e-12));
  }

  // ─── 4. Xavier init: deterministic given a seed; weights inside the
  //     expected uniform bound; biases zero. ────────────────────────────
  {
    Linear a(8, 4), b(8, 4);
    a.init_xavier(42);
    b.init_xavier(42);
    assert(a.weights().size() == 32);
    assert(a.bias().size() == 4);
    for (std::size_t i = 0; i < a.weights().size(); ++i) {
      assert(a.weights()[i] == b.weights()[i]);
    }
    const double bound = std::sqrt(6.0 / (8.0 + 4.0));
    for (double w : a.weights()) {
      assert(w >= -bound && w <= bound);
    }
    for (double bs : a.bias()) assert(bs == 0.0);
  }

  // ─── 5. Sequential with Xavier init runs deterministically and
  //     produces finite output. ────────────────────────────────────────
  {
    Sequential mlp(5, 3, {16, 16}, Activation::ReLU);
    mlp.init_xavier(2026);
    std::vector<double> x{0.1, -0.2, 0.3, -0.4, 0.5};
    auto y = mlp.forward(x);
    assert(y.size() == 3);
    for (double v : y) assert(std::isfinite(v));

    // Same seed -> same output.
    Sequential mlp2(5, 3, {16, 16}, Activation::ReLU);
    mlp2.init_xavier(2026);
    auto y2 = mlp2.forward(x);
    for (std::size_t i = 0; i < y.size(); ++i) assert(y[i] == y2[i]);

    // Different seed -> different output.
    Sequential mlp3(5, 3, {16, 16}, Activation::ReLU);
    mlp3.init_xavier(2027);
    auto y3 = mlp3.forward(x);
    bool any_diff = false;
    for (std::size_t i = 0; i < y.size(); ++i) {
      if (y[i] != y3[i]) { any_diff = true; break; }
    }
    assert(any_diff);
  }

  // ─── 6. Sequential layout: in/out dims and layer count. ────────────
  {
    Sequential mlp(4, 2, {8, 8, 8}, Activation::Tanh);
    assert(mlp.in_dim()   == 4);
    assert(mlp.out_dim()  == 2);
    assert(mlp.n_layers() == 4);   // 3 hidden + 1 output
    assert(mlp.layer(0).in_dim()  == 4);
    assert(mlp.layer(0).out_dim() == 8);
    assert(mlp.layer(3).in_dim()  == 8);
    assert(mlp.layer(3).out_dim() == 2);
  }

  // ─── 7. Bad inputs throw. ───────────────────────────────────────────
  {
    bool threw = false;
    try { Linear(0, 5); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    threw = false;
    try { Linear(5, 0); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    Linear lin(2, 2);
    threw = false;
    try { lin.set_weights({1.0, 2.0, 3.0}); }   // wrong size
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    threw = false;
    try { lin.set_bias({1.0}); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    threw = false;
    try { lin.forward({1.0}); }   // wrong input dim
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    threw = false;
    try { Sequential(2, 2, {0}, Activation::Tanh); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    Sequential mlp(2, 2, {4}, Activation::Tanh);
    threw = false;
    try { (void)mlp.layer(99); }
    catch (const std::out_of_range&) { threw = true; }
    assert(threw);
  }

  return 0;
}
