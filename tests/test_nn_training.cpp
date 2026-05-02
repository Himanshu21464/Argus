// End-to-end NN training: a single hidden-layer MLP fits a synthetic
// regression target using reverse-mode autograd + Adam. Verifies that
// the loss decreases substantially and final predictions track truth.
//
// This is the "amortized SBI" foundation: the same training pattern
// scales to normalizing flows once the user picks a richer architecture.

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
  using ad::Tape;
  using ad::Var;
  using ad::Adam;

  // ─── 1. ad::linear computes the same answer as nn::Linear::forward ─
  {
    nn::Linear lin(3, 2);
    lin.set_weights({1.0, 2.0, 3.0, -1.0, 0.0, 4.0});
    lin.set_bias({10.0, -5.0});
    std::vector<double> x{1.0, 1.0, 1.0};
    auto y_ref = lin.forward(x);

    Tape t;
    auto W = ad::to_vars(t, lin.weights());
    auto b = ad::to_vars(t, lin.bias());
    auto xv = ad::to_vars(t, x);
    auto y = ad::linear(t, W, b, xv, 3, 2);
    assert(y.size() == 2);
    assert(close(y[0].val, y_ref[0], 1.0e-12));
    assert(close(y[1].val, y_ref[1], 1.0e-12));
  }

  // ─── 2. ad::mse on equal pred and target → 0; on offset → variance.
  {
    Tape t;
    auto a = ad::to_vars(t, {1.0, 2.0, 3.0});
    auto b = ad::to_vars(t, {1.0, 2.0, 3.0});
    Var loss = ad::mse(a, b);
    assert(close(loss.val, 0.0, 1.0e-12));

    Tape t2;
    auto a2 = ad::to_vars(t2, {1.0, 2.0, 3.0});
    auto b2 = ad::to_vars(t2, {2.0, 3.0, 4.0});  // off by 1
    Var loss2 = ad::mse(a2, b2);
    assert(close(loss2.val, 1.0, 1.0e-12));
  }

  // ─── 3. End-to-end MLP training: fit f(x) = sin(2*x) on x ∈ [-1, 1].
  //     Use a 1-input, 16-hidden, 1-output MLP with tanh activation.
  //     Verify: final loss < 5% of initial loss; predictions track truth.
  {
    const std::size_t in_dim = 1;
    const std::size_t hid_dim = 16;
    const std::size_t out_dim = 1;

    // Initialise weights + biases (Xavier, single seed for determinism).
    nn::Linear l1(in_dim, hid_dim);
    nn::Linear l2(hid_dim, out_dim);
    l1.init_xavier(2026);
    l2.init_xavier(2027);

    // Mutable parameter copies that Adam updates.
    std::vector<double> W1 = l1.weights();
    std::vector<double> b1 = l1.bias();
    std::vector<double> W2 = l2.weights();
    std::vector<double> b2 = l2.bias();

    Adam adam(W1.size() + b1.size() + W2.size() + b2.size(), /*lr=*/0.01);

    // Synthetic data: 64 points from sin(2x) on [-1, 1].
    std::vector<double> xs, ys;
    for (int i = 0; i < 64; ++i) {
      const double x = -1.0 + (2.0 / 63.0) * static_cast<double>(i);
      xs.push_back(x);
      ys.push_back(std::sin(2.0 * x));
    }

    auto loss_at_params = [&](double& initial_loss_out) {
      Tape t;
      auto W1v = ad::to_vars(t, W1);
      auto b1v = ad::to_vars(t, b1);
      auto W2v = ad::to_vars(t, W2);
      auto b2v = ad::to_vars(t, b2);

      // Build per-sample predictions, accumulate MSE.
      Var loss = t.input(0.0);
      for (std::size_t i = 0; i < xs.size(); ++i) {
        auto x_i  = ad::to_vars(t, {xs[i]});
        auto h    = ad::linear(t, W1v, b1v, x_i, in_dim, hid_dim);
        auto h_a  = ad::tanh_vec(h);
        auto out  = ad::linear(t, W2v, b2v, h_a, hid_dim, out_dim);
        Var resid = out[0] - t.input(ys[i]);
        loss = loss + resid * resid;
      }
      loss = loss / static_cast<double>(xs.size());
      initial_loss_out = loss.val;
      t.backward(loss);

      auto gW1 = ad::grads_of(t, W1v);
      auto gb1 = ad::grads_of(t, b1v);
      auto gW2 = ad::grads_of(t, W2v);
      auto gb2 = ad::grads_of(t, b2v);

      // Concatenate parameters + gradients.
      std::vector<double> params, grads;
      params.insert(params.end(), W1.begin(), W1.end());
      params.insert(params.end(), b1.begin(), b1.end());
      params.insert(params.end(), W2.begin(), W2.end());
      params.insert(params.end(), b2.begin(), b2.end());
      grads.insert(grads.end(), gW1.begin(), gW1.end());
      grads.insert(grads.end(), gb1.begin(), gb1.end());
      grads.insert(grads.end(), gW2.begin(), gW2.end());
      grads.insert(grads.end(), gb2.begin(), gb2.end());

      adam.step(params, grads);

      // Splat back.
      std::size_t off = 0;
      for (std::size_t k = 0; k < W1.size(); ++k) W1[k] = params[off++];
      for (std::size_t k = 0; k < b1.size(); ++k) b1[k] = params[off++];
      for (std::size_t k = 0; k < W2.size(); ++k) W2[k] = params[off++];
      for (std::size_t k = 0; k < b2.size(); ++k) b2[k] = params[off++];
    };

    double init_loss = 0.0;
    loss_at_params(init_loss);
    double final_loss = init_loss;
    for (int epoch = 0; epoch < 800; ++epoch) {
      loss_at_params(final_loss);
    }
    // Final loss should be a small fraction of initial loss.
    assert(final_loss < 0.05 * init_loss);
    assert(final_loss < 0.02);   // absolute target

    // Verify the trained MLP predicts well at a held-out point.
    nn::Linear l1_trained(in_dim, hid_dim);
    nn::Linear l2_trained(hid_dim, out_dim);
    l1_trained.set_weights(W1);
    l1_trained.set_bias(b1);
    l2_trained.set_weights(W2);
    l2_trained.set_bias(b2);

    auto h = l1_trained.forward({0.42});
    for (auto& v : h) v = std::tanh(v);
    auto y = l2_trained.forward(h);
    const double truth = std::sin(2.0 * 0.42);
    assert(std::fabs(y[0] - truth) < 0.10);
  }

  // ─── 4. Bad inputs throw. ───────────────────────────────────────────
  {
    Tape t;
    auto W = ad::to_vars(t, {1.0, 2.0});
    auto b = ad::to_vars(t, {0.0});
    auto x = ad::to_vars(t, {1.0, 1.0});
    bool threw = false;
    try { (void)ad::linear(t, W, b, x, /*in_dim=*/3, /*out_dim=*/1); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);

    threw = false;
    try { (void)ad::mse({}, {}); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);

    threw = false;
    try {
      Tape t1;
      auto a = ad::to_vars(t1, {1.0, 2.0, 3.0});
      auto bvec = ad::to_vars(t1, {1.0, 2.0});
      (void)ad::mse(a, bvec);
    } catch (const std::runtime_error&) { threw = true; }
    assert(threw);
  }

  return 0;
}
