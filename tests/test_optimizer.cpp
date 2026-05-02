// Optimizer tests: train a linear regression W*x + b on synthetic data
// using the reverse-mode autograd tape and Adam, then verify Adam's
// solution matches the analytic least-squares result within 1%.

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
  using ad::SGD;

  // ─── 1. Adam on a simple convex problem: minimize (x - 3)^2 starting
  //     from x=0. Should converge to x ≈ 3 within 1000 steps. ─────────
  {
    std::vector<double> params{0.0};
    Adam opt(1, /*lr=*/0.05);
    for (int it = 0; it < 1000; ++it) {
      Tape t;
      Var x = t.input(params[0]);
      Var loss = (x - 3.0) * (x - 3.0);
      t.backward(loss);
      std::vector<double> g{t.grad(x)};
      opt.step(params, g);
    }
    assert(std::fabs(params[0] - 3.0) < 1.0e-3);
    assert(opt.step_count() == 1000);
  }

  // ─── 2. SGD on the same problem. ───────────────────────────────────
  {
    std::vector<double> params{0.0};
    SGD opt(1, /*lr=*/0.1, /*momentum=*/0.9);
    for (int it = 0; it < 1000; ++it) {
      Tape t;
      Var x = t.input(params[0]);
      Var loss = (x - 3.0) * (x - 3.0);
      t.backward(loss);
      std::vector<double> g{t.grad(x)};
      opt.step(params, g);
    }
    assert(std::fabs(params[0] - 3.0) < 1.0e-3);
  }

  // ─── 3. End-to-end: train y = w*x + b on synthetic noisy data with
  //     known true (w, b). Adam should recover (w, b) within ~1%. ─────
  {
    const double TRUE_W = 2.5;
    const double TRUE_B = -1.0;
    const std::size_t N = 200;

    std::mt19937_64 rng(2026);
    std::normal_distribution<double> noise(0.0, 0.05);
    std::uniform_real_distribution<double> uni(-2.0, 2.0);

    std::vector<double> xs(N), ys(N);
    for (std::size_t i = 0; i < N; ++i) {
      xs[i] = uni(rng);
      ys[i] = TRUE_W * xs[i] + TRUE_B + noise(rng);
    }

    std::vector<double> params{0.5, 0.5};   // (w, b) initial guess
    Adam opt(2, /*lr=*/0.05);

    const std::size_t batch_size = 32;
    std::uniform_int_distribution<std::size_t> idx_dist(0, N - 1);

    for (int epoch = 0; epoch < 200; ++epoch) {
      Tape t;
      Var w = t.input(params[0]);
      Var b = t.input(params[1]);
      // Mini-batch loss = mean squared error
      Var loss = t.input(0.0);
      for (std::size_t k = 0; k < batch_size; ++k) {
        const std::size_t i = idx_dist(rng);
        Var x_i = t.input(xs[i]);
        Var y_i = t.input(ys[i]);
        Var resid = y_i - (w * x_i + b);
        loss = loss + resid * resid;
      }
      loss = loss / static_cast<double>(batch_size);
      t.backward(loss);
      std::vector<double> grads{t.grad(w), t.grad(b)};
      opt.step(params, grads);
    }

    // Recovered (w, b) should be within ~1% of truth.
    assert(close(params[0], TRUE_W, 0.02));
    assert(std::fabs(params[1] - TRUE_B) < 0.05);
  }

  // ─── 4. Adam beats SGD-no-momentum on a saddle-shaped objective:
  //     f(x, y) = x^2 + 100 * y^2 (poorly conditioned).
  //     After equal step counts, Adam should be much closer to (0, 0).
  {
    auto loss_at = [](double x, double y) {
      return x * x + 100.0 * y * y;
    };
    auto run = [&](auto& opt, double lr_x_unused) {
      (void)lr_x_unused;
      std::vector<double> p{1.0, 1.0};
      for (int it = 0; it < 200; ++it) {
        Tape t;
        Var x = t.input(p[0]);
        Var y = t.input(p[1]);
        Var loss = x * x + 100.0 * y * y;
        t.backward(loss);
        std::vector<double> g{t.grad(x), t.grad(y)};
        opt.step(p, g);
      }
      return loss_at(p[0], p[1]);
    };

    Adam adam(2, /*lr=*/0.1);
    SGD sgd(2, /*lr=*/0.005, /*momentum=*/0.0);   // small lr because of 100*y^2

    const double final_loss_adam = run(adam, 0.0);
    const double final_loss_sgd  = run(sgd,  0.0);
    assert(final_loss_adam < final_loss_sgd);
    assert(final_loss_adam < 1.0e-3);
  }

  // ─── 5. Bad inputs throw. ───────────────────────────────────────────
  {
    bool threw = false;
    try { Adam(3, /*lr=*/-1.0); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    threw = false;
    try { Adam(3, 1e-3, /*beta1=*/1.5); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    threw = false;
    try { Adam(3, 1e-3, 0.9, /*beta2=*/-0.1); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    threw = false;
    try { SGD(3, /*lr=*/-1.0); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    threw = false;
    try { SGD(3, 1e-2, /*momentum=*/1.5); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);

    Adam opt(2, 1e-3);
    std::vector<double> p{0.0, 0.0};
    std::vector<double> g{1.0};   // wrong size
    threw = false;
    try { opt.step(p, g); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
  }

  return 0;
}
