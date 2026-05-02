// Train a learnable scale+shift 1-D transform y = x · exp(s) + t via
// maximum-likelihood NLL against samples from N(2.5, 0.7²). Verify the
// trained (s, t) match the closed-form optimum:
//
//     s* = -log(σ), t* = -μ/σ
//
// This is the simplest 1-parameter "normalizing flow" — the same
// training loop (forward through transform on the tape → compute
// log-density → -log p as loss → backward → Adam) scales to the full
// `nn::NormalizingFlow` with stacked AffineCoupling layers in
// M3.5+.

#include <cassert>
#include <cmath>
#include <numeric>
#include <random>
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

  // Sample 500 points from N(2.5, 0.7²).
  const double TRUE_MU    = 2.5;
  const double TRUE_SIGMA = 0.7;
  std::mt19937_64 rng(2026);
  std::normal_distribution<double> nd(TRUE_MU, TRUE_SIGMA);
  std::vector<double> samples;
  samples.reserve(500);
  for (int i = 0; i < 500; ++i) samples.push_back(nd(rng));

  // Trainable scalars for the affine flow y = x·exp(s) + t.
  std::vector<double> params{0.0, 0.0};   // (s, t)
  Adam opt(2, /*lr=*/0.05);

  constexpr double kLog2Pi = 1.8378770664093453;

  // Mini-batch training loop.
  std::uniform_int_distribution<std::size_t> idx_dist(0, samples.size() - 1);
  const std::size_t batch_size = 64;

  double initial_loss = 0.0;
  double final_loss   = 0.0;

  for (int epoch = 0; epoch < 1500; ++epoch) {
    Tape t;
    Var sv = t.input(params[0]);
    Var tv = t.input(params[1]);

    Var nll = t.input(0.0);
    for (std::size_t b = 0; b < batch_size; ++b) {
      const double x_raw = samples[idx_dist(rng)];
      Var x  = t.input(x_raw);
      Var y  = x * ad::exp(sv) + tv;
      // log p_x(x) = log p_z(y) + log|det df/dx|
      //            = (-0.5 y² - 0.5 log(2π)) + s
      Var log_p = sv - 0.5 * (y * y) - 0.5 * kLog2Pi;
      nll = nll - log_p;
    }
    nll = nll / static_cast<double>(batch_size);

    if (epoch == 0) initial_loss = nll.val;
    final_loss = nll.val;

    t.backward(nll);
    std::vector<double> grads{t.grad(sv), t.grad(tv)};
    opt.step(params, grads);
  }

  const double S_OPT = -std::log(TRUE_SIGMA);
  const double T_OPT = -TRUE_MU / TRUE_SIGMA;

  // Loss should drop substantially.
  assert(final_loss < initial_loss - 0.5);

  // Recovered parameters should be within 5% of the analytic optimum.
  assert(close(params[0], S_OPT, 0.05, 0.05));
  assert(close(params[1], T_OPT, 0.05, 0.05));

  // Sanity check by transforming the data and verifying it's
  // approximately N(0, 1). Tolerance allows for finite-sample noise:
  // 500 draws from N(0, 1) have a sample-mean stddev of 1/√500 ≈ 0.045.
  double mean_y = 0.0, sq_y = 0.0;
  for (double x : samples) {
    const double y = x * std::exp(params[0]) + params[1];
    mean_y += y; sq_y += y * y;
  }
  mean_y /= static_cast<double>(samples.size());
  const double var_y = sq_y / static_cast<double>(samples.size()) - mean_y * mean_y;
  assert(std::fabs(mean_y) < 0.15);
  assert(close(var_y, 1.0, 0.15));

  return 0;
}
