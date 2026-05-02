#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace argus::ad {

// Reverse-mode automatic differentiation via a Wengert tape.
//
// Compared to forward-mode `Dual<T>` (which costs D forward passes
// for a D-parameter gradient), reverse-mode costs ONE forward pass
// plus one backward pass, regardless of D. This is the algorithm
// PyTorch/JAX/TensorFlow use; it's the only practical way to train
// neural networks (and normalizing flows) in C++.
//
// Usage:
//   Tape tape;
//   Var x = tape.input(2.0);
//   Var y = tape.input(3.0);
//   Var z = x * x + y * sin(x);
//   tape.backward(z);
//   double dz_dx = tape.grad(x);   // = 2x + y*cos(x)
//   double dz_dy = tape.grad(y);   // = sin(x)
//
// Vars are tiny handles (pointer + index + cached value); the tape
// owns all storage. The tape must outlive every Var derived from it.

class Tape;

struct Var {
  Tape* tape = nullptr;
  std::size_t idx = 0;
  double val = 0.0;

  Var() = default;
  Var(Tape* t, std::size_t i, double v) : tape(t), idx(i), val(v) {}

  // Arithmetic operators record nodes on the tape and return a new
  // Var pointing at the result.
  Var operator+(const Var& o) const;
  Var operator-(const Var& o) const;
  Var operator*(const Var& o) const;
  Var operator/(const Var& o) const;
  Var operator-() const;

  // Mixed Var / double overloads.
  Var operator+(double s) const;
  Var operator-(double s) const;
  Var operator*(double s) const;
  Var operator/(double s) const;

  // In-place
  Var& operator+=(const Var& o);
  Var& operator-=(const Var& o);
  Var& operator*=(const Var& o);
  Var& operator/=(const Var& o);
};

// Free-form math overloads (call into Tape).
Var exp (const Var& a);
Var log (const Var& a);
Var sqrt(const Var& a);
Var pow (const Var& a, double p);
Var sin (const Var& a);
Var cos (const Var& a);
Var tanh(const Var& a);

// Mixed double op Var (left scalar).
Var operator+(double s, const Var& a);
Var operator-(double s, const Var& a);
Var operator*(double s, const Var& a);
Var operator/(double s, const Var& a);

class Tape {
 public:
  Tape() = default;
  Tape(const Tape&) = delete;
  Tape& operator=(const Tape&) = delete;

  // Create a leaf node (no parents). Returns a Var holding its index.
  Var input(double value);

  // Reset the tape — invalidates all existing Vars.
  void reset() noexcept;

  std::size_t size() const noexcept { return values_.size(); }

  // ─── recording API used by Var operators (kept public so free-
  //     standing math functions can call them) ────────────────────────
  Var record_unary(double value,
                   std::size_t parent,
                   double local_grad);
  Var record_binary(double value,
                    std::size_t parent_a, std::size_t parent_b,
                    double local_grad_a, double local_grad_b);

  // Compute gradients of `output` with respect to every leaf on the
  // tape. After this call, `grad(var)` returns d(output)/d(var).
  void backward(const Var& output);

  // Gradient lookup by Var or by index.
  double grad(const Var& v) const;
  double grad(std::size_t idx) const;

 private:
  // Sentinel for "no parent" on leaf nodes.
  static constexpr std::size_t kNoParent =
      static_cast<std::size_t>(-1);

  std::vector<double>      values_;
  std::vector<std::size_t> parent_a_;
  std::vector<std::size_t> parent_b_;
  std::vector<double>      local_grad_a_;
  std::vector<double>      local_grad_b_;
  std::vector<double>      gradients_;       // populated by backward()
};

// ─── Optimizers ──────────────────────────────────────────────────────
//
// These operate on raw `double` parameter arrays + their gradients.
// The training loop pattern is:
//
//   std::vector<double> params = ...;
//   Adam opt(params.size(), 1e-2);
//   for (epoch = ...) {
//     Tape t;
//     std::vector<Var> ps;
//     for (double p : params) ps.push_back(t.input(p));
//     Var loss = compute_loss(t, ps, batch);
//     t.backward(loss);
//     std::vector<double> grads(params.size());
//     for (i = ...) grads[i] = t.grad(ps[i]);
//     opt.step(params, grads);
//   }

class Adam {
 public:
  Adam(std::size_t n_params,
       double lr = 1.0e-3,
       double beta1 = 0.9,
       double beta2 = 0.999,
       double eps = 1.0e-8);

  // Apply one optimizer step: params -= lr * (m_hat / (sqrt(v_hat) + eps))
  // m and v are the bias-corrected first/second-moment estimates.
  void step(std::vector<double>& params,
            const std::vector<double>& grads);

  std::size_t step_count() const noexcept { return t_; }

 private:
  std::vector<double> m_;       // 1st moment
  std::vector<double> v_;       // 2nd moment
  double lr_, beta1_, beta2_, eps_;
  std::size_t t_ = 0;
};

// SGD with momentum. Useful baseline; Adam usually outperforms.
class SGD {
 public:
  SGD(std::size_t n_params, double lr = 1.0e-2, double momentum = 0.0);

  void step(std::vector<double>& params,
            const std::vector<double>& grads);

 private:
  std::vector<double> velocity_;
  double lr_, momentum_;
};

// ─── NN forward helpers on the autograd tape ──────────────────────────
//
// These allow trainable neural-network layers to be built from
// `argus::ad::Var` inputs so backward() through the loss yields
// gradients on the weights and biases. The pattern is:
//
//   Tape t;
//   auto W = ad::to_vars(t, weight_data);
//   auto b = ad::to_vars(t, bias_data);
//   auto x = ad::to_vars(t, input);
//   auto h = ad::linear(t, W, b, x, in_dim, out_dim);
//   auto y_hat = ad::tanh_vec(h);
//   Var loss = mse(y_hat, target);
//   t.backward(loss);
//   auto W_grads = ad::grads_of(t, W);  // -> apply Adam.step
//
// These free functions compose with anything that produces `Var`s.

// Convert a vector of doubles into leaf Vars on a tape.
std::vector<Var> to_vars(Tape& t, const std::vector<double>& xs);

// Read gradients of every Var in `vs` after backward() has been called.
std::vector<double> grads_of(const Tape& t, const std::vector<Var>& vs);

// Linear layer forward: y = W * x + b.
//   weights row-major [out_dim x in_dim]
//   bias    [out_dim]
std::vector<Var> linear(Tape& t,
                        const std::vector<Var>& weights,
                        const std::vector<Var>& bias,
                        const std::vector<Var>& input,
                        std::size_t in_dim, std::size_t out_dim);

// Element-wise activations on a vector<Var>.
std::vector<Var> relu_vec(const std::vector<Var>& xs);
std::vector<Var> tanh_vec(const std::vector<Var>& xs);
std::vector<Var> sigmoid_vec(const std::vector<Var>& xs);

// Mean squared error: (1/N) * Σ (pred - target)^2 as a single Var.
Var mse(const std::vector<Var>& pred, const std::vector<Var>& target);

}  // namespace argus::ad
