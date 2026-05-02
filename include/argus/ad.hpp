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

}  // namespace argus::ad
