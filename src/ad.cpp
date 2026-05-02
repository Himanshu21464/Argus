#include "argus/ad.hpp"

#include <cmath>
#include <stdexcept>

namespace argus::ad {

namespace {
[[noreturn]] void die(const char* msg) {
  throw std::runtime_error(std::string("argus::ad: ") + msg);
}
}  // namespace

// ─── Tape ─────────────────────────────────────────────────────────────

Var Tape::input(double value) {
  values_.push_back(value);
  parent_a_.push_back(kNoParent);
  parent_b_.push_back(kNoParent);
  local_grad_a_.push_back(0.0);
  local_grad_b_.push_back(0.0);
  return Var{this, values_.size() - 1, value};
}

void Tape::reset() noexcept {
  values_.clear();
  parent_a_.clear();
  parent_b_.clear();
  local_grad_a_.clear();
  local_grad_b_.clear();
  gradients_.clear();
}

Var Tape::record_unary(double value,
                       std::size_t parent,
                       double local_grad) {
  if (parent >= values_.size()) die("record_unary: parent out of range");
  values_.push_back(value);
  parent_a_.push_back(parent);
  parent_b_.push_back(kNoParent);
  local_grad_a_.push_back(local_grad);
  local_grad_b_.push_back(0.0);
  return Var{this, values_.size() - 1, value};
}

Var Tape::record_binary(double value,
                        std::size_t parent_a, std::size_t parent_b,
                        double local_grad_a, double local_grad_b) {
  if (parent_a >= values_.size() || parent_b >= values_.size()) {
    die("record_binary: parent out of range");
  }
  values_.push_back(value);
  parent_a_.push_back(parent_a);
  parent_b_.push_back(parent_b);
  local_grad_a_.push_back(local_grad_a);
  local_grad_b_.push_back(local_grad_b);
  return Var{this, values_.size() - 1, value};
}

void Tape::backward(const Var& output) {
  if (output.tape != this) die("backward: Var belongs to a different tape");
  if (output.idx >= values_.size()) die("backward: output index OOB");

  gradients_.assign(values_.size(), 0.0);
  gradients_[output.idx] = 1.0;

  // Reverse-topological traversal: parents are always recorded before
  // their children, so walking indices from `output.idx` downward
  // visits each node after all children that consume it.
  for (std::size_t i = output.idx + 1; i-- > 0;) {
    const double g = gradients_[i];
    if (g == 0.0) continue;
    const std::size_t pa = parent_a_[i];
    const std::size_t pb = parent_b_[i];
    if (pa != kNoParent) gradients_[pa] += g * local_grad_a_[i];
    if (pb != kNoParent) gradients_[pb] += g * local_grad_b_[i];
  }
}

double Tape::grad(const Var& v) const {
  if (v.tape != this) die("grad: Var belongs to a different tape");
  return grad(v.idx);
}

double Tape::grad(std::size_t idx) const {
  if (idx >= gradients_.size()) {
    die("grad: index out of range (did you call backward()?)");
  }
  return gradients_[idx];
}

// ─── Var operators ────────────────────────────────────────────────────

Var Var::operator+(const Var& o) const {
  return tape->record_binary(val + o.val, idx, o.idx, 1.0, 1.0);
}
Var Var::operator-(const Var& o) const {
  return tape->record_binary(val - o.val, idx, o.idx, 1.0, -1.0);
}
Var Var::operator*(const Var& o) const {
  return tape->record_binary(val * o.val, idx, o.idx, o.val, val);
}
Var Var::operator/(const Var& o) const {
  const double inv = 1.0 / o.val;
  return tape->record_binary(val * inv, idx, o.idx,
                             inv, -val * inv * inv);
}
Var Var::operator-() const {
  return tape->record_unary(-val, idx, -1.0);
}

Var Var::operator+(double s) const {
  return tape->record_unary(val + s, idx, 1.0);
}
Var Var::operator-(double s) const {
  return tape->record_unary(val - s, idx, 1.0);
}
Var Var::operator*(double s) const {
  return tape->record_unary(val * s, idx, s);
}
Var Var::operator/(double s) const {
  const double inv = 1.0 / s;
  return tape->record_unary(val * inv, idx, inv);
}

Var& Var::operator+=(const Var& o) { *this = *this + o; return *this; }
Var& Var::operator-=(const Var& o) { *this = *this - o; return *this; }
Var& Var::operator*=(const Var& o) { *this = *this * o; return *this; }
Var& Var::operator/=(const Var& o) { *this = *this / o; return *this; }

// ─── Free-standing math ───────────────────────────────────────────────

Var exp(const Var& a) {
  const double e = std::exp(a.val);
  return a.tape->record_unary(e, a.idx, e);
}
Var log(const Var& a) {
  return a.tape->record_unary(std::log(a.val), a.idx, 1.0 / a.val);
}
Var sqrt(const Var& a) {
  const double s = std::sqrt(a.val);
  return a.tape->record_unary(s, a.idx, 0.5 / s);
}
Var pow(const Var& a, double p) {
  const double v = std::pow(a.val, p);
  const double dv = p * std::pow(a.val, p - 1.0);
  return a.tape->record_unary(v, a.idx, dv);
}
Var sin(const Var& a) {
  return a.tape->record_unary(std::sin(a.val), a.idx, std::cos(a.val));
}
Var cos(const Var& a) {
  return a.tape->record_unary(std::cos(a.val), a.idx, -std::sin(a.val));
}
Var tanh(const Var& a) {
  const double t = std::tanh(a.val);
  return a.tape->record_unary(t, a.idx, 1.0 - t * t);
}

// scalar op Var (left scalar)
Var operator+(double s, const Var& a) { return a + s; }
Var operator-(double s, const Var& a) {
  return a.tape->record_unary(s - a.val, a.idx, -1.0);
}
Var operator*(double s, const Var& a) { return a * s; }
Var operator/(double s, const Var& a) {
  const double v = s / a.val;
  return a.tape->record_unary(v, a.idx, -s / (a.val * a.val));
}

}  // namespace argus::ad
