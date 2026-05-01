#pragma once

#include <cmath>
#include <iosfwd>

namespace argus {

// Forward-mode automatic differentiation via dual numbers.
//
// A Dual<T> is a pair (v, d) where v is the value and d is the derivative
// of v with respect to a single seed parameter. Arithmetic on Dual<T>
// propagates derivatives through the chain rule.
//
// To compute d f(x) / dx at x0, evaluate f with `Dual<T>{x0, T(1)}`. The
// returned Dual carries the derivative in `.d`.
//
// M2 ships scalar Dual<double>. M3 will add Dual<Vector> for batched
// gradients and a reverse-mode tape for many-parameter cases.

template <typename T>
struct Dual {
  T v{};
  T d{};

  constexpr Dual() = default;
  constexpr Dual(T value) : v(value), d(T(0)) {}
  constexpr Dual(T value, T deriv) : v(value), d(deriv) {}
};

// ─── arithmetic ────────────────────────────────────────────────────────

template <typename T>
constexpr Dual<T> operator+(const Dual<T>& a, const Dual<T>& b) {
  return {a.v + b.v, a.d + b.d};
}
template <typename T>
constexpr Dual<T> operator-(const Dual<T>& a, const Dual<T>& b) {
  return {a.v - b.v, a.d - b.d};
}
template <typename T>
constexpr Dual<T> operator*(const Dual<T>& a, const Dual<T>& b) {
  return {a.v * b.v, a.v * b.d + a.d * b.v};
}
template <typename T>
constexpr Dual<T> operator/(const Dual<T>& a, const Dual<T>& b) {
  return {a.v / b.v, (a.d * b.v - a.v * b.d) / (b.v * b.v)};
}
template <typename T>
constexpr Dual<T> operator-(const Dual<T>& a) { return {-a.v, -a.d}; }

// scalar / dual mixed
template <typename T> constexpr Dual<T> operator+(T s, const Dual<T>& a) { return Dual<T>{s} + a; }
template <typename T> constexpr Dual<T> operator+(const Dual<T>& a, T s) { return a + Dual<T>{s}; }
template <typename T> constexpr Dual<T> operator-(T s, const Dual<T>& a) { return Dual<T>{s} - a; }
template <typename T> constexpr Dual<T> operator-(const Dual<T>& a, T s) { return a - Dual<T>{s}; }
template <typename T> constexpr Dual<T> operator*(T s, const Dual<T>& a) { return Dual<T>{s} * a; }
template <typename T> constexpr Dual<T> operator*(const Dual<T>& a, T s) { return a * Dual<T>{s}; }
template <typename T> constexpr Dual<T> operator/(T s, const Dual<T>& a) { return Dual<T>{s} / a; }
template <typename T> constexpr Dual<T> operator/(const Dual<T>& a, T s) { return a / Dual<T>{s}; }

// in-place
template <typename T> Dual<T>& operator+=(Dual<T>& a, const Dual<T>& b) { a = a + b; return a; }
template <typename T> Dual<T>& operator-=(Dual<T>& a, const Dual<T>& b) { a = a - b; return a; }
template <typename T> Dual<T>& operator*=(Dual<T>& a, const Dual<T>& b) { a = a * b; return a; }
template <typename T> Dual<T>& operator/=(Dual<T>& a, const Dual<T>& b) { a = a / b; return a; }

// ─── std math overloads ────────────────────────────────────────────────

template <typename T>
inline Dual<T> exp(const Dual<T>& a) {
  const T e = std::exp(a.v);
  return {e, e * a.d};
}
template <typename T>
inline Dual<T> log(const Dual<T>& a) {
  return {std::log(a.v), a.d / a.v};
}
template <typename T>
inline Dual<T> sqrt(const Dual<T>& a) {
  const T s = std::sqrt(a.v);
  return {s, a.d / (T(2) * s)};
}
template <typename T>
inline Dual<T> pow(const Dual<T>& a, T p) {
  // d/dx [a^p] = p * a^(p-1) * a'
  const T pw = std::pow(a.v, p);
  const T dpw = p * std::pow(a.v, p - T(1)) * a.d;
  return {pw, dpw};
}
template <typename T>
inline Dual<T> sin(const Dual<T>& a) { return {std::sin(a.v),  std::cos(a.v) * a.d}; }
template <typename T>
inline Dual<T> cos(const Dual<T>& a) { return {std::cos(a.v), -std::sin(a.v) * a.d}; }
template <typename T>
inline T fabs(const Dual<T>& a) { return std::fabs(a.v); }  // value only

// ─── value-extraction helpers (for branchy code that needs to inspect
// a value while keeping derivatives composable) ───────────────────────

inline double value_of(double x) noexcept { return x; }
template <typename T>
inline double value_of(const Dual<T>& x) noexcept { return value_of(x.v); }

}  // namespace argus

// std::sqrt / std::exp / std::log inside templated callees pick up the
// argus::* overloads via ADL. These using-aliases inside namespace std are
// avoided intentionally; clients use ADL or call argus::voigt with a
// `using std::*; using argus::*;` block as needed.
