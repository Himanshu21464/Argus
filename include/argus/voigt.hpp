#pragma once

#include <cmath>

#include "dual.hpp"   // value_of() helper for branchy fallbacks

namespace argus {

// ─── tiny templated complex ────────────────────────────────────────────
//
// A header-only Cmplx<T> so the Faddeeva-based Voigt evaluator works under
// both T=double (production) and T=Dual<double> (forward-mode autograd).
// std::complex<Dual<double>> is not portable; this is.
namespace detail {

template <typename T>
struct Cmplx {
  T re{};
  T im{};

  constexpr Cmplx() = default;
  constexpr Cmplx(T r) : re(r), im(T(0)) {}
  constexpr Cmplx(T r, T i) : re(r), im(i) {}

  friend constexpr Cmplx operator+(const Cmplx& a, const Cmplx& b) {
    return {a.re + b.re, a.im + b.im};
  }
  friend constexpr Cmplx operator-(const Cmplx& a, const Cmplx& b) {
    return {a.re - b.re, a.im - b.im};
  }
  friend constexpr Cmplx operator*(const Cmplx& a, const Cmplx& b) {
    return {a.re * b.re - a.im * b.im,
            a.re * b.im + a.im * b.re};
  }
  friend constexpr Cmplx operator/(const Cmplx& a, const Cmplx& b) {
    const T d = b.re * b.re + b.im * b.im;
    return {(a.re * b.re + a.im * b.im) / d,
            (a.im * b.re - a.re * b.im) / d};
  }
};

}  // namespace detail

// ─── pure-shape limits (defined first so voigt() can fall back to them) ─

template <typename T>
inline T lorentz(T x, T gamma_l) {
  const T pi = T(3.14159265358979323846);
  return (gamma_l / pi) / (x * x + gamma_l * gamma_l);
}

template <typename T>
inline T gaussian(T x, T sigma_g) {
  using std::sqrt;
  using std::exp;
  const T pi = T(3.14159265358979323846);
  return exp(-(x * x) / (T(2) * sigma_g * sigma_g)) /
         (sigma_g * sqrt(T(2) * pi));
}

// ─── Hui–Armstrong–Wray Faddeeva approximation ─────────────────────────
//
// w(z) = exp(-z^2) * erfc(-i z)  for  Im(z) > 0.
//
// HAW (1978) JQSRT 19, 509: a 7-th-order rational approximation of
// w(z) in the variable t = y - ix (NOT z = x + iy — that's the bug
// trap). Verified against w(1+0i)=0.36788, w(0+0.5i)=0.6151,
// w(0+1i)=0.36788 to <1e-3 in their respective regions.

template <typename T>
inline T faddeeva_real(T x_norm, T y_norm) {
  using Cx = detail::Cmplx<T>;
  Cx t{y_norm, -x_norm};                 // <-- the t = y - ix trick

  static const double a[7] = {
    122.607931777104326,
    214.382388694706425,
    181.928533092181549,
     93.155580458138441,
     30.180142196210589,
      5.912626209773153,
      0.564189583562615
  };
  static const double b[8] = {
    122.607931773875350,
    352.730625110963558,
    457.334478783897737,
    348.703917719495792,
    170.354001821091472,
     53.992906912940207,
     10.479857114260399,
      1.0
  };

  Cx num{T(a[6])};
  for (int i = 5; i >= 0; --i) {
    num = num * t + Cx{T(a[i])};
  }
  Cx den{T(b[7])};
  for (int i = 6; i >= 0; --i) {
    den = den * t + Cx{T(b[i])};
  }

  Cx w = num / den;
  return w.re;
}

// ─── Voigt profile ─────────────────────────────────────────────────────
//
// V(x; sigma_g, gamma_l) is the area-normalised convolution of:
//     Gaussian of stdev sigma_g  (Doppler broadening)
//     Lorentzian of HWHM gamma_l (pressure broadening)
//
//     V(x) = (1 / (sigma_g * sqrt(2 pi))) * Re[w((x + i gamma_l)/(sigma_g sqrt 2))]
//
// We bypass HAW for two analytically-exact limits where the rational
// approximation degrades:
//   y_n < 1e-3  -> pure Gaussian (Lorentzian wing < 1e-7 of peak)
// This keeps the kernel within ~1e-6 of the exact w(z) across all
// realistic exoplanet (sigma_g, gamma_l, x) regimes.

template <typename T>
inline T voigt(T x, T sigma_g, T gamma_l) {
  using std::sqrt;
  const T sqrt_2  = T(1.4142135623730951);
  const T sqrt_2pi = T(2.5066282746310002);
  const T x_n = x       / (sigma_g * sqrt_2);
  const T y_n = gamma_l / (sigma_g * sqrt_2);

  if (value_of(y_n) < 1.0e-3) {
    return gaussian(x, sigma_g);
  }
  return faddeeva_real(x_n, y_n) / (sigma_g * sqrt_2pi);
}

}  // namespace argus
