#pragma once

#include <cmath>

namespace argus {

// Voigt line-shape evaluation.
//
// V(x; sigma_g, gamma_l) is the convolution of:
//     a Gaussian of standard deviation sigma_g (Doppler broadening)
//     a Lorentzian of half-width gamma_l       (pressure broadening)
//
// We use the Thompson-Cox-Hastings (1987) pseudo-Voigt approximation:
// a weighted sum of a Gaussian and a Lorentzian whose width matches a
// closed-form expression for the Voigt FWHM. The approximation is good to
// ~1% across the typical Doppler/Lorentz ratio range encountered in
// exoplanet atmospheres, and is cheap enough to put in the inner loop.
//
// All arguments are area-normalised so that integral V dx = 1.
//
// Templated so the same code path can be evaluated under
// argus::Dual<double> for forward-mode autograd.
//
// Inputs:
//   x        offset from the line centre (cm^-1)
//   sigma_g  Gaussian standard deviation (cm^-1) — Doppler width / sqrt(2)
//   gamma_l  Lorentzian half-width at half maximum (cm^-1) — pressure width

template <typename T>
inline T voigt(T x, T sigma_g, T gamma_l) {
  using std::sqrt;
  using std::log;
  using std::exp;

  // Gaussian and Lorentzian FWHM
  const T two_sqrt_2ln2 = T(2.354820045030949);   // 2*sqrt(2*ln 2)
  const T pi = T(3.14159265358979323846);
  const T ln2 = T(0.6931471805599453);

  const T fG = sigma_g * two_sqrt_2ln2;          // Gaussian FWHM
  const T fL = gamma_l * T(2);                   // Lorentzian FWHM

  // Olivero-Longbothum FWHM combination
  const T fG2 = fG * fG;
  const T fG3 = fG2 * fG;
  const T fG4 = fG2 * fG2;
  const T fG5 = fG4 * fG;
  const T fL2 = fL * fL;
  const T fL3 = fL2 * fL;
  const T fL4 = fL2 * fL2;
  const T fL5 = fL4 * fL;

  const T fV5 = fG5
              + T(2.69269) * fG4 * fL
              + T(2.42843) * fG3 * fL2
              + T(4.47163) * fG2 * fL3
              + T(0.07842) * fG  * fL4
              + fL5;
  // fifth root via exp(log(.)/5) — avoids needing a pow(Dual,Dual) overload
  // when T is Dual<double>. ADL picks the right log/exp for both branches.
  const T fV = exp(log(fV5) * T(0.2));           // Voigt FWHM

  const T r = fL / fV;
  const T eta = T(1.36603) * r
              - T(0.47719) * r * r
              + T(0.11116) * r * r * r;          // Lorentzian fraction

  // Build a Gaussian and a Lorentzian both with FWHM fV, area-normalised,
  // then mix.
  const T sig = fV / two_sqrt_2ln2;              // Gaussian std-dev with FWHM=fV
  const T g_ampl = sqrt(ln2 / pi) / (T(0.5) * fV); // peak of normalised Gaussian
  // simpler: f_G(x) = sqrt(ln2/pi) * (2/fV) * exp(-ln2 * (2x/fV)^2)
  const T u = T(2) * x / fV;
  const T fG_x = (T(2) / fV) * sqrt(ln2 / pi) * exp(-ln2 * u * u);

  const T half_fV = T(0.5) * fV;
  const T fL_x = (half_fV / pi) / (x * x + half_fV * half_fV);

  (void)sig; (void)g_ampl;  // kept for clarity; not in final formula
  return eta * fL_x + (T(1) - eta) * fG_x;
}

// Pure-Lorentz fallback (diagnostic / sanity).
template <typename T>
inline T lorentz(T x, T gamma_l) {
  const T pi = T(3.14159265358979323846);
  return (gamma_l / pi) / (x * x + gamma_l * gamma_l);
}

// Pure-Gaussian fallback (diagnostic / sanity).
template <typename T>
inline T gaussian(T x, T sigma_g) {
  using std::sqrt;
  using std::exp;
  const T pi = T(3.14159265358979323846);
  return exp(-(x * x) / (T(2) * sigma_g * sigma_g)) /
         (sigma_g * sqrt(T(2) * pi));
}

}  // namespace argus
