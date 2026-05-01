#include "argus/line_list.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "argus/voigt.hpp"

namespace argus {

namespace {

// Reference temperature for HITRAN intensities and broadening.
constexpr double kTref_k    = 296.0;
constexpr double kPref_atm  = 1.0;
constexpr double kSpeedOfLightCm = 2.99792458e10;     // cm / s
constexpr double kBoltzmannSI    = 1.380649e-23;       // J / K
constexpr double kAmuKg          = 1.66053906660e-27;  // kg / amu
constexpr double kBarPerAtm      = 1.01325;            // bar in 1 atm

// Doppler standard deviation (cm^-1) for line at nu0_cm in K, mass amu.
double doppler_sigma_cm(double nu0_cm, double T_k, double mass_amu) {
  const double mass_kg = mass_amu * kAmuKg;
  const double sig_v = std::sqrt(kBoltzmannSI * T_k / mass_kg);  // m / s
  // sigma in cm^-1 = nu0 * sig_v / c (with c in cm / s)
  return nu0_cm * (sig_v * 100.0) / kSpeedOfLightCm;
}

// Pressure-broadening Lorentzian HWHM for one (T, P) state.
double lorentz_hwhm_cm(const Line& line, double T_k, double P_bar) {
  const double P_atm = P_bar / kBarPerAtm;
  const double tratio = std::pow(kTref_k / T_k, line.n_air);
  return tratio * line.gamma_air_cm * P_atm;
}

// Temperature dependence of line intensity (M2 keeps the simple
// E_lower-only form; M3 adds full partition-function ratios).
double intensity_T(const Line& line, double T_k) {
  // S(T) = S(Tref) * exp(-c2 * E_lower * (1/T - 1/Tref)) — partition
  // function ratio omitted at M2.
  constexpr double kC2 = 1.4387768775039338;  // h*c/k in cm * K
  const double exponent = -kC2 * line.E_lower_cm * (1.0 / T_k - 1.0 / kTref_k);
  return line.intensity_cm * std::exp(exponent);
}

}  // namespace

LineListOpacity::LineListOpacity(std::string key, std::vector<Line> lines,
                                 double molar_mass_amu)
    : key_(std::move(key)),
      lines_(std::move(lines)),
      molar_mass_amu_(molar_mass_amu) {
  if (molar_mass_amu_ <= 0.0) {
    throw std::invalid_argument(
        "LineListOpacity: molar_mass_amu must be positive");
  }
}

Tensor LineListOpacity::cross_section(
    const std::vector<double>& wavenumber_cm,
    const std::vector<double>& T_k,
    const std::vector<double>& P_bar) const {
  const std::size_t nT = T_k.size();
  const std::size_t nP = P_bar.size();
  const std::size_t nW = wavenumber_cm.size();

  Tensor out({nT, nP, nW});
  // Lines outside this many HWHMs of any wavenumber sample are skipped.
  // Conservative: 50 HWHMs ~ 1e-4 of the line peak, more than adequate
  // for M2 prototyping.
  constexpr double kCutoffHWHMs = 50.0;

  for (std::size_t iT = 0; iT < nT; ++iT) {
    const double T = T_k[iT];
    for (std::size_t iP = 0; iP < nP; ++iP) {
      const double P = P_bar[iP];
      // Pre-compute (sigma_g, gamma_l) for every line at this (T, P).
      for (const Line& line : lines_) {
        const double sigma_g = doppler_sigma_cm(line.nu0_cm, T,
                                                molar_mass_amu_);
        const double gamma_l = lorentz_hwhm_cm(line, T, P);
        const double S_T     = intensity_T(line, T);
        const double cutoff  = kCutoffHWHMs *
                               std::max(gamma_l,
                                        sigma_g * 2.354820045030949);

        for (std::size_t w = 0; w < nW; ++w) {
          const double dx = wavenumber_cm[w] - line.nu0_cm;
          if (std::fabs(dx) > cutoff) continue;
          const double phi = voigt(dx, sigma_g, gamma_l);
          // sigma(nu) [cm^2 / molecule] = S(T) [cm/molecule] * phi(nu) [cm]
          const std::size_t flat = (iT * nP + iP) * nW + w;
          out[flat] += S_T * phi;
        }
      }
    }
  }
  return out;
}

}  // namespace argus
