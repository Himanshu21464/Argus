#include "argus/line_list.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "argus/partition.hpp"
#include "argus/voigt.hpp"

namespace argus {

namespace {

// Reference temperature for HITRAN intensities and broadening.
constexpr double kTref_k    = 296.0;
constexpr double kSpeedOfLightCm = 2.99792458e10;     // cm / s
constexpr double kBoltzmannSI    = 1.380649e-23;       // J / K
constexpr double kAmuKg          = 1.66053906660e-27;  // kg / amu
constexpr double kBarPerAtm      = 1.01325;            // bar in 1 atm
constexpr double kC2             = 1.4387768775039338; // h*c/k in cm·K

// Doppler standard deviation (cm^-1) for line at nu0_cm in K, mass amu.
double doppler_sigma_cm(double nu0_cm, double T_k, double mass_amu) {
  const double mass_kg = mass_amu * kAmuKg;
  const double sig_v = std::sqrt(kBoltzmannSI * T_k / mass_kg);  // m / s
  return nu0_cm * (sig_v * 100.0) / kSpeedOfLightCm;
}

// Pressure-broadening Lorentzian HWHM for one (T, P) state.
double lorentz_hwhm_cm(const Line& line, double T_k, double P_bar) {
  const double P_atm = P_bar / kBarPerAtm;
  const double tratio = std::pow(kTref_k / T_k, line.n_air);
  return tratio * line.gamma_air_cm * P_atm;
}

// Full HITRAN-style temperature scaling of line intensity:
//   S(T) = S(Tref) * (Q(Tref) / Q(T))
//                  * exp(-c2 * E_lower * (1/T - 1/Tref))
//                  * (1 - exp(-c2 * nu0 / T))
//                  / (1 - exp(-c2 * nu0 / Tref))
// `q_ratio = Q(Tref)/Q(T)` is computed once per (T, species) by the caller.
double intensity_T(const Line& line, double T_k, double q_ratio) {
  const double boltz = std::exp(-kC2 * line.E_lower_cm *
                                (1.0 / T_k - 1.0 / kTref_k));
  const double induced = (1.0 - std::exp(-kC2 * line.nu0_cm / T_k)) /
                         (1.0 - std::exp(-kC2 * line.nu0_cm / kTref_k));
  return line.intensity_cm * q_ratio * boltz * induced;
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
    // Q-ratio is shared across all lines at this temperature when the
    // line list belongs to one species (the M2 case). If the species key
    // is not in the partition table, fall back to q_ratio = 1 — the
    // intensity scaling becomes Boltzmann-only.
    double q_ratio = 1.0;
    try {
      q_ratio = Partition::Q_ref(key_) / Partition::Q(key_, T);
    } catch (const std::invalid_argument&) {
      // unknown species; q_ratio stays 1
    }

    for (std::size_t iP = 0; iP < nP; ++iP) {
      const double P = P_bar[iP];
      for (const Line& line : lines_) {
        const double sigma_g = doppler_sigma_cm(line.nu0_cm, T,
                                                molar_mass_amu_);
        const double gamma_l = lorentz_hwhm_cm(line, T, P);
        const double S_T     = intensity_T(line, T, q_ratio);
        const double cutoff  = kCutoffHWHMs *
                               std::max(gamma_l,
                                        sigma_g * 2.354820045030949);

        for (std::size_t w = 0; w < nW; ++w) {
          const double dx = wavenumber_cm[w] - line.nu0_cm;
          if (std::fabs(dx) > cutoff) continue;
          const double phi = voigt(dx, sigma_g, gamma_l);
          const std::size_t flat = (iT * nP + iP) * nW + w;
          out[flat] += S_T * phi;
        }
      }
    }
  }
  return out;
}

}  // namespace argus
