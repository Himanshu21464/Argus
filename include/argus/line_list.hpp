#pragma once

#include <string>
#include <vector>

#include "opacity.hpp"
#include "tensor.hpp"

namespace argus {

// One molecular line — the M2 subset of the HITRAN .par columns we need
// to evaluate cross-section. M3+ will extend to the full HITRAN row.
//
// Units follow HITRAN convention:
//   nu0_cm           line centre wavenumber (cm^-1)
//   intensity_cm     line intensity at 296 K (cm^-1 / (molec cm^-2))
//   gamma_air_cm     air-broadened HWHM at 296 K, 1 atm (cm^-1)
//   gamma_self_cm    self-broadened HWHM at 296 K, 1 atm (cm^-1)
//   n_air            temperature exponent of gamma_air (dimensionless)
//   delta_air_cm     air pressure-induced line shift (cm^-1 / atm)
//   E_lower_cm       lower-state energy (cm^-1)
struct Line {
  double nu0_cm        = 0.0;
  double intensity_cm  = 0.0;
  double gamma_air_cm  = 0.05;
  double gamma_self_cm = 0.0;
  double n_air         = 0.5;
  double delta_air_cm  = 0.0;
  double E_lower_cm    = 0.0;
};

// LineListOpacity — sum a list of Voigt-shaped lines into a cross-section
// over a (T, P) grid.
//
// M2 ships a hand-rolled list (used by the example and tests). M3 will add
// HITRAN .par loaders and ExoMol bundle support behind the same interface.
class LineListOpacity final : public OpacityKernel {
 public:
  LineListOpacity(std::string species_key, std::vector<Line> lines,
                  double molar_mass_amu);

  const std::string& species_key() const noexcept override { return key_; }

  Tensor cross_section(const std::vector<double>& wavenumber_cm,
                       const std::vector<double>& T_k,
                       const std::vector<double>& P_bar) const override;

  Tensor cross_section_with_self(
      const std::vector<double>& wavenumber_cm,
      const std::vector<double>& T_k,
      const std::vector<double>& P_bar,
      const std::vector<double>& VMR_self_at_TP) const override;

  std::size_t num_lines() const noexcept { return lines_.size(); }

 private:
  std::string key_;
  std::vector<Line> lines_;
  double molar_mass_amu_;   // used by Doppler width
};

}  // namespace argus
