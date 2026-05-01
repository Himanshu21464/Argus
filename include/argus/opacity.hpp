#pragma once

#include <string>
#include <vector>

#include "atmosphere.hpp"
#include "tensor.hpp"

namespace argus {

// Source of molecular opacity. The production kernel will ship loaders for
// HITRAN, HITEMP, and ExoMol; this header is the stable surface they implement.
class OpacityKernel {
 public:
  virtual ~OpacityKernel() = default;

  virtual const std::string& species_key() const noexcept = 0;

  // Cross-section in cm^2 / molecule, broadcast over (T, P, wavenumber).
  // Inputs:
  //   wavenumber_cm  [n_wave]
  //   T_k            [n_T]
  //   P_bar          [n_P]
  // Output shape: [n_T, n_P, n_wave].
  virtual Tensor cross_section(const std::vector<double>& wavenumber_cm,
                               const std::vector<double>& T_k,
                               const std::vector<double>& P_bar) const = 0;

  // Cross-section accounting for self-broadening: the species' own VMR
  // contributes to pressure broadening via γ_self instead of γ_air.
  //   γ_total(layer) = γ_air * (1 - VMR_self) + γ_self * VMR_self
  // VMR_self_at_TP[i,j] gives the self-VMR for the (T_k[i], P_bar[j])
  // sample. Default implementation forwards to cross_section() (no
  // self-broadening), so existing kernels keep working.
  virtual Tensor cross_section_with_self(
      const std::vector<double>& wavenumber_cm,
      const std::vector<double>& T_k,
      const std::vector<double>& P_bar,
      const std::vector<double>& VMR_self_at_TP) const {
    (void)VMR_self_at_TP;
    return cross_section(wavenumber_cm, T_k, P_bar);
  }
};

// Placeholder Voigt-profile opacity used for tests and the M1 example.
// Real M3 implementation will replace this with a HITRAN line-list backed
// CUDA Voigt evaluator.
class GreyOpacity final : public OpacityKernel {
 public:
  GreyOpacity(std::string species_key, double sigma_cm2);

  const std::string& species_key() const noexcept override { return key_; }

  Tensor cross_section(const std::vector<double>& wavenumber_cm,
                       const std::vector<double>& T_k,
                       const std::vector<double>& P_bar) const override;

 private:
  std::string key_;
  double sigma_cm2_;
};

}  // namespace argus
