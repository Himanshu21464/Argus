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

// Rayleigh-scattering opacity. Cross-section scales as λ⁻⁴ — strongly
// peaked toward short wavelengths. The reference value at λ = 1 μm is
// configurable; for an H2-dominated background (typical hot Jupiter)
// use σ_1μm ≈ 8.49e-29 cm²; for an H2/He mix ~6e-29 cm². At ν in cm⁻¹:
//
//     σ(ν) = σ_1μm · (ν / 10000)^4         where 10000 cm⁻¹ ≡ 1 μm
//
// This kernel doesn't care about (T, P) — Rayleigh is a continuum
// scatterer. The TransmissionModel multiplies σ by VMR · n_total in
// the standard way, so to apply Rayleigh to the bulk gas just set the
// associated species' VMR to 1.0 in every layer.
class RayleighOpacity final : public OpacityKernel {
 public:
  RayleighOpacity(std::string species_key, double sigma_at_1um_cm2);

  const std::string& species_key() const noexcept override { return key_; }

  Tensor cross_section(const std::vector<double>& wavenumber_cm,
                       const std::vector<double>& T_k,
                       const std::vector<double>& P_bar) const override;

  double sigma_at_1um_cm2() const noexcept { return sigma_1um_; }

 private:
  std::string key_;
  double sigma_1um_;
};

// Canonical "gray cloud deck" model used in exoplanet retrievals.
// Above the cloud-top pressure P_cloud_bar the atmosphere is fully
// opaque (huge sigma per cloud molecule), below it the cloud
// contributes nothing. The species_key is a canonical bookkeeping
// label ("CLOUD_DECK") and the kernel attaches its opacity to the
// fictitious species's VMR (typically set to 1.0 for layers below
// the cloud and 0 above).
//
// The expected wiring: caller adds a "CLOUD_DECK" species to the
// Atmosphere with VMR = 1.0 in every layer, then sets up the
// CloudDeckOpacity. For each (T,P) sample, this kernel returns:
//   sigma_cm2_per_molecule  if P > P_cloud_bar
//   0                       otherwise
// The transmission model multiplies sigma * VMR * n_total to get
// cm⁻¹, which becomes very large when sigma is large.
class CloudDeckOpacity final : public OpacityKernel {
 public:
  // P_cloud_bar : cloud-top pressure (atmosphere becomes opaque at higher P)
  // sigma_cm2   : effective per-cloud-particle cross section (use a
  //               large value, e.g. 1e-18 cm², to make the deck opaque)
  CloudDeckOpacity(std::string species_key,
                   double P_cloud_bar, double sigma_cm2);

  const std::string& species_key() const noexcept override { return key_; }

  Tensor cross_section(const std::vector<double>& wavenumber_cm,
                       const std::vector<double>& T_k,
                       const std::vector<double>& P_bar) const override;

  double P_cloud_bar() const noexcept { return P_cloud_bar_; }

 private:
  std::string key_;
  double P_cloud_bar_;
  double sigma_cm2_;
};

}  // namespace argus
