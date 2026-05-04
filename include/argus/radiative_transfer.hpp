#pragma once

#include <memory>
#include <vector>

#include "atmosphere.hpp"
#include "opacity.hpp"

namespace argus {

// Result of a single forward call — sample values on the wavenumber grid.
// Reverse-mode autograd lives behind `argus/ad.hpp` (Wengert tape) for
// HMC and other gradient-driven samplers.
struct Spectrum {
  std::vector<double> wavenumber_cm;  // [n_wave]
  std::vector<double> values;         // [n_wave] — transit depth or flux
};

// Linearly-spaced wavenumber grid, low and high inclusive. Keeps the
// public-API code samples (and the per-band convenience tests) one line
// shorter than open-coding the loop, and makes the integration contract
// explicit: spectra in argus are sampled at equispaced wavenumbers.
std::vector<double> make_grid(double low_cm, double high_cm, std::size_t n);

// Differentiable transmission-spectrum forward model.
// Computes (R_lambda / R_star)^2 from the layered atmosphere, integrating the
// chord opacity through each impact parameter.
class TransmissionModel {
 public:
  TransmissionModel() = default;

  void add_opacity(std::shared_ptr<OpacityKernel> kernel);

  Spectrum forward(const Atmosphere& atmosphere,
                   const std::vector<double>& wavenumber_cm) const;

 private:
  std::vector<std::shared_ptr<OpacityKernel>> opacities_;
};

}  // namespace argus
