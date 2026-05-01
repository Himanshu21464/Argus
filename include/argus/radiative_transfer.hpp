#pragma once

#include <memory>
#include <vector>

#include "atmosphere.hpp"
#include "opacity.hpp"

namespace argus {

// Result of a single forward call. Production version carries a tape handle
// for autograd; the M1 stub carries spectrum samples only.
struct Spectrum {
  std::vector<double> wavenumber_cm;  // [n_wave]
  std::vector<double> values;         // [n_wave] — transit depth or flux
};

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
