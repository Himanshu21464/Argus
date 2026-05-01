#include "argus/radiative_transfer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace argus {

namespace {

// Boltzmann constant in SI units, for converting (P, T) -> number density.
constexpr double kBoltzmannSI = 1.380649e-23;  // J / K

// Number density of gas at (P_bar, T_k) in molecules / cm^3.
double number_density_cm3(double P_bar, double T_k) {
  const double P_pa = P_bar * 1.0e5;
  const double n_per_m3 = P_pa / (kBoltzmannSI * T_k);
  return n_per_m3 * 1.0e-6;
}

}  // namespace

void TransmissionModel::add_opacity(std::shared_ptr<OpacityKernel> kernel) {
  if (!kernel) {
    throw std::invalid_argument("TransmissionModel::add_opacity: null kernel");
  }
  opacities_.push_back(std::move(kernel));
}

Spectrum TransmissionModel::forward(
    const Atmosphere& atm, const std::vector<double>& wavenumber_cm) const {
  atm.validate();

  // M1 stub: chord-integrated optical depth in the simplest possible
  // approximation (uniform path length per layer = scale-height proxy).
  // M3 will replace this with the proper geometric chord integral and
  // the autograd tape.
  const std::size_t nW = wavenumber_cm.size();
  const std::size_t nL = atm.num_layers();

  // Spawn dummy single-point T/P views so existing OpacityKernel signature
  // works unchanged across layers.
  Spectrum out;
  out.wavenumber_cm = wavenumber_cm;
  out.values.assign(nW, 0.0);

  for (std::size_t i_layer = 0; i_layer < nL; ++i_layer) {
    const double T = atm.temperature_k[i_layer];
    const double P = atm.pressure_bar[i_layer];
    const double n_total = number_density_cm3(P, T);
    const std::vector<double> T_v{T};
    const std::vector<double> P_v{P};

    for (std::size_t k = 0; k < opacities_.size(); ++k) {
      const auto& kernel = *opacities_[k];

      // Find species index in the atmosphere by key.
      std::size_t species_idx = atm.num_species();
      for (std::size_t s = 0; s < atm.num_species(); ++s) {
        if (atm.species[s].key == kernel.species_key()) {
          species_idx = s;
          break;
        }
      }
      if (species_idx == atm.num_species()) {
        continue;  // species not present in this atmosphere
      }
      const double vmr = atm.mixing_ratios.at(i_layer, species_idx);
      const double n_species = vmr * n_total;
      const Tensor sigma = kernel.cross_section(wavenumber_cm, T_v, P_v);

      // Path length proxy: 1 km per layer in cm. M3 derives this from
      // hydrostatic geometry.
      constexpr double path_cm = 1.0e5;
      for (std::size_t w = 0; w < nW; ++w) {
        out.values[w] += sigma[w] * n_species * path_cm;
      }
    }
  }

  // Convert optical depth -> transit depth proxy (1 - exp(-tau)).
  for (std::size_t w = 0; w < nW; ++w) {
    out.values[w] = 1.0 - std::exp(-out.values[w]);
  }
  return out;
}

}  // namespace argus
