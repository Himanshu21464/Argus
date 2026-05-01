#include "argus/atmosphere.hpp"

#include <cmath>
#include <stdexcept>

namespace argus {

void Atmosphere::validate() const {
  if (pressure_bar.size() != temperature_k.size()) {
    throw std::invalid_argument(
        "Atmosphere: pressure and temperature must have equal layer count");
  }
  if (mixing_ratios.rank() != 2) {
    throw std::invalid_argument(
        "Atmosphere: mixing_ratios must be a 2-D tensor");
  }
  if (mixing_ratios.shape()[0] != pressure_bar.size()) {
    throw std::invalid_argument(
        "Atmosphere: mixing_ratios first dim must equal n_layers");
  }
  if (mixing_ratios.shape()[1] != species.size()) {
    throw std::invalid_argument(
        "Atmosphere: mixing_ratios second dim must equal n_species");
  }
  for (std::size_t i = 0; i < mixing_ratios.size(); ++i) {
    const double v = mixing_ratios[i];
    if (v < 0.0 || v > 1.0) {
      throw std::invalid_argument(
          "Atmosphere: mixing ratio outside [0, 1]");
    }
  }
}

Atmosphere isothermal(double T_k, double P_top_bar, double P_bot_bar,
                      std::size_t n_layers, Species species,
                      double mixing_ratio) {
  if (n_layers < 2) {
    throw std::invalid_argument("isothermal: n_layers must be >= 2");
  }
  if (P_top_bar <= 0.0 || P_bot_bar <= P_top_bar) {
    throw std::invalid_argument("isothermal: require 0 < P_top < P_bot");
  }

  Atmosphere a;
  a.species = {std::move(species)};
  a.pressure_bar.resize(n_layers);
  a.temperature_k.assign(n_layers, T_k);
  a.mixing_ratios = Tensor(std::vector<std::size_t>{n_layers, 1});

  const double log_top = std::log(P_top_bar);
  const double log_bot = std::log(P_bot_bar);
  const double dlog = (log_bot - log_top) / static_cast<double>(n_layers - 1);
  for (std::size_t i = 0; i < n_layers; ++i) {
    a.pressure_bar[i] = std::exp(log_top + dlog * static_cast<double>(i));
    a.mixing_ratios.at(i, 0) = mixing_ratio;
  }
  return a;
}

}  // namespace argus
