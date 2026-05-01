#pragma once

#include <string>
#include <vector>

#include "tensor.hpp"

namespace argus {

// One absorbing/emitting molecular species in the atmosphere.
// `key` is the canonical HITRAN-style identifier ("H2O", "CO2", "CH4", ...).
struct Species {
  std::string key;
  double molecular_weight_amu = 0.0;
};

// Layered 1-D atmosphere. The production version will support 3-D GCM-style
// columns; the IR is the same.
struct Atmosphere {
  // [n_layers] — pressure at the centre of each layer, in bar.
  std::vector<double> pressure_bar;
  // [n_layers] — temperature at the centre of each layer, in K.
  std::vector<double> temperature_k;
  // [n_layers, n_species] — volume mixing ratios.
  Tensor mixing_ratios;
  // ordered list of species; rows of mixing_ratios match this order.
  std::vector<Species> species;
  // bulk planet properties used by the geometry pass.
  double planet_radius_rj = 1.0;   // R_Jupiter
  double planet_gravity_si = 25.0; // m/s^2

  std::size_t num_layers() const { return pressure_bar.size(); }
  std::size_t num_species() const { return species.size(); }

  // Sanity checks: layer counts agree, mixing ratios in [0,1], etc.
  // Throws std::invalid_argument on failure.
  void validate() const;
};

// Convenience: build an isothermal hydrostatic atmosphere with a single
// species, useful for benchmarks and tests.
Atmosphere isothermal(double T_k, double P_top_bar, double P_bot_bar,
                      std::size_t n_layers, Species species,
                      double mixing_ratio);

}  // namespace argus
