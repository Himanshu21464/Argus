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
  double planet_radius_rj = 1.0;     // R_Jupiter
  double planet_gravity_si = 25.0;   // m/s^2
  // bulk mean molecular weight of the background gas (amu).
  // 2.3 is the standard H2/He hot-Jupiter value.
  double bulk_mmw_amu = 2.3;
  // host-star radius in solar radii, used to convert R_eff -> transit depth.
  double star_radius_rsun = 1.0;

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

// Guillot 2010 analytic temperature-pressure profile for an irradiated
// hot Jupiter. Based on Guillot (2010), A&A 520, A27, eq. (29):
//
//   T⁴(τ) = (3/4) T_int⁴ · (2/3 + τ)
//         + (3/4) T_irr⁴ · [2/3 + 1/(γ√3) + (γ/√3 - 1/(γ√3)) · exp(-γτ√3)]
//
// where τ is the IR optical depth (∝ pressure for a constant-opacity
// atmosphere): τ = κ_IR · P / g.
//
// For a typical hot Jupiter: T_int ≈ 100–500 K, T_irr ≈ 1000–2500 K,
// γ ≈ 0.1–1 (thermal IR / visible opacity ratio), kappa_IR ≈ 10⁻³–10⁻¹
// cm²/g.
//
// `pressures_bar` must be sorted top-first (ascending pressure).
// Returns T(K) at each layer, same shape as the input.
std::vector<double> guillot_profile(const std::vector<double>& pressures_bar,
                                    double T_int_k, double T_irr_k,
                                    double gamma_thermal_visible,
                                    double kappa_IR_cm2_per_g = 1.0e-2,
                                    double gravity_si = 25.0);

// Build a complete Atmosphere with a Guillot T-P profile and a single
// species. Logarithmic pressure grid from P_top to P_bot.
Atmosphere guillot(double T_int_k, double T_irr_k,
                   double gamma_thermal_visible,
                   double P_top_bar, double P_bot_bar,
                   std::size_t n_layers,
                   Species species, double mixing_ratio,
                   double kappa_IR_cm2_per_g = 1.0e-2,
                   double planet_gravity_si = 25.0);

}  // namespace argus
