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

std::vector<double> guillot_profile(const std::vector<double>& pressures_bar,
                                    double T_int_k, double T_irr_k,
                                    double gamma,
                                    double kappa_IR_cm2_per_g,
                                    double gravity_si) {
  if (T_int_k < 0.0 || T_irr_k < 0.0 || gamma <= 0.0 ||
      kappa_IR_cm2_per_g <= 0.0 || gravity_si <= 0.0) {
    throw std::invalid_argument("guillot_profile: parameters must be positive "
                                 "(T_int and T_irr may be 0)");
  }
  // Guillot 2010 eq. 29:
  //   T⁴(τ) = (3/4) T_int⁴ · (2/3 + τ)
  //         + (3/4) T_irr⁴ · [2/3 + 1/(γ√3) + (γ/√3 - 1/(γ√3)) e^(-γτ√3)]
  // with τ = κ_IR · P / g.  Convert P from bar to SI and κ from cm²/g
  // to m²/kg (multiply by 0.1).
  const double kappa_si = kappa_IR_cm2_per_g * 0.1;     // m²/kg
  const double sqrt3    = std::sqrt(3.0);
  const double T_int4   = std::pow(T_int_k, 4.0);
  const double T_irr4   = std::pow(T_irr_k, 4.0);

  std::vector<double> T(pressures_bar.size());
  for (std::size_t i = 0; i < pressures_bar.size(); ++i) {
    const double P_si = pressures_bar[i] * 1.0e5;        // bar -> Pa
    const double tau  = kappa_si * P_si / gravity_si;
    const double a    = 2.0 / 3.0;
    const double b    = 1.0 / (gamma * sqrt3);
    const double c    = (gamma / sqrt3) - b;
    const double f_tau = a + b + c * std::exp(-gamma * tau * sqrt3);
    const double T4   = 0.75 * T_int4 * (a + tau)
                      + 0.75 * T_irr4 * f_tau;
    T[i] = std::pow(T4, 0.25);
  }
  return T;
}

Atmosphere guillot(double T_int_k, double T_irr_k, double gamma,
                   double P_top_bar, double P_bot_bar,
                   std::size_t n_layers,
                   Species species, double mixing_ratio,
                   double kappa_IR_cm2_per_g,
                   double planet_gravity_si) {
  if (n_layers < 2) {
    throw std::invalid_argument("guillot: n_layers must be >= 2");
  }
  if (P_top_bar <= 0.0 || P_bot_bar <= P_top_bar) {
    throw std::invalid_argument("guillot: require 0 < P_top < P_bot");
  }

  Atmosphere a;
  a.species = {std::move(species)};
  a.pressure_bar.resize(n_layers);
  a.mixing_ratios = Tensor(std::vector<std::size_t>{n_layers, 1});
  a.planet_gravity_si = planet_gravity_si;

  const double log_top = std::log(P_top_bar);
  const double log_bot = std::log(P_bot_bar);
  const double dlog = (log_bot - log_top) / static_cast<double>(n_layers - 1);
  for (std::size_t i = 0; i < n_layers; ++i) {
    a.pressure_bar[i] = std::exp(log_top + dlog * static_cast<double>(i));
    a.mixing_ratios.at(i, 0) = mixing_ratio;
  }
  a.temperature_k = guillot_profile(a.pressure_bar, T_int_k, T_irr_k, gamma,
                                    kappa_IR_cm2_per_g, planet_gravity_si);
  return a;
}

}  // namespace argus
