#include "argus/radiative_transfer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "argus/geometry.hpp"

namespace argus {

namespace {

// Number density of gas at (P_bar, T_k) in molecules / cm^3.
double number_density_cm3(double P_bar, double T_k) {
  const double P_pa = P_bar * 1.0e5;
  const double n_per_m3 = P_pa / (kBoltzmannSI * T_k);
  return n_per_m3 * 1.0e-6;
}

// Locate a species in an atmosphere by its HITRAN-style key. Returns
// `atm.num_species()` if absent.
std::size_t find_species_index(const Atmosphere& atm, const std::string& key) {
  for (std::size_t s = 0; s < atm.num_species(); ++s) {
    if (atm.species[s].key == key) return s;
  }
  return atm.num_species();
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

  const Geometry geom = build_geometry(atm);
  const std::size_t n  = atm.num_layers();
  const std::size_t nW = wavenumber_cm.size();

  // Layer edges (m). edge[i] is the upper boundary of layer i / the lower
  // boundary of layer i-1. Layer 0 is the top of the atmosphere, layer n-1
  // is the bottom. edge[0] is extrapolated above the topmost layer; edge[n]
  // is the planet "solid" radius.
  std::vector<double> edge(n + 1);
  edge[0] = 1.5 * geom.radius_m[0] - 0.5 * geom.radius_m[1];
  for (std::size_t i = 1; i < n; ++i) {
    edge[i] = 0.5 * (geom.radius_m[i - 1] + geom.radius_m[i]);
  }
  edge[n] = geom.planet_radius_m;

  // Pre-compute the per-layer absorption coefficient alpha[j][w] (cm^-1):
  //   alpha[j][w] = sum over species ( sigma_species(w) * n_species )
  // Cross-sections are evaluated once per layer per kernel.
  std::vector<std::vector<double>> alpha(n, std::vector<double>(nW, 0.0));
  for (std::size_t j = 0; j < n; ++j) {
    const double T = atm.temperature_k[j];
    const double P = atm.pressure_bar[j];
    const double n_total = number_density_cm3(P, T);

    for (const auto& kernel_ptr : opacities_) {
      const auto& kernel = *kernel_ptr;
      const std::size_t species_idx =
          find_species_index(atm, kernel.species_key());
      if (species_idx == atm.num_species()) continue;
      const double vmr = atm.mixing_ratios.at(j, species_idx);
      if (vmr <= 0.0) continue;
      const double n_species = vmr * n_total;

      const std::vector<double> T_v{T};
      const std::vector<double> P_v{P};
      const Tensor sig = kernel.cross_section(wavenumber_cm, T_v, P_v);
      for (std::size_t w = 0; w < nW; ++w) {
        alpha[j][w] += sig[w] * n_species;       // cm^-1
      }
    }
  }

  // Transit-radius integration. For each impact parameter b (sampled at
  // each layer centre), compute the optical depth tau(b, w) along the chord
  // and accumulate the effective area:
  //   R_eff^2(w) = R_p^2 + 2 * sum_b ( 1 - exp(-tau(b, w)) ) * b * db
  const double R_p     = geom.planet_radius_m;
  const double R_star  = atm.star_radius_rsun * kSolarRadiusM;
  const double R_star2 = R_star * R_star;

  Spectrum out;
  out.wavenumber_cm = wavenumber_cm;
  out.values.assign(nW, 0.0);

  for (std::size_t w = 0; w < nW; ++w) {
    double R_eff_sq = R_p * R_p;
    for (std::size_t i = 0; i < n; ++i) {
      const double b  = geom.radius_m[i];
      const double db = geom.thickness_m[i];
      double tau = 0.0;
      for (std::size_t j = 0; j < n; ++j) {
        const double r_outer = edge[j];
        const double r_inner = edge[j + 1];
        const double ds_m = chord_path_length(r_inner, r_outer, b);
        if (ds_m <= 0.0) continue;
        const double ds_cm = ds_m * 100.0;
        // Factor of 2 for in+out chord traversal.
        tau += 2.0 * ds_cm * alpha[j][w];
      }
      R_eff_sq += 2.0 * (1.0 - std::exp(-tau)) * b * db;
    }
    out.values[w] = R_eff_sq / R_star2;
  }
  return out;
}

std::vector<double> make_grid(double low_cm, double high_cm, std::size_t n) {
  if (n < 2) {
    throw std::invalid_argument("make_grid: n must be >= 2");
  }
  if (!(high_cm > low_cm)) {
    throw std::invalid_argument("make_grid: high_cm must be > low_cm");
  }
  std::vector<double> grid(n);
  const double step = (high_cm - low_cm) / static_cast<double>(n - 1);
  for (std::size_t i = 0; i < n; ++i) {
    grid[i] = low_cm + static_cast<double>(i) * step;
  }
  grid.back() = high_cm;  // avoid floating-point drift on the last sample
  return grid;
}

}  // namespace argus
