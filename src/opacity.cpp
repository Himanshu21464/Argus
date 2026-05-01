#include "argus/opacity.hpp"

#include <stdexcept>

namespace argus {

GreyOpacity::GreyOpacity(std::string species_key, double sigma_cm2)
    : key_(std::move(species_key)), sigma_cm2_(sigma_cm2) {}

Tensor GreyOpacity::cross_section(const std::vector<double>& wavenumber_cm,
                                  const std::vector<double>& T_k,
                                  const std::vector<double>& P_bar) const {
  const std::size_t nT = T_k.size();
  const std::size_t nP = P_bar.size();
  const std::size_t nW = wavenumber_cm.size();
  Tensor out({nT, nP, nW});
  for (std::size_t i = 0; i < out.size(); ++i) {
    out[i] = sigma_cm2_;
  }
  return out;
}

// ─── CloudDeckOpacity ─────────────────────────────────────────────────

CloudDeckOpacity::CloudDeckOpacity(std::string species_key,
                                   double P_cloud_bar, double sigma_cm2)
    : key_(std::move(species_key)),
      P_cloud_bar_(P_cloud_bar),
      sigma_cm2_(sigma_cm2) {
  if (P_cloud_bar_ <= 0.0) {
    throw std::invalid_argument(
        "CloudDeckOpacity: P_cloud_bar must be positive");
  }
  if (sigma_cm2_ < 0.0) {
    throw std::invalid_argument(
        "CloudDeckOpacity: sigma_cm2 must be >= 0");
  }
}

Tensor CloudDeckOpacity::cross_section(
    const std::vector<double>& wavenumber_cm,
    const std::vector<double>& T_k,
    const std::vector<double>& P_bar) const {
  const std::size_t nT = T_k.size();
  const std::size_t nP = P_bar.size();
  const std::size_t nW = wavenumber_cm.size();
  Tensor out({nT, nP, nW});
  for (std::size_t iT = 0; iT < nT; ++iT) {
    for (std::size_t iP = 0; iP < nP; ++iP) {
      const double sigma_layer =
          (P_bar[iP] >= P_cloud_bar_) ? sigma_cm2_ : 0.0;
      for (std::size_t w = 0; w < nW; ++w) {
        out[(iT * nP + iP) * nW + w] = sigma_layer;
      }
    }
  }
  return out;
}

}  // namespace argus
