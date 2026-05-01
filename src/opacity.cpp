#include "argus/opacity.hpp"

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

}  // namespace argus
