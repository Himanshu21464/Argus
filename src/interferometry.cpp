#include "argus/interferometry.hpp"

#include <cmath>
#include <stdexcept>

namespace argus::interferometry {

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;

}  // namespace

Visibility predict_visibility(const PointSource& src, UVPoint uv) {
  const double phase = -kTwoPi * (uv.u * src.l + uv.v * src.m);
  return {src.flux * std::cos(phase), src.flux * std::sin(phase)};
}

Visibility predict_visibility(const GaussianSource& src, UVPoint uv) {
  // FT of a normalised 2-D Gaussian (∫I = F, σ in image plane) is a
  // Gaussian in the UV plane with width 1/(2π σ):
  //   V_g(u,v) = F · exp(-2π² σ² (u² + v²)) · exp(-2πi (u·l + v·m))
  constexpr double kHalfTwoPiSq = 0.5 * kTwoPi * kTwoPi;     // 2π²
  const double r2 = uv.u * uv.u + uv.v * uv.v;
  const double envelope = src.flux *
      std::exp(-kHalfTwoPiSq * src.sigma * src.sigma * r2);
  const double phase = -kTwoPi * (uv.u * src.l + uv.v * src.m);
  return {envelope * std::cos(phase), envelope * std::sin(phase)};
}

std::vector<Visibility> predict_visibilities(
    const std::vector<PointSource>& sources,
    const std::vector<UVPoint>&     uv_points) {
  std::vector<Visibility> out(uv_points.size(), Visibility{0.0, 0.0});
  for (std::size_t i = 0; i < uv_points.size(); ++i) {
    Visibility sum{0.0, 0.0};
    for (const auto& s : sources) {
      const Visibility v = predict_visibility(s, uv_points[i]);
      sum.real += v.real;
      sum.imag += v.imag;
    }
    out[i] = sum;
  }
  return out;
}

std::vector<Visibility> predict_visibilities(
    const std::vector<GaussianSource>& sources,
    const std::vector<UVPoint>&        uv_points) {
  for (const auto& s : sources) {
    if (s.sigma < 0.0) {
      throw std::invalid_argument(
          "predict_visibilities: GaussianSource sigma must be >= 0");
    }
  }
  std::vector<Visibility> out(uv_points.size(), Visibility{0.0, 0.0});
  for (std::size_t i = 0; i < uv_points.size(); ++i) {
    Visibility sum{0.0, 0.0};
    for (const auto& s : sources) {
      const Visibility v = predict_visibility(s, uv_points[i]);
      sum.real += v.real;
      sum.imag += v.imag;
    }
    out[i] = sum;
  }
  return out;
}

std::vector<UVPoint> uv_coverage_snapshot(
    const std::vector<double>& east_m,
    const std::vector<double>& north_m,
    double                     wavelength_m) {
  if (east_m.size() != north_m.size()) {
    throw std::invalid_argument(
        "uv_coverage_snapshot: east/north antenna lists must match in length");
  }
  if (east_m.size() < 2) {
    throw std::invalid_argument(
        "uv_coverage_snapshot: need at least 2 antennas");
  }
  if (!(wavelength_m > 0.0)) {
    throw std::invalid_argument(
        "uv_coverage_snapshot: wavelength_m must be positive");
  }

  const std::size_t n = east_m.size();
  std::vector<UVPoint> out;
  out.reserve(n * (n - 1) / 2);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      const double bx = east_m[j]  - east_m[i];
      const double by = north_m[j] - north_m[i];
      out.push_back({bx / wavelength_m, by / wavelength_m});
    }
  }
  return out;
}

}  // namespace argus::interferometry
