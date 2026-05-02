#include "argus/lensing.hpp"

#include <cmath>
#include <stdexcept>

namespace argus::lensing {

// ─── SIS ──────────────────────────────────────────────────────────────

SIS::SIS(double einstein_radius_arcsec, Vec2 centre)
    : theta_E_(einstein_radius_arcsec), centre_(centre) {
  if (!(theta_E_ > 0.0)) {
    throw std::invalid_argument("SIS: einstein_radius must be positive");
  }
}

Vec2 SIS::deflection(Vec2 theta) const {
  const Vec2 d = theta - centre_;
  const double r = std::sqrt(d.x * d.x + d.y * d.y);
  if (r < 1.0e-15) return {0.0, 0.0};
  return {theta_E_ * d.x / r, theta_E_ * d.y / r};
}

// ─── Lens equation ────────────────────────────────────────────────────

Vec2 lens_equation(const Lens& lens, Vec2 theta) {
  return theta - lens.deflection(theta);
}

// ─── SIS image solver ─────────────────────────────────────────────────
//
// In a 1-D coordinate aligned with the source-lens direction:
//   β > 0 (source above the lens), let θ be the image-plane radial
//   position measured from lens centre along that line.
//   Lens equation: β = θ - θ_E · sign(θ).
//
// Solutions:
//     θ_+ = β + θ_E     (outer image, always present)
//     θ_- = β - θ_E     (inner image if β < θ_E; appears on opposite
//                        side of the lens from the source)
//
// Magnifications for SIS:
//     μ_± = |θ_± / (θ_± - θ_E·sign(θ_±))|
//         = |θ_± / β|        (since β = θ - θ_E·sign(θ))
//
// Total magnification across both images equals |θ_+ + θ_-| / |β|
// = 2 θ_+ / β = (β + θ_E)/β + |β - θ_E|/|β| (sum of absolute values).

std::vector<Image> sis_images(const SIS& lens, Vec2 beta) {
  const Vec2  c     = lens.centre();
  const double te   = lens.einstein_radius();
  const Vec2  d     = beta - c;
  const double br   = std::sqrt(d.x * d.x + d.y * d.y);

  std::vector<Image> out;

  if (br < 1.0e-15) {
    // Source exactly at the lens: degenerate Einstein ring.
    return out;
  }

  // Unit vector toward the source from the lens centre.
  const double ux = d.x / br;
  const double uy = d.y / br;

  // Outer image: same side as source, distance (br + te) from centre.
  const double r_plus = br + te;
  out.push_back(Image{
    .theta = {c.x + r_plus * ux, c.y + r_plus * uy},
    .magnification = r_plus / br
  });

  // Inner image: opposite side, distance |te - br| from centre.
  // Exists only when source is inside the Einstein radius (br < te).
  if (br < te) {
    const double r_minus = te - br;
    out.push_back(Image{
      .theta = {c.x - r_minus * ux, c.y - r_minus * uy},
      .magnification = r_minus / br
    });
  }

  return out;
}

}  // namespace argus::lensing
