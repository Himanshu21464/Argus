#pragma once

#include <memory>
#include <utility>
#include <vector>

namespace argus::lensing {

// Strong gravitational lensing forward models.
//
// The lens equation in the thin-lens approximation:
//
//     β = θ - α(θ)
//
// where β is the (unobservable) source position on the sky, θ is an
// image position on the lens plane, and α(θ) is the deflection angle
// produced by the lens mass distribution.
//
// M4's substrate claim: the same Argus IR + Retrieval API used for
// exoplanet atmospheric retrieval also applies to lens-parameter
// inference from observed image positions / arc shapes — only the
// physics layer (this file) is different.

struct Vec2 {
  double x = 0.0;
  double y = 0.0;
};

inline Vec2 operator+(Vec2 a, Vec2 b) { return {a.x + b.x, a.y + b.y}; }
inline Vec2 operator-(Vec2 a, Vec2 b) { return {a.x - b.x, a.y - b.y}; }
inline Vec2 operator*(double s, Vec2 a) { return {s * a.x, s * a.y}; }
inline Vec2 operator*(Vec2 a, double s) { return {s * a.x, s * a.y}; }

class Lens {
 public:
  virtual ~Lens() = default;
  // Deflection angle α(θ) (radians, in the same coordinate convention
  // as θ). For axisymmetric lenses α points radially outward.
  virtual Vec2 deflection(Vec2 theta) const = 0;
};

// Singular Isothermal Sphere — analytic; the workhorse first-order
// lens model. Constant velocity dispersion produces a flat rotation
// curve; the projected mass interior to angular radius θ scales
// linearly with θ.
//
// Deflection magnitude is constant: |α| = θ_E (the Einstein radius).
// Lens equation reduces to a 1-D problem in the source-plane direction
// of β.
class SIS final : public Lens {
 public:
  // einstein_radius_arcsec is θ_E; centre is the lens position on
  // the sky.
  SIS(double einstein_radius_arcsec,
      Vec2 centre = {0.0, 0.0});

  Vec2 deflection(Vec2 theta) const override;

  double einstein_radius() const noexcept { return theta_E_; }
  Vec2   centre()          const noexcept { return centre_; }

 private:
  double theta_E_;
  Vec2   centre_;
};

// β = θ - α(θ). Forward map.
Vec2 lens_equation(const Lens& lens, Vec2 theta);

// One image of a point source.
struct Image {
  Vec2   theta;          // image position on the lens plane
  double magnification;  // |β / dβ| — total magnification factor
};

// Closed-form image solver for an SIS. Returns 1 image if β > θ_E
// (source outside the Einstein radius), 2 images otherwise. The
// caller's β is in the SAME coordinate frame as the lens centre.
std::vector<Image> sis_images(const SIS& lens, Vec2 beta);

}  // namespace argus::lensing
