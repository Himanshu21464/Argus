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

  // Lensing (scaled gravitational) potential ψ(θ) such that
  // α(θ) = ∇ψ(θ). Required for the Fermat potential and time delays.
  // Default implementation returns 0 — concrete lenses override.
  virtual double potential(Vec2 theta) const { (void)theta; return 0.0; }
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

  Vec2   deflection(Vec2 theta) const override;
  double potential (Vec2 theta) const override;     // ψ_SIS = θ_E · |θ - centre|

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

// Singular Isothermal Ellipsoid — the standard elliptical-galaxy /
// cluster-scale lens model. Same surface-mass-density profile as SIS
// (κ ∝ 1/ξ where ξ is the elliptical radius) but with axis ratio
// q ∈ (0, 1] and a position-angle rotation φ.
//
// q = 1 reduces exactly to SIS. For q < 1 the lens produces 2- or
// 4-image configurations depending on whether the source lies inside
// the tangential caustic (4 images) or outside (2 images). The
// closed-form deflection is from Kormann, Schneider & Bartelmann
// (1994) eqs (40-41), evaluated in the major-axis-aligned frame:
//
//     ψ(x, y) = √(q² x² + y²)
//     α_x =  (θ_E √q / √(1-q²)) · arctan( √(1-q²)·x / ψ )
//     α_y =  (θ_E √q / √(1-q²)) · arctanh(√(1-q²)·y / ψ )
//
// Numerically stable via series expansion as q → 1.
class SIE final : public Lens {
 public:
  // einstein_radius_arcsec: SIS-equivalent θ_E (the q=1 limit).
  // axis_ratio: q ∈ (0, 1]; 1 ⇒ circular SIS.
  // position_angle_rad: angle of the major axis from +x, CCW.
  // centre: lens position on the sky.
  SIE(double einstein_radius_arcsec,
      double axis_ratio,
      double position_angle_rad = 0.0,
      Vec2 centre = {0.0, 0.0});

  Vec2   deflection(Vec2 theta) const override;
  double potential (Vec2 theta) const override;

  double einstein_radius() const noexcept { return theta_E_; }
  double axis_ratio()     const noexcept { return q_; }
  double position_angle() const noexcept { return phi_; }
  Vec2   centre()         const noexcept { return centre_; }

 private:
  double theta_E_;
  double q_;
  double phi_;
  Vec2   centre_;
};

// Numerical image solver for an arbitrary `Lens` via a coarse grid
// search followed by Newton refinement. Suitable for SIE and other
// non-axisymmetric models that have no closed-form image solver.
//
// Searches a square box of size 2·search_radius centred on the lens
// at `grid_n × grid_n` resolution; any cell whose lens-equation
// residual lies near zero is refined with Newton iteration until
// |β_model − β_target| < tol. Duplicates within `dedup_tol` are
// merged.
//
// Magnification per image is the inverse of the lens-equation
// Jacobian determinant at that image position, computed via central
// differences of the deflection.
//
// Typical configuration for an SIE with θ_E ≈ 1: search_radius = 3,
// grid_n = 60, tol = 1.0e-9, dedup_tol = 1.0e-4.
std::vector<Image> find_images(const Lens& lens,
                               Vec2 beta,
                               double search_radius,
                               std::size_t grid_n = 60,
                               double tol = 1.0e-9,
                               double dedup_tol = 1.0e-4,
                               std::size_t newton_max_iter = 50);

// External tidal shear field — the standard ingredient in real
// strong-lensing models. Combines with an SIE (or any other lens)
// via CompoundLens to model perturbations from large-scale structure
// or nearby galaxies.
//
// Convention: shear matrix [[γ_1, γ_2], [γ_2, -γ_1]]. The deflection
// is then α_x = γ_1 θ_x + γ_2 θ_y, α_y = γ_2 θ_x - γ_1 θ_y, and the
// scalar potential is ψ = 0.5 [γ_1 (θ_x² - θ_y²) + 2 γ_2 θ_x θ_y].
// The shear is centred at the origin by convention; compose with a
// translated lens for a different effective shear pivot.
class ExternalShear final : public Lens {
 public:
  ExternalShear(double gamma1, double gamma2);

  Vec2   deflection(Vec2 theta) const override;
  double potential (Vec2 theta) const override;

  double gamma1() const noexcept { return gamma1_; }
  double gamma2() const noexcept { return gamma2_; }

 private:
  double gamma1_;
  double gamma2_;
};

// Compound lens: deflection / potential of the sum of multiple lens
// components. The owning std::shared_ptr design makes the components
// safe to share between models (e.g. fix a galaxy halo while varying
// shear and source).
class CompoundLens final : public Lens {
 public:
  CompoundLens() = default;
  explicit CompoundLens(std::vector<std::shared_ptr<const Lens>> components);

  void add(std::shared_ptr<const Lens> lens);
  std::size_t size() const noexcept { return components_.size(); }

  Vec2   deflection(Vec2 theta) const override;
  double potential (Vec2 theta) const override;

 private:
  std::vector<std::shared_ptr<const Lens>> components_;
};

// Fermat potential (the time-arrival surface, modulo a global
// constant): τ(θ; β) = 0.5 · |θ - β|² - ψ(θ). Critical points of τ
// are images of β; the values of τ at images give image arrival times
// up to a multiplicative cosmological factor.
double fermat_potential(const Lens& lens, Vec2 theta, Vec2 beta);

// Dimensionless time-delay surface difference Δτ between two images
// of the same source β. For two SIS images Δτ = 2 · θ_E · β
// (closed form). To convert to a physical time delay in days:
//
//     Δt[days] = (1 + z_l) · D_Δt[Mpc] / c[Mpc/day] · Δτ[arcsec²]
//                · (4.848e-6)²
//
// where D_Δt is the time-delay distance and the (4.848e-6)² factor
// converts arcsec² to rad². The cosmology is intentionally pushed
// out to the caller so this kernel is unit-clean.
double time_delay_arcsec2(const Lens& lens, Vec2 theta_a, Vec2 theta_b,
                          Vec2 beta);

}  // namespace argus::lensing
