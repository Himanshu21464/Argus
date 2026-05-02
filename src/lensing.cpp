#include "argus/lensing.hpp"

#include <algorithm>
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

// ─── SIE ──────────────────────────────────────────────────────────────

SIE::SIE(double einstein_radius_arcsec, double axis_ratio,
         double position_angle_rad, Vec2 centre)
    : theta_E_(einstein_radius_arcsec),
      q_(axis_ratio),
      phi_(position_angle_rad),
      centre_(centre) {
  if (!(theta_E_ > 0.0)) {
    throw std::invalid_argument("SIE: einstein_radius must be positive");
  }
  if (!(q_ > 0.0 && q_ <= 1.0)) {
    throw std::invalid_argument("SIE: axis_ratio must be in (0, 1]");
  }
}

Vec2 SIE::deflection(Vec2 theta) const {
  // Translate to lens-centred frame.
  const double dx = theta.x - centre_.x;
  const double dy = theta.y - centre_.y;

  // Rotate into the major-axis-aligned frame (x' along major axis).
  const double cphi = std::cos(phi_);
  const double sphi = std::sin(phi_);
  const double xp =  cphi * dx + sphi * dy;
  const double yp = -sphi * dx + cphi * dy;

  Vec2 ap;  // deflection in primed frame.

  // Series fall-back for q very close to 1 (avoid division by √(1-q²)).
  // Taylor to leading order: arctan(z)/z → 1, arctanh(z)/z → 1, so
  // α'_x/y → θ_E · √q · x'/y' / ψ; with q=1 this is θ_E · x'/y' / r.
  if (std::fabs(1.0 - q_) < 1.0e-6) {
    const double r = std::sqrt(xp * xp + yp * yp);
    if (r < 1.0e-15) return {0.0, 0.0};
    ap = {theta_E_ * xp / r, theta_E_ * yp / r};
  } else {
    const double q2  = q_ * q_;
    const double psi = std::sqrt(q2 * xp * xp + yp * yp);
    if (psi < 1.0e-15) return {0.0, 0.0};
    const double qp  = std::sqrt(1.0 - q2);
    const double pre = theta_E_ * std::sqrt(q_) / qp;
    ap.x = pre * std::atan( qp * xp / psi);
    ap.y = pre * std::atanh(qp * yp / psi);
  }

  // Rotate back to sky frame.
  return {cphi * ap.x - sphi * ap.y,
          sphi * ap.x + cphi * ap.y};
}

// ─── Generic image finder ─────────────────────────────────────────────
//
// Lens equation: F(θ) = θ - α(θ) - β = 0. Find every root of F in a
// search box. Strategy:
//   1. Evaluate |F|² on a `grid_n × grid_n` mesh.
//   2. For every cell whose corner residual is below a coarse cut,
//      seed a Newton iteration from the cell centre. The Jacobian
//      ∂F/∂θ is computed via central differences of α.
//   3. Newton converges quadratically near a root; cap iterations and
//      reject if residual stays above tol.
//   4. Deduplicate accepted roots within `dedup_tol`.
//   5. Magnification at each image is 1 / |det(∂F/∂θ)|.

namespace {

double det2(double a, double b, double c, double d) {
  return a * d - b * c;
}

// Central-difference Jacobian of α at θ. Returns the 2×2 entries
// (αx_x, αx_y, αy_x, αy_y) where αx_x = ∂α_x/∂θ_x etc.
struct Jac { double xx, xy, yx, yy; };

Jac alpha_jacobian(const Lens& lens, Vec2 theta, double h) {
  const Vec2 ax_p = lens.deflection({theta.x + h, theta.y});
  const Vec2 ax_m = lens.deflection({theta.x - h, theta.y});
  const Vec2 ay_p = lens.deflection({theta.x, theta.y + h});
  const Vec2 ay_m = lens.deflection({theta.x, theta.y - h});
  return {
    (ax_p.x - ax_m.x) / (2.0 * h),
    (ay_p.x - ay_m.x) / (2.0 * h),
    (ax_p.y - ax_m.y) / (2.0 * h),
    (ay_p.y - ay_m.y) / (2.0 * h),
  };
}

}  // namespace

std::vector<Image> find_images(const Lens& lens,
                               Vec2 beta,
                               double search_radius,
                               std::size_t grid_n,
                               double tol,
                               double dedup_tol,
                               std::size_t newton_max_iter) {
  if (!(search_radius > 0.0)) {
    throw std::invalid_argument("find_images: search_radius must be positive");
  }
  if (grid_n < 2) {
    throw std::invalid_argument("find_images: grid_n must be >= 2");
  }
  if (!(tol > 0.0)) {
    throw std::invalid_argument("find_images: tol must be positive");
  }
  if (!(dedup_tol > 0.0)) {
    throw std::invalid_argument("find_images: dedup_tol must be positive");
  }
  if (newton_max_iter == 0) {
    throw std::invalid_argument("find_images: newton_max_iter must be >= 1");
  }

  const double dx = 2.0 * search_radius / static_cast<double>(grid_n);
  const double h_jac = std::max(1.0e-7, 1.0e-4 * search_radius);

  // 1. Coarse grid — record |F|² and per-cell-corner sign of (F_x, F_y).
  std::vector<double> Fx(grid_n * grid_n);
  std::vector<double> Fy(grid_n * grid_n);
  std::vector<Vec2>   theta_grid(grid_n * grid_n);
  for (std::size_t j = 0; j <= grid_n; ++j) {
    for (std::size_t i = 0; i <= grid_n; ++i) {
      // Note: (grid_n+1)² corner mesh would be cleaner; use (grid_n)²
      // node-centred plus cell-centred seeds for simplicity.
      if (i >= grid_n || j >= grid_n) continue;
      const Vec2 t{
        -search_radius + dx * (static_cast<double>(i) + 0.5),
        -search_radius + dx * (static_cast<double>(j) + 0.5)
      };
      const Vec2 a = lens.deflection(t);
      const Vec2 F = {t.x - a.x - beta.x, t.y - a.y - beta.y};
      const std::size_t idx = j * grid_n + i;
      Fx[idx] = F.x;
      Fy[idx] = F.y;
      theta_grid[idx] = t;
    }
  }

  // 2. Seed Newton at every cell whose F-magnitude is "near zero"
  //    relative to the local grid scale. Loose threshold catches all
  //    candidate basins of attraction; deduplication after convergence
  //    collapses duplicates.
  const double seed_threshold = 4.0 * dx;
  std::vector<Image> roots;

  auto try_newton = [&](Vec2 seed) -> bool {
    Vec2 t = seed;
    for (std::size_t it = 0; it < newton_max_iter; ++it) {
      const Vec2 a = lens.deflection(t);
      const Vec2 F = {t.x - a.x - beta.x, t.y - a.y - beta.y};
      const double resid = std::sqrt(F.x * F.x + F.y * F.y);
      if (resid < tol) {
        // Magnification = 1 / |det(∂β/∂θ)|, where ∂β/∂θ = I - ∂α/∂θ.
        const Jac J = alpha_jacobian(lens, t, h_jac);
        const double det_beta = det2(1.0 - J.xx, -J.xy, -J.yx, 1.0 - J.yy);
        const double mu = 1.0 / std::fabs(det_beta);
        // Deduplicate.
        for (const auto& r : roots) {
          const double d = std::sqrt((r.theta.x - t.x) * (r.theta.x - t.x) +
                                     (r.theta.y - t.y) * (r.theta.y - t.y));
          if (d < dedup_tol) return false;
        }
        roots.push_back(Image{t, mu});
        return true;
      }
      // Newton step: J_F · Δθ = -F  with  J_F = I - ∂α/∂θ.
      const Jac J = alpha_jacobian(lens, t, h_jac);
      const double a11 = 1.0 - J.xx;
      const double a12 = -J.xy;
      const double a21 = -J.yx;
      const double a22 = 1.0 - J.yy;
      const double det = det2(a11, a12, a21, a22);
      if (std::fabs(det) < 1.0e-14) return false;  // Singular Jacobian.
      const double dxs = (-a22 * F.x + a12 * F.y) / det;
      const double dys = ( a21 * F.x - a11 * F.y) / det;
      t.x += dxs;
      t.y += dys;
      // Damp wild excursions to keep the iteration bounded.
      if (std::fabs(t.x) > 10.0 * search_radius ||
          std::fabs(t.y) > 10.0 * search_radius) return false;
    }
    return false;
  };

  for (std::size_t j = 0; j < grid_n; ++j) {
    for (std::size_t i = 0; i < grid_n; ++i) {
      const std::size_t idx = j * grid_n + i;
      const double resid = std::sqrt(Fx[idx] * Fx[idx] + Fy[idx] * Fy[idx]);
      if (resid < seed_threshold) {
        try_newton(theta_grid[idx]);
      }
    }
  }

  // Sort by image x for stable identification by callers.
  std::sort(roots.begin(), roots.end(),
            [](const Image& a, const Image& b) {
              if (a.theta.x != b.theta.x) return a.theta.x < b.theta.x;
              return a.theta.y < b.theta.y;
            });
  return roots;
}

}  // namespace argus::lensing
