// Tests for the M4 strong-lensing pass.
//
// Validates the SIS lens model + image solver against analytic
// closed-form expressions. The substrate claim is that the same
// Argus IR + Retrieval API used for atmospheric retrieval ALSO
// works for lens-parameter inference; this file proves the physics
// layer is correct first.

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "argus/argus.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

}  // namespace

int main() {
  using namespace argus::lensing;

  // ─── 1. SIS deflection magnitude is constant = θ_E for any radius. ──
  {
    SIS lens(/*theta_E=*/1.0);
    for (double r : {0.5, 1.0, 2.0, 10.0}) {
      Vec2 a = lens.deflection({r, 0.0});
      const double mag = std::sqrt(a.x * a.x + a.y * a.y);
      assert(close(mag, 1.0, 1.0e-12));
      assert(close(a.x, 1.0, 1.0e-12));    // points outward in +x
      assert(std::fabs(a.y) < 1.0e-15);
    }
  }

  // ─── 2. Deflection points radially outward in arbitrary direction. ──
  {
    SIS lens(2.0);
    Vec2 theta{3.0, 4.0};                    // radius 5
    Vec2 a = lens.deflection(theta);
    // Magnitude = θ_E = 2; direction = (3, 4) / 5
    assert(close(a.x, 2.0 * 3.0 / 5.0, 1.0e-12));
    assert(close(a.y, 2.0 * 4.0 / 5.0, 1.0e-12));
  }

  // ─── 3. Source inside Einstein radius gives 2 images with the lens
  //     equation closing on each. β = θ - α(θ). ────────────────────────
  {
    SIS lens(1.0);
    Vec2 beta{0.3, 0.0};                     // inside Einstein radius
    auto imgs = sis_images(lens, beta);
    assert(imgs.size() == 2);

    for (const auto& im : imgs) {
      Vec2 b = lens_equation(lens, im.theta);
      assert(close(b.x, beta.x, 1.0e-10));
      assert(close(b.y, beta.y, 1.0e-10));
    }

    // Outer image at θ=1.3, magnification = 1.3/0.3
    assert(close(imgs[0].theta.x, 1.3, 1.0e-12));
    assert(close(imgs[0].theta.y, 0.0, 1.0e-12));
    assert(close(imgs[0].magnification, 1.3 / 0.3, 1.0e-12));

    // Inner image at θ=-0.7, magnification = 0.7/0.3
    assert(close(imgs[1].theta.x, -0.7, 1.0e-12));
    assert(close(imgs[1].theta.y, 0.0, 1.0e-12));
    assert(close(imgs[1].magnification, 0.7 / 0.3, 1.0e-12));
  }

  // ─── 4. Source outside Einstein radius → 1 image, also closing. ────
  {
    SIS lens(1.0);
    Vec2 beta{2.0, 0.0};
    auto imgs = sis_images(lens, beta);
    assert(imgs.size() == 1);
    Vec2 b = lens_equation(lens, imgs[0].theta);
    assert(close(b.x, 2.0, 1.0e-12));
    assert(close(b.y, 0.0, 1.0e-12));
    assert(close(imgs[0].theta.x, 3.0, 1.0e-12));
    assert(close(imgs[0].magnification, 3.0 / 2.0, 1.0e-12));
  }

  // ─── 5. Off-axis source: image positions trace a line through the
  //     lens centre + source point. ───────────────────────────────────
  {
    SIS lens(1.0);
    Vec2 beta{0.4, 0.3};                       // |β| = 0.5
    auto imgs = sis_images(lens, beta);
    assert(imgs.size() == 2);

    // For each image, θ is collinear with β and the lens centre (origin).
    for (const auto& im : imgs) {
      // Cross-product θ × β should be zero (collinearity test).
      assert(std::fabs(im.theta.x * beta.y - im.theta.y * beta.x) < 1.0e-12);
    }

    // Lens-equation closure at each image.
    for (const auto& im : imgs) {
      Vec2 b = lens_equation(lens, im.theta);
      assert(close(b.x, beta.x, 1.0e-10));
      assert(close(b.y, beta.y, 1.0e-10));
    }

    // Total magnification for SIS at β < θ_E:
    //   μ_total = (β+θ_E)/β + (θ_E-β)/β = 2 θ_E / β = 2 / 0.5 = 4
    const double mu_total = imgs[0].magnification + imgs[1].magnification;
    assert(close(mu_total, 4.0, 1.0e-12));
  }

  // ─── 6. Off-centre lens: images are computed in the same global
  //     coordinate frame as β, with the lens centre as the offset. ────
  {
    SIS lens(1.0, /*centre=*/{5.0, -3.0});
    Vec2 beta{5.4, -3.0};                       // 0.4 east of centre
    auto imgs = sis_images(lens, beta);
    assert(imgs.size() == 2);

    // Outer image at lens.centre + (1.4, 0)
    assert(close(imgs[0].theta.x, 5.0 + 1.4, 1.0e-10));
    assert(close(imgs[0].theta.y, -3.0,        1.0e-10));
    // Inner image at lens.centre + (-0.6, 0)
    assert(close(imgs[1].theta.x, 5.0 - 0.6, 1.0e-10));
    assert(close(imgs[1].theta.y, -3.0,        1.0e-10));
  }

  // ─── 7. Source exactly on the lens centre: empty (degenerate ring). ─
  {
    SIS lens(1.0);
    auto imgs = sis_images(lens, {0.0, 0.0});
    assert(imgs.empty());
  }

  // ─── 8. Bad inputs throw. ───────────────────────────────────────────
  {
    bool threw = false;
    try { SIS(0.0); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { SIS(-1.0); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  return 0;
}
