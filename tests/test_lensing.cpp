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

  // ─── 9. SIE with q=1 reduces to SIS bit-exactly (or to 1e-12). ──────
  {
    SIS sis(1.5);
    SIE sie(1.5, /*q=*/1.0);
    for (Vec2 t : {Vec2{1.0, 0.5}, Vec2{-0.3, 2.0},
                   Vec2{0.5, -1.2}, Vec2{2.0, 0.0}}) {
      Vec2 a_sis = sis.deflection(t);
      Vec2 a_sie = sie.deflection(t);
      assert(close(a_sis.x, a_sie.x, 1.0e-12, 1.0e-12));
      assert(close(a_sis.y, a_sie.y, 1.0e-12, 1.0e-12));
    }
  }

  // ─── 10. SIE deflection on the major axis (y=0) — closed form. ─────
  // At (x, 0), ψ = q|x|, so:
  //   α_x = (θ_E √q / √(1-q²)) · arctan(√(1-q²)·sign(x)/q)
  //   α_y = 0
  {
    const double te = 1.2;
    const double q  = 0.7;
    SIE sie(te, q);
    Vec2 a = sie.deflection({1.0, 0.0});
    const double qp = std::sqrt(1.0 - q * q);
    const double exp_ax = te * std::sqrt(q) / qp * std::atan(qp / q);
    assert(close(a.x, exp_ax, 1.0e-12));
    assert(std::fabs(a.y) < 1.0e-15);
  }

  // ─── 11. SIE deflection on the minor axis (x=0) — closed form. ─────
  //   α_x = 0
  //   α_y = (θ_E √q / √(1-q²)) · arctanh(√(1-q²)·sign(y))
  {
    const double te = 1.2;
    const double q  = 0.7;
    SIE sie(te, q);
    Vec2 a = sie.deflection({0.0, 1.0});
    const double qp = std::sqrt(1.0 - q * q);
    const double exp_ay = te * std::sqrt(q) / qp * std::atanh(qp);
    assert(close(a.y, exp_ay, 1.0e-12));
    assert(std::fabs(a.x) < 1.0e-15);
  }

  // ─── 12. SIE point-symmetry: α(-θ) = -α(θ). ─────────────────────────
  {
    SIE sie(1.0, 0.6, /*phi=*/0.7);
    for (Vec2 t : {Vec2{0.7, 0.4}, Vec2{1.5, -0.3}, Vec2{-0.2, 0.9}}) {
      Vec2 ap = sie.deflection(t);
      Vec2 an = sie.deflection({-t.x, -t.y});
      assert(std::fabs(ap.x + an.x) < 1.0e-12);
      assert(std::fabs(ap.y + an.y) < 1.0e-12);
    }
  }

  // ─── 13. SIE rotation covariance: rotating the lens by π/2 and
  //     rotating θ by π/2 yields the π/2-rotated deflection. ──────────
  {
    SIE sie_0 (1.0, 0.7, /*phi=*/0.0);
    SIE sie_90(1.0, 0.7, /*phi=*/M_PI / 2.0);
    Vec2 t {1.0, 0.5};
    Vec2 t90{-0.5, 1.0};                          // 90° CCW rotation
    Vec2 a0  = sie_0.deflection(t);
    Vec2 a90 = sie_90.deflection(t90);
    Vec2 a0_rot{-a0.y, a0.x};                      // 90° CCW rotation
    assert(close(a0_rot.x, a90.x, 1.0e-12));
    assert(close(a0_rot.y, a90.y, 1.0e-12));
  }

  // ─── 14. SIE off-centre lens: deflection translation-invariant. ────
  {
    Vec2 c{0.4, -0.3};
    SIE sie_off(1.2, 0.7, 0.5, c);
    SIE sie_at0(1.2, 0.7, 0.5, {0.0, 0.0});
    Vec2 t_rel{0.8, 1.1};
    Vec2 a_off = sie_off.deflection({t_rel.x + c.x, t_rel.y + c.y});
    Vec2 a_at0 = sie_at0.deflection(t_rel);
    assert(close(a_off.x, a_at0.x, 1.0e-12));
    assert(close(a_off.y, a_at0.y, 1.0e-12));
  }

  // ─── 15. SIE 4-image config via find_images: source inside the
  //     tangential caustic produces 4 images, each closing the lens
  //     equation. ───────────────────────────────────────────────────
  {
    SIE sie(/*theta_E=*/1.0, /*q=*/0.7);
    Vec2 beta{0.05, 0.03};                          // near-cusp source
    auto imgs = find_images(sie, beta, /*radius=*/2.0, /*grid=*/100);
    assert(imgs.size() == 4);
    for (const auto& im : imgs) {
      Vec2 b = lens_equation(sie, im.theta);
      assert(close(b.x, beta.x, 1.0e-8, 1.0e-8));
      assert(close(b.y, beta.y, 1.0e-8, 1.0e-8));
      // Magnification finite and > 1 near a caustic.
      assert(std::isfinite(im.magnification));
      assert(im.magnification > 1.0);
    }
  }

  // ─── 16. SIE q=1 (i.e. SIS via SIE) with find_images recovers the
  //     2-image solution that sis_images returns. ────────────────────
  {
    SIE sie(1.0, 1.0);
    Vec2 beta{0.3, 0.0};
    auto imgs = find_images(sie, beta, 2.0, 80);
    assert(imgs.size() == 2);
    // Sorted by x: inner (-0.7, 0) then outer (1.3, 0).
    assert(close(imgs[0].theta.x, -0.7, 1.0e-6, 1.0e-6));
    assert(close(imgs[1].theta.x,  1.3, 1.0e-6, 1.0e-6));
    assert(std::fabs(imgs[0].theta.y) < 1.0e-6);
    assert(std::fabs(imgs[1].theta.y) < 1.0e-6);
    // SIS magnifications: 0.7/0.3, 1.3/0.3.
    assert(close(imgs[0].magnification, 0.7 / 0.3, 1.0e-3));
    assert(close(imgs[1].magnification, 1.3 / 0.3, 1.0e-3));
  }

  // ─── 17. SIE bad inputs throw. ──────────────────────────────────────
  {
    bool threw = false;
    try { SIE(0.0, 0.7); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { SIE(1.0, 0.0); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { SIE(1.0, 1.5); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { SIE(1.0, -0.5); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // ─── 18. find_images bad inputs throw. ──────────────────────────────
  {
    SIS sis(1.0);
    bool threw = false;
    try { (void)find_images(sis, {0.0, 0.0}, /*radius=*/0.0); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { (void)find_images(sis, {0.0, 0.0}, 1.0, /*grid_n=*/1); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { (void)find_images(sis, {0.0, 0.0}, 1.0, 10, /*tol=*/0.0); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // ─── 19. ExternalShear deflection closed form. ──────────────────────
  // α_x = γ_1 θ_x + γ_2 θ_y, α_y = γ_2 θ_x - γ_1 θ_y.
  {
    ExternalShear sh(/*g1=*/0.05, /*g2=*/-0.03);
    for (Vec2 t : {Vec2{1.0, 0.5}, Vec2{0.0, 1.0}, Vec2{-2.0, 0.7}}) {
      Vec2 a = sh.deflection(t);
      assert(close(a.x,  0.05 * t.x + (-0.03) * t.y, 1.0e-12));
      assert(close(a.y, -0.03 * t.x - ( 0.05) * t.y, 1.0e-12));
    }
  }

  // ─── 20. ExternalShear potential is gradient-consistent with α. ────
  {
    ExternalShear sh(0.06, 0.04);
    const double h = 1.0e-5;
    for (Vec2 t : {Vec2{0.7, -0.2}, Vec2{1.5, 1.0}, Vec2{-0.9, 0.4}}) {
      const double pxp = sh.potential({t.x + h, t.y});
      const double pxm = sh.potential({t.x - h, t.y});
      const double pyp = sh.potential({t.x, t.y + h});
      const double pym = sh.potential({t.x, t.y - h});
      const double gx = (pxp - pxm) / (2.0 * h);
      const double gy = (pyp - pym) / (2.0 * h);
      Vec2 a = sh.deflection(t);
      assert(close(gx, a.x, 1.0e-7, 1.0e-9));
      assert(close(gy, a.y, 1.0e-7, 1.0e-9));
    }
  }

  // ─── 21. CompoundLens sums deflection and potential. ───────────────
  {
    auto sis = std::make_shared<SIS>(1.0);
    auto sh  = std::make_shared<ExternalShear>(0.05, -0.02);
    CompoundLens c{{sis, sh}};
    for (Vec2 t : {Vec2{1.0, 0.5}, Vec2{-0.3, 1.5}, Vec2{0.7, -0.7}}) {
      Vec2 a_sum = sis->deflection(t);
      Vec2 a_sh  = sh->deflection(t);
      Vec2 a_c   = c.deflection(t);
      assert(close(a_c.x, a_sum.x + a_sh.x, 1.0e-12));
      assert(close(a_c.y, a_sum.y + a_sh.y, 1.0e-12));
      const double p_c = c.potential(t);
      const double p_sum = sis->potential(t) + sh->potential(t);
      assert(close(p_c, p_sum, 1.0e-12));
    }
  }

  // ─── 22. SIS + strong shear: 4-image quad config via find_images. ──
  // Shear breaks the SIS axisymmetry → up to 4 images for source
  // inside the tangential caustic. Verify lens-equation closure on
  // every image.
  {
    auto sis = std::make_shared<SIS>(1.0);
    auto sh  = std::make_shared<ExternalShear>(0.10, 0.0);   // axis-aligned
    CompoundLens c{{sis, sh}};
    Vec2 beta{0.05, 0.02};
    auto imgs = find_images(c, beta, /*radius=*/2.0, /*grid=*/120);
    assert(imgs.size() >= 2);                                // at minimum 2
    for (const auto& im : imgs) {
      Vec2 b = lens_equation(c, im.theta);
      assert(close(b.x, beta.x, 1.0e-7, 1.0e-7));
      assert(close(b.y, beta.y, 1.0e-7, 1.0e-7));
      assert(std::isfinite(im.magnification));
    }
  }

  // ─── 23. CompoundLens null component throws. ───────────────────────
  {
    bool threw = false;
    try {
      std::vector<std::shared_ptr<const Lens>> bad{nullptr};
      CompoundLens c{bad};
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try {
      CompoundLens c;
      c.add(nullptr);
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // ─── 24. NFW deflection at x = 1 → α_s · (1 + ln(1/2)). ────────────
  {
    NFW nfw(/*alpha_s=*/2.0, /*theta_s=*/1.0);
    Vec2 a = nfw.deflection({1.0, 0.0});
    const double mag = std::sqrt(a.x * a.x + a.y * a.y);
    const double expected = 2.0 * (1.0 + std::log(0.5));
    assert(close(mag, expected, 1.0e-10));
    assert(std::fabs(a.y) < 1.0e-15);
  }

  // ─── 25. NFW small-x: |α| → 0 (scale-radius-interior cancellation).
  {
    NFW nfw(2.0, 1.0);
    Vec2 a = nfw.deflection({0.001, 0.0});
    const double mag = std::sqrt(a.x * a.x + a.y * a.y);
    // α magnitude should be small (∝ x · ln-like at small x).
    assert(mag < 0.05);
    assert(mag > 0.0);
  }

  // ─── 26. NFW deflection points radially outward (axisymmetric). ───
  {
    NFW nfw(2.0, 1.0);
    for (Vec2 t : {Vec2{1.0, 0.0}, Vec2{0.0, 2.0},
                   Vec2{1.5, 1.5}, Vec2{-0.7, 0.4}}) {
      Vec2 a = nfw.deflection(t);
      // Cross product θ × α = 0 ⇒ collinear with θ.
      assert(std::fabs(t.x * a.y - t.y * a.x) < 1.0e-12);
      // Same sign as θ ⇒ outward not inward.
      assert(t.x * a.x + t.y * a.y > 0.0);
    }
  }

  // ─── 27. NFW translation invariance under off-centre shift. ───────
  {
    Vec2 c{0.4, -0.3};
    NFW nfw_off(2.0, 1.0, c);
    NFW nfw_at0(2.0, 1.0);
    for (Vec2 t_rel : {Vec2{0.7, 0.4}, Vec2{1.5, -0.3}}) {
      Vec2 a_off = nfw_off.deflection({t_rel.x + c.x, t_rel.y + c.y});
      Vec2 a_at0 = nfw_at0.deflection(t_rel);
      assert(close(a_off.x, a_at0.x, 1.0e-12));
      assert(close(a_off.y, a_at0.y, 1.0e-12));
    }
  }

  // ─── 28. NFW + ExternalShear via CompoundLens: find_images
  //     recovers a multi-image config; lens equation closes. ─────────
  {
    auto nfw = std::make_shared<NFW>(2.0, 1.0);
    auto sh  = std::make_shared<ExternalShear>(0.05, 0.0);
    CompoundLens c{{nfw, sh}};
    Vec2 beta{0.05, 0.02};
    auto imgs = find_images(c, beta, /*radius=*/3.0, /*grid=*/120);
    assert(imgs.size() >= 2);                  // multi-image (NFW often has 3)
    for (const auto& im : imgs) {
      Vec2 b = lens_equation(c, im.theta);
      assert(close(b.x, beta.x, 1.0e-7, 1.0e-7));
      assert(close(b.y, beta.y, 1.0e-7, 1.0e-7));
      assert(std::isfinite(im.magnification));
    }
  }

  // ─── 29. NFW bad inputs throw. ──────────────────────────────────────
  {
    bool threw = false;
    try { NFW(0.0, 1.0); } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { NFW(-1.0, 1.0); } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { NFW(1.0, 0.0); } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { NFW(1.0, -0.5); } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // ─── 30. NFW small-x precision: the explicit Wright-Brainerd
  //     formula has catastrophic cancellation at x ≪ 1
  //     (ln(x/2) + arccosh(1/x)/√(1-x²) is the difference of two
  //     ~ln(2/x) terms). The Taylor expansion h(x) ≈ ½ x²(ln(2/x) - ½)
  //     should agree with the kernel at small x to ~ 1e-3 relative. ──
  {
    NFW nfw(/*alpha_s=*/1.0, /*theta_s=*/1.0);
    auto h_via_alpha = [&](double r) {
      // |α| = α_s · h(x) / x  ⇒  h(x) = |α| · x / α_s = |α| · r / θ_s.
      Vec2 a = nfw.deflection({r, 0.0});
      return std::sqrt(a.x * a.x + a.y * a.y) * r;
    };
    auto h_ref = [](double x) {
      return 0.5 * x * x * (std::log(2.0 / x) - 0.5);
    };
    for (double r : {1.0e-4, 5.0e-4, 1.0e-5, 1.0e-7, 1.0e-10}) {
      const double got      = h_via_alpha(r);
      const double expected = h_ref(r);
      assert(close(got, expected, 1.0e-3, 1.0e-15));
      // Direction sanity: α purely along r̂.
      Vec2 a = nfw.deflection({r, 0.0});
      assert(std::fabs(a.y) < 1.0e-15);
    }
    // Continuity around the patch boundary x = 1e-3 (no jump). Use
    // a tiny Δx so the natural derivative-driven shift is negligible
    // — what we're really checking is that the two branches agree
    // at the seam, not that h is constant across a finite window.
    Vec2 a_in  = nfw.deflection({1.0e-3 - 1.0e-9, 0.0});  // taylor branch
    Vec2 a_out = nfw.deflection({1.0e-3 + 1.0e-9, 0.0});  // full branch
    const double m_in  = std::sqrt(a_in.x  * a_in.x  + a_in.y  * a_in.y);
    const double m_out = std::sqrt(a_out.x * a_out.x + a_out.y * a_out.y);
    assert(close(m_in, m_out, 1.0e-5, 1.0e-12));
  }

  return 0;
}
