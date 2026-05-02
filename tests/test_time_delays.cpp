// Tests for the M4 time-delay surface (Fermat potential).
//
// Time delays between multiply-imaged quasar images are the
// observable that drives H0 measurement via the TDCOSMO/H0LiCOW
// methodology. Argus's substrate claim covers this end-to-end:
// once the gravitational potential ψ(θ) is defined, the Fermat
// potential τ = 0.5|θ-β|² - ψ(θ) and the time-delay difference
// between images follow mechanically.
//
// This file verifies:
//   * SIS potential matches the closed form ψ_SIS = θ_E · |θ - centre|.
//   * SIE potential is gradient-consistent with the deflection (so
//     ∇ψ_SIE = α_SIE to ~1e-6 via finite differences).
//   * SIE potential at q=1 reduces to SIS to ~1e-12.
//   * Fermat potential is critical at images (∇τ = 0 at θ that
//     solves the lens equation).
//   * Closed-form SIS time delay Δτ = 2 · θ_E · |β - centre|
//     between the two on-axis images of an aligned source.
//   * Off-axis source: the Argus time-delay machinery agrees with a
//     direct evaluation of (τ_2 - τ_1).
//   * Translation invariance under lens-centre shift.

#include <cassert>
#include <cmath>
#include <vector>

#include "argus/argus.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

}  // namespace

int main() {
  using namespace argus::lensing;

  // ─── 1. SIS potential = θ_E · |θ - centre| (closed form). ──────────
  {
    SIS sis(1.5, /*centre=*/{0.3, -0.2});
    for (Vec2 t : {Vec2{1.0, 0.5}, Vec2{-0.4, 1.7}, Vec2{0.3, -0.2}}) {
      const double dx = t.x - 0.3;
      const double dy = t.y + 0.2;
      const double exp = 1.5 * std::sqrt(dx * dx + dy * dy);
      assert(close(sis.potential(t), exp, 1.0e-12));
    }
  }

  // ─── 2. ∇ψ_SIE = α_SIE — gradient consistency via finite diffs. ────
  {
    SIE sie(/*theta_E=*/1.0, /*q=*/0.7, /*phi=*/0.5);
    const double h = 1.0e-5;
    for (Vec2 t : {Vec2{1.2, 0.4}, Vec2{-0.6, 0.9}, Vec2{0.3, -1.1}}) {
      const double psi_xp = sie.potential({t.x + h, t.y});
      const double psi_xm = sie.potential({t.x - h, t.y});
      const double psi_yp = sie.potential({t.x, t.y + h});
      const double psi_ym = sie.potential({t.x, t.y - h});
      const double gx = (psi_xp - psi_xm) / (2.0 * h);
      const double gy = (psi_yp - psi_ym) / (2.0 * h);
      Vec2 alpha = sie.deflection(t);
      // FD with h=1e-5 on an arctan/arctanh field gets ~1e-7 absolute.
      assert(close(gx, alpha.x, 1.0e-5, 1.0e-7));
      assert(close(gy, alpha.y, 1.0e-5, 1.0e-7));
    }
  }

  // ─── 3. SIE q=1 reduces to SIS for the potential. ──────────────────
  {
    SIS sis(1.2);
    SIE sie(1.2, /*q=*/1.0);
    for (Vec2 t : {Vec2{1.0, 0.5}, Vec2{-0.3, 2.0}, Vec2{0.5, -1.2}}) {
      assert(close(sis.potential(t), sie.potential(t), 1.0e-12, 1.0e-12));
    }
  }

  // ─── 4. Fermat potential is critical at lens-equation roots. ───────
  // The condition ∇_θ τ = 0 is exactly the lens equation θ - α(θ) = β.
  // Verify directly: for an SIS-with-source-inside the two images
  // satisfy ∇τ = 0 to ≈ 1e-7 (FD precision).
  {
    SIS sis(1.0);
    Vec2 beta{0.4, 0.0};
    auto imgs = sis_images(sis, beta);
    assert(imgs.size() == 2);
    const double h = 1.0e-5;
    for (const auto& im : imgs) {
      const double tau_xp = fermat_potential(sis, {im.theta.x + h, im.theta.y}, beta);
      const double tau_xm = fermat_potential(sis, {im.theta.x - h, im.theta.y}, beta);
      const double tau_yp = fermat_potential(sis, {im.theta.x, im.theta.y + h}, beta);
      const double tau_ym = fermat_potential(sis, {im.theta.x, im.theta.y - h}, beta);
      const double gx = (tau_xp - tau_xm) / (2.0 * h);
      const double gy = (tau_yp - tau_ym) / (2.0 * h);
      assert(std::fabs(gx) < 1.0e-7);
      assert(std::fabs(gy) < 1.0e-7);
    }
  }

  // ─── 5. SIS closed-form time delay Δτ = 2 · θ_E · β (on axis). ─────
  // For a source at β > 0 on the x-axis through an SIS centred at the
  // origin, the two images are at θ_+ = β + θ_E and θ_- = β - θ_E.
  // The Fermat-potential difference Δτ = τ(θ_-) - τ(θ_+) = 2 θ_E β.
  {
    const double te = 1.0;
    const double b  = 0.3;
    SIS sis(te);
    Vec2 beta{b, 0.0};
    Vec2 theta_plus{b + te, 0.0};
    Vec2 theta_minus{b - te, 0.0};
    const double dt = time_delay_arcsec2(sis, theta_plus, theta_minus, beta);
    const double exp = 2.0 * te * b;
    assert(close(dt, exp, 1.0e-12));
  }

  // ─── 6. Off-axis source: numerical agreement with direct τ diff. ───
  {
    SIS sis(1.0);
    Vec2 beta{0.25, 0.18};
    auto imgs = sis_images(sis, beta);
    assert(imgs.size() == 2);
    const double tau_a = fermat_potential(sis, imgs[0].theta, beta);
    const double tau_b = fermat_potential(sis, imgs[1].theta, beta);
    const double dt    = time_delay_arcsec2(sis, imgs[0].theta, imgs[1].theta, beta);
    assert(close(dt, tau_b - tau_a, 1.0e-12));
    // Sign: time-delay depends on which image is "image_a" vs "image_b";
    // swapping flips the sign.
    const double dt_swap = time_delay_arcsec2(sis, imgs[1].theta, imgs[0].theta, beta);
    assert(close(dt_swap, -dt, 1.0e-12));
  }

  // ─── 7. Translation invariance: shifting lens, source, and images
  //     by the same vector leaves Δτ unchanged. ──────────────────────
  {
    Vec2 shift{1.7, -0.8};
    SIS sis_at0(1.0);
    SIS sis_off(1.0, /*centre=*/shift);
    Vec2 beta_at0{0.3, 0.1};
    Vec2 beta_off{beta_at0.x + shift.x, beta_at0.y + shift.y};
    auto imgs = sis_images(sis_at0, beta_at0);
    Vec2 t1_at0 = imgs[0].theta;
    Vec2 t2_at0 = imgs[1].theta;
    Vec2 t1_off{t1_at0.x + shift.x, t1_at0.y + shift.y};
    Vec2 t2_off{t2_at0.x + shift.x, t2_at0.y + shift.y};
    const double dt_at0 = time_delay_arcsec2(sis_at0, t1_at0, t2_at0, beta_at0);
    const double dt_off = time_delay_arcsec2(sis_off, t1_off, t2_off, beta_off);
    assert(close(dt_at0, dt_off, 1.0e-12));
  }

  // ─── 8. SIE: time delays between the 4 images of a cusp source. ────
  // Each pairwise delay must equal the difference of Fermat potentials.
  // The 4 delays should not all be zero (no degeneracy here).
  {
    SIE sie(1.0, /*q=*/0.7);
    Vec2 beta{0.05, 0.03};
    auto imgs = find_images(sie, beta, /*radius=*/2.0, /*grid=*/100);
    assert(imgs.size() == 4);
    std::vector<double> taus;
    for (const auto& im : imgs) {
      taus.push_back(fermat_potential(sie, im.theta, beta));
    }
    // Pairwise consistency: time_delay must match τ-diff to 1e-12.
    for (std::size_t i = 0; i < imgs.size(); ++i) {
      for (std::size_t j = 0; j < imgs.size(); ++j) {
        const double dt = time_delay_arcsec2(sie, imgs[i].theta, imgs[j].theta, beta);
        assert(close(dt, taus[j] - taus[i], 1.0e-12, 1.0e-12));
      }
    }
    // Not all delays are equal (otherwise the configuration would be
    // degenerate and time-delay cosmography wouldn't work).
    bool any_nonzero = false;
    for (std::size_t i = 0; i < imgs.size(); ++i) {
      for (std::size_t j = i + 1; j < imgs.size(); ++j) {
        if (std::fabs(taus[i] - taus[j]) > 1.0e-6) any_nonzero = true;
      }
    }
    assert(any_nonzero);
  }

  return 0;
}
