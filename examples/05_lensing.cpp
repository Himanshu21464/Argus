// End-to-end M4 demo: an SIE galaxy with external shear lensing a
// background quasar produces 4 quad images. Argus computes:
//   * the 4 image positions via the generic Newton solver
//   * each image's magnification
//   * the 6 pairwise Fermat-potential time-delay differences
//
// This is the same kernel surface that the test_sie_retrieval substrate
// proof MCMC-fits to recover lens parameters from observed image
// positions.

#include <cstdio>
#include <memory>
#include <vector>

#include "argus/argus.hpp"

int main() {
  using namespace argus::lensing;

  // ─── Lens model: SIE + external shear (the workhorse used in
  //     real H0LiCOW/TDCOSMO analyses). ───────────────────────────────
  auto sie = std::make_shared<SIE>(/*theta_E=*/1.0,
                                   /*q=*/0.7,
                                   /*phi=*/0.3,
                                   /*centre=*/Vec2{0.05, -0.02});
  auto shear = std::make_shared<ExternalShear>(/*g1=*/0.05, /*g2=*/-0.02);
  CompoundLens lens{{sie, shear}};

  // ─── Source position (inside the tangential caustic for a quad). ──
  Vec2 source{0.04, 0.03};

  printf("=== M4 lensing demo ===\n\n");
  printf("Lens: SIE (θ_E=1.0\", q=0.7, φ=0.3 rad, centre=(0.05, -0.02)\")\n");
  printf("    + ExternalShear (γ₁=0.05, γ₂=-0.02)\n");
  printf("Source: (%+.3f, %+.3f) arcsec\n\n", source.x, source.y);

  // ─── Image-finding via the generic Newton solver. ──────────────────
  auto imgs = find_images(lens, source, /*radius=*/2.0, /*grid=*/120);
  printf("found %zu images:\n", imgs.size());
  printf("  %2s   %10s   %10s   %12s   %14s\n",
         "ID", "θ_x", "θ_y", "magnification", "fermat τ");
  std::vector<double> taus;
  for (std::size_t i = 0; i < imgs.size(); ++i) {
    const auto& im = imgs[i];
    const double tau = fermat_potential(lens, im.theta, source);
    taus.push_back(tau);
    printf("  %2zu   %+10.4f   %+10.4f   %12.3f   %+14.6f\n",
           i + 1, im.theta.x, im.theta.y, im.magnification, tau);
  }

  // ─── Pairwise time delays Δτ in arcsec². ──────────────────────────
  printf("\npairwise time-delay differences Δτ (arcsec²):\n");
  printf("  ");
  for (std::size_t j = 0; j < imgs.size(); ++j) printf("    img %zu  ", j + 1);
  printf("\n");
  for (std::size_t i = 0; i < imgs.size(); ++i) {
    printf("  img %zu  ", i + 1);
    for (std::size_t j = 0; j < imgs.size(); ++j) {
      printf("%+10.6f  ", taus[j] - taus[i]);
    }
    printf("\n");
  }

  printf("\n");
  printf("Δτ × (1+z_l)·D_Δt/c gives the physical time delay in days.\n");
  printf("This is the observable measured by COSMOGRAIL / TDCOSMO and\n");
  printf("inverted to constrain H_0.\n");

  return 0;
}
