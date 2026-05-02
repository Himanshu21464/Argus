// End-to-end M5 demo: a compact-core + extended-jet 2-Gaussian radio
// source observed by a 7-antenna VLA-like Y-array over a 5-hour
// Earth-rotation track at 21 cm (HI line). Argus computes:
//   * the UV coverage (105 baselines × 2 polarisation-like real/imag
//     components)
//   * the predicted complex visibility at every UV sample
//   * a brief ASCII visibility magnitude summary
//
// This is the same kernel surface that the test_multi_component_retrieval
// substrate proof MCMC-fits to recover the 8 source parameters from
// noisy synthetic visibilities.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "argus/argus.hpp"

int main() {
  using namespace argus::interferometry;

  // ─── 7-antenna VLA-like Y-array (positions in metres). ─────────────
  std::vector<double> east { 0.0, -250.0,  250.0,    0.0, -500.0,  500.0,    0.0};
  std::vector<double> north{ 0.0, -433.0, -433.0,  500.0, -866.0, -866.0, 1000.0};

  // ─── 5-hour Earth-rotation track at dec=30°, lat=34°. ──────────────
  std::vector<double> hour_angles;
  for (int k = -2; k <= 2; ++k) hour_angles.push_back(k * 0.5);  // 5 HAs
  const double declination = 30.0 * M_PI / 180.0;
  const double latitude    = 34.0 * M_PI / 180.0;
  const double wavelength  = 0.21;                                // HI 21 cm

  auto uv = uv_coverage_track(east, north, latitude, hour_angles,
                              declination, wavelength);
  printf("=== M5 interferometry demo ===\n\n");
  printf("Array: 7-antenna VLA-like Y configuration\n");
  printf("Track: 5 HAs spanning ±2.5 h, source dec=30°, station lat=34°\n");
  printf("λ:     %.3f m (HI 21 cm)\n", wavelength);
  printf("UV coverage: %zu samples (21 baselines × 5 HAs)\n\n", uv.size());

  // Find max baseline length in wavelengths.
  double max_uv = 0.0;
  for (const auto& p : uv) {
    const double r = std::sqrt(p.u * p.u + p.v * p.v);
    if (r > max_uv) max_uv = r;
  }
  printf("max |uv| = %.0f wavelengths → angular resolution ≈ %.2f arcsec\n\n",
         max_uv, 1.0 / max_uv * 206265.0);

  // ─── Sky model: compact core + extended jet. ───────────────────────
  std::vector<GaussianSource> sources{
    /*core=*/ {0.0,    0.0,    1.0, 1.0e-5},
    /*jet=*/  {1.5e-5, 5.0e-6, 0.3, 3.0e-5},
  };
  printf("Sky model:\n");
  printf("  core  : (l, m) = (%+8.2e, %+8.2e), F=%.2f Jy, σ=%8.2e rad\n",
         sources[0].l, sources[0].m, sources[0].flux, sources[0].sigma);
  printf("  jet   : (l, m) = (%+8.2e, %+8.2e), F=%.2f Jy, σ=%8.2e rad\n",
         sources[1].l, sources[1].m, sources[1].flux, sources[1].sigma);
  printf("\n");

  auto vis = predict_visibilities(sources, uv);

  // ─── Summary: |V| histogram, sorted by baseline length. ────────────
  // Bin baselines into 8 |uv| bands and report mean |V| in each.
  const int n_bins = 8;
  std::vector<double> sum_v(n_bins, 0.0);
  std::vector<int>    cnt  (n_bins, 0);
  const double dr = max_uv / static_cast<double>(n_bins);
  for (std::size_t i = 0; i < uv.size(); ++i) {
    const double r = std::sqrt(uv[i].u * uv[i].u + uv[i].v * uv[i].v);
    int bin = static_cast<int>(r / dr);
    if (bin >= n_bins) bin = n_bins - 1;
    const double mag = std::sqrt(vis[i].real * vis[i].real +
                                 vis[i].imag * vis[i].imag);
    sum_v[static_cast<std::size_t>(bin)] += mag;
    cnt  [static_cast<std::size_t>(bin)] += 1;
  }
  printf("Visibility |V| vs baseline length (bin = %.0f wavelengths):\n", dr);
  printf("  bin  |uv| range (λ)           N    <|V|>   bar\n");
  for (int b = 0; b < n_bins; ++b) {
    if (cnt[static_cast<std::size_t>(b)] == 0) continue;
    const double mean = sum_v[static_cast<std::size_t>(b)] /
                        cnt[static_cast<std::size_t>(b)];
    const int bar_len = static_cast<int>(40.0 * mean / sources[0].flux);
    printf("  %2d   %6.0f .. %-6.0f   %3d   %.3f   ",
           b, b * dr, (b + 1) * dr, cnt[static_cast<std::size_t>(b)], mean);
    for (int k = 0; k < bar_len; ++k) printf("█");
    printf("\n");
  }

  printf("\nThe core (compact) contributes a flat-amplitude visibility ≈ 1 Jy\n");
  printf("everywhere; the jet (extended) contributes only at short baselines\n");
  printf("(< 1/(2π σ_jet) ≈ %.0f λ). The decay of |V| with baseline length\n",
         1.0 / (2.0 * M_PI * sources[1].sigma));
  printf("encodes the extended structure — exactly what the substrate-proof\n");
  printf("retrieval inverts to recover (l, m, F, σ) for both components.\n");

  return 0;
}
