// Tests for the M5 interferometry forward model.
//
// Validates the closed-form visibility formulas for point sources +
// circular Gaussians + their composition. The substrate claim is that
// the Argus IR + Retrieval pattern works on interferometric
// visibilities the same way it works on JWST spectra and lensed-image
// positions; this file proves the physics layer is correct first.

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "argus/argus.hpp"

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

}  // namespace

int main() {
  using namespace argus::interferometry;

  // ─── 1. Point source at origin: V(u, v) = F (no phase). ────────────
  {
    PointSource src{0.0, 0.0, /*flux=*/2.5};
    for (UVPoint uv : {UVPoint{0.0, 0.0}, UVPoint{100.0, 0.0},
                       UVPoint{-50.0, 80.0}, UVPoint{1.0e4, 1.0e4}}) {
      Visibility v = predict_visibility(src, uv);
      assert(close(v.real, 2.5, 1.0e-12));
      assert(std::fabs(v.imag) < 1.0e-12);
    }
  }

  // ─── 2. Off-origin point: |V|=F, phase = -2π(u·l + v·m). ───────────
  {
    const double l = 1.0e-5, m = -3.0e-6, F = 1.0;
    PointSource src{l, m, F};
    UVPoint uv{1.0e5, 4.0e4};
    Visibility v = predict_visibility(src, uv);
    const double mag2 = v.real * v.real + v.imag * v.imag;
    assert(close(std::sqrt(mag2), F, 1.0e-12));
    const double phase = -kTwoPi * (uv.u * l + uv.v * m);
    assert(close(v.real, F * std::cos(phase), 1.0e-12));
    assert(close(v.imag, F * std::sin(phase), 1.0e-12));
  }

  // ─── 3. Visibilities add over multiple components. ─────────────────
  {
    std::vector<PointSource> srcs{
      {1.0e-5, 0.0, 1.0},
      {0.0, 1.0e-5, 0.5},
      {-2.0e-5, 1.5e-5, 0.3},
    };
    std::vector<UVPoint> uvs{{1.0e5, 0.0}, {0.0, 1.0e5}, {5.0e4, 5.0e4}};
    auto vs = predict_visibilities(srcs, uvs);
    for (std::size_t i = 0; i < uvs.size(); ++i) {
      Visibility expected{0.0, 0.0};
      for (const auto& s : srcs) {
        Visibility v = predict_visibility(s, uvs[i]);
        expected.real += v.real;
        expected.imag += v.imag;
      }
      assert(close(vs[i].real, expected.real, 1.0e-12));
      assert(close(vs[i].imag, expected.imag, 1.0e-12));
    }
  }

  // ─── 4. Translation theorem: shifting the source by Δ adds a
  //     phase -2π(u·Δ_l + v·Δ_m) — i.e. multiplies V by exp(-2πi(...))
  {
    PointSource s0{0.0, 0.0, 1.0};
    PointSource s1{1.0e-5, -3.0e-6, 1.0};   // shift by (Δ_l, Δ_m)
    UVPoint uv{2.5e5, -1.0e5};
    Visibility v0 = predict_visibility(s0, uv);
    Visibility v1 = predict_visibility(s1, uv);
    const double phase = -kTwoPi * (uv.u * s1.l + uv.v * s1.m);
    // v1 = v0 · exp(i·phase)
    const double exp_re = v0.real * std::cos(phase) - v0.imag * std::sin(phase);
    const double exp_im = v0.real * std::sin(phase) + v0.imag * std::cos(phase);
    assert(close(v1.real, exp_re, 1.0e-12));
    assert(close(v1.imag, exp_im, 1.0e-12));
  }

  // ─── 5. Conjugate symmetry: for a real sky, V(-u,-v) = conj(V(u,v))
  {
    std::vector<PointSource> srcs{
      {1.0e-5, 0.0, 0.7},
      {-2.0e-5, 1.5e-5, 0.3},
    };
    std::vector<UVPoint> uvs{{2.0e5, 1.0e5}, {-2.0e5, -1.0e5}};
    auto vs = predict_visibilities(srcs, uvs);
    assert(close(vs[0].real,  vs[1].real, 1.0e-12));
    assert(close(vs[0].imag, -vs[1].imag, 1.0e-12));
  }

  // ─── 6. Gaussian: V at origin = F, decays as exp(-2π² σ² r²). ──────
  {
    GaussianSource src{0.0, 0.0, /*F=*/1.0, /*sigma=*/1.0e-5};
    Visibility v0 = predict_visibility(src, {0.0, 0.0});
    assert(close(v0.real, 1.0, 1.0e-12));
    assert(std::fabs(v0.imag) < 1.0e-12);

    // At a UV separation r:
    //   |V| = F · exp(-2π² σ² r²)
    for (double r : {1.0e3, 5.0e4, 2.0e5}) {
      Visibility v = predict_visibility(src, {r, 0.0});
      const double expected = std::exp(-2.0 * M_PI * M_PI *
                                       src.sigma * src.sigma * r * r);
      assert(close(v.real, expected, 1.0e-12));
      assert(std::fabs(v.imag) < 1.0e-12);
    }
  }

  // ─── 7. Gaussian σ=0 reduces to a point source. ─────────────────────
  {
    GaussianSource g{1.0e-5, -3.0e-6, 0.7, /*sigma=*/0.0};
    PointSource    p{1.0e-5, -3.0e-6, 0.7};
    for (UVPoint uv : {UVPoint{1.0e5, 0.0}, UVPoint{2.0e4, 8.0e4}}) {
      Visibility vg = predict_visibility(g, uv);
      Visibility vp = predict_visibility(p, uv);
      assert(close(vg.real, vp.real, 1.0e-12));
      assert(close(vg.imag, vp.imag, 1.0e-12));
    }
  }

  // ─── 8. Mixed Gaussian + point — visibilities sum component-wise. ──
  {
    std::vector<GaussianSource> srcs{
      {1.0e-5, 0.0, 1.0, 1.0e-5},
      {0.0,    0.0, 0.5, 0.0},          // point via σ=0
    };
    std::vector<UVPoint> uvs{{1.0e4, 0.0}, {1.0e5, 1.0e5}};
    auto vs = predict_visibilities(srcs, uvs);
    for (std::size_t i = 0; i < uvs.size(); ++i) {
      Visibility expected{0.0, 0.0};
      for (const auto& s : srcs) {
        Visibility v = predict_visibility(s, uvs[i]);
        expected.real += v.real;
        expected.imag += v.imag;
      }
      assert(close(vs[i].real, expected.real, 1.0e-12, 1.0e-15));
      assert(close(vs[i].imag, expected.imag, 1.0e-12, 1.0e-15));
    }
  }

  // ─── 9. UV coverage from a tiny array: 3 antennas in a triangle. ──
  // Antennas at (0,0), (1000, 0), (500, 866) m at λ=0.21 m (HI).
  // Expected: 3 baselines = (1000, 0), (500, 866), (-500, 866) in
  // metres → divide by λ for UV (in wavelengths).
  {
    std::vector<double> east{0.0, 1000.0, 500.0};
    std::vector<double> north{0.0, 0.0, 866.0};
    auto uvs = uv_coverage_snapshot(east, north, /*wavelength=*/0.21);
    assert(uvs.size() == 3);
    // Pair (0,1): bx=1000, by=0
    assert(close(uvs[0].u, 1000.0 / 0.21, 1.0e-10));
    assert(close(uvs[0].v, 0.0,            1.0e-10));
    // Pair (0,2): bx=500, by=866
    assert(close(uvs[1].u, 500.0 / 0.21, 1.0e-10));
    assert(close(uvs[1].v, 866.0 / 0.21, 1.0e-10));
    // Pair (1,2): bx=-500, by=866
    assert(close(uvs[2].u, -500.0 / 0.21, 1.0e-10));
    assert(close(uvs[2].v,  866.0 / 0.21, 1.0e-10));
  }

  // ─── 10. Bad inputs throw. ───────────────────────────────────────
  {
    bool threw = false;
    try { (void)uv_coverage_snapshot({0.0, 1.0}, {0.0}, 1.0); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { (void)uv_coverage_snapshot({0.0}, {0.0}, 1.0); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try { (void)uv_coverage_snapshot({0.0, 1.0}, {0.0, 1.0}, 0.0); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    threw = false;
    try {
      std::vector<GaussianSource> bad{{0, 0, 1, /*sigma=*/-1.0}};
      (void)predict_visibilities(bad, std::vector<UVPoint>{{0.0, 0.0}});
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
  }

  // ─── 11. Determinism: same inputs → bit-equal outputs. ────────────
  {
    std::vector<PointSource> srcs{
      {1.0e-5, 2.0e-5, 1.5},
      {-3.0e-5, 1.0e-5, 0.7},
    };
    std::vector<UVPoint> uvs{{1.0e5, 0.0}, {0.0, 5.0e4}, {2.0e5, 2.0e5}};
    auto a = predict_visibilities(srcs, uvs);
    auto b = predict_visibilities(srcs, uvs);
    assert(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
      assert(a[i].real == b[i].real);
      assert(a[i].imag == b[i].imag);
    }
  }

  return 0;
}
