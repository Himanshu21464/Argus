// Performance baseline for the M5 interferometry kernels. Records
// median wall time over 5 runs of:
//   * predict_visibility for one PointSource at one UV point
//   * predict_visibilities summed over multi-Gaussian / Earth-rotation track
//   * uv_coverage_track for a 27-antenna VLA-sized array, 30 HAs
//
// Asserts loose per-call bounds — catches order-of-magnitude
// regressions without being CI-flaky.

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include "argus/argus.hpp"

int main() {
  using namespace argus::interferometry;
  using clock = std::chrono::steady_clock;

  // ─── 1. predict_visibility per-call (PointSource): sub-µs target. ──
  {
    PointSource src{1.0e-5, -3.0e-6, 1.0};
    constexpr int N = 1000000;
    double acc = 0.0;
    const auto t0 = clock::now();
    for (int i = 0; i < N; ++i) {
      Visibility v = predict_visibility(src, {1.0e3 + i, 5.0e2 - i});
      acc += v.real;
    }
    const auto t1 = clock::now();
    const double per_ns =
        std::chrono::duration<double, std::nano>(t1 - t0).count() / N;
    assert(std::isfinite(acc));
    std::cerr << "[bench] predict_visibility(point) per-call ≈ "
              << per_ns << " ns\n";
    assert(per_ns < 1000.0);                        // < 1 µs
  }

  // ─── 2. predict_visibilities for 2-Gaussian × VLA-track UV. ────────
  // 27-antenna VLA D-config, 30 HAs spanning 5 hours — 350 baselines × 30
  // = 10500 UV points × 2 Gaussian = 21000 visibility evaluations.
  {
    std::vector<double> east(27), north(27);
    for (std::size_t i = 0; i < 27; ++i) {
      const double t = 2.0 * M_PI * static_cast<double>(i) / 27.0;
      // Spiral arms (toy VLA).
      east [i] = 500.0 * std::cos(t) * (1.0 + static_cast<double>(i) / 27.0);
      north[i] = 500.0 * std::sin(t) * (1.0 + static_cast<double>(i) / 27.0);
    }
    std::vector<double> ha;
    for (int k = 0; k < 30; ++k) ha.push_back(-2.5 + k * 5.0 / 29.0);
    auto uv = uv_coverage_track(east, north, /*lat=*/0.6, ha,
                                /*dec=*/0.5, /*λ=*/0.21);
    assert(uv.size() == 27 * 26 / 2 * 30);

    std::vector<GaussianSource> srcs{
      {0.0,    0.0,    1.0, 1.0e-5},
      {1.5e-5, 5.0e-6, 0.3, 3.0e-5},
    };
    // Warm-up.
    (void)predict_visibilities(srcs, uv);

    std::vector<double> ms;
    ms.reserve(5);
    for (int i = 0; i < 5; ++i) {
      const auto t0 = clock::now();
      auto vs = predict_visibilities(srcs, uv);
      const auto t1 = clock::now();
      ms.push_back(
          std::chrono::duration<double, std::milli>(t1 - t0).count());
      assert(vs.size() == uv.size());
    }
    std::sort(ms.begin(), ms.end());
    const double median_ms = ms[ms.size() / 2];
    const double per_eval_ns =
        median_ms * 1.0e6 /
        (static_cast<double>(uv.size()) * static_cast<double>(srcs.size()));
    std::cerr << "[bench] predict_visibilities(2-Gaussian, "
              << uv.size() << " UV) median = " << median_ms
              << " ms (" << per_eval_ns << " ns/eval)\n";
    assert(median_ms < 500.0);                      // CI-safe bound
    assert(per_eval_ns < 200.0);                    // sin/cos+exp ≈ 30-100 ns
  }

  // ─── 3. uv_coverage_track scaling: 27 antennas × 30 HAs. ──────────
  {
    std::vector<double> east(27), north(27);
    for (std::size_t i = 0; i < 27; ++i) {
      east [i] = static_cast<double>(i) * 100.0;
      north[i] = static_cast<double>(i) * 50.0;
    }
    std::vector<double> ha;
    for (int k = 0; k < 30; ++k) ha.push_back(0.1 * k);

    // Warm-up.
    (void)uv_coverage_track(east, north, 0.6, ha, 0.5, 0.21);

    std::vector<double> ms;
    for (int i = 0; i < 5; ++i) {
      const auto t0 = clock::now();
      auto uv = uv_coverage_track(east, north, 0.6, ha, 0.5, 0.21);
      const auto t1 = clock::now();
      ms.push_back(
          std::chrono::duration<double, std::milli>(t1 - t0).count());
      assert(uv.size() == 27 * 26 / 2 * 30);
    }
    std::sort(ms.begin(), ms.end());
    const double median_ms = ms[ms.size() / 2];
    std::cerr << "[bench] uv_coverage_track(27 ant, 30 HA) median ≈ "
              << median_ms << " ms\n";
    assert(median_ms < 100.0);
  }

  return 0;
}
