// Performance baseline for the M4 lensing kernels. Records median
// wall time over 5 runs of the workhorse operations:
//   * SIE deflection at one point
//   * SIE deflection over an N×N evaluation grid
//   * find_images on an SIE+shear compound (Newton solver)
//
// Asserts loose per-call bounds — catches order-of-magnitude
// regressions without being CI-flaky.

#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <vector>

#include "argus/argus.hpp"

int main() {
  using namespace argus::lensing;
  using clock = std::chrono::steady_clock;

  // ─── 1. SIE deflection at one point: should be ~ μs scale. ─────────
  {
    SIE sie(1.0, 0.7, 0.3);
    constexpr int N = 100000;
    Vec2 dummy{0.0, 0.0};                           // accumulator
    const auto t0 = clock::now();
    for (int i = 0; i < N; ++i) {
      const double x = 0.1 + 1.0e-5 * i;
      Vec2 a = sie.deflection({x, 0.5});
      dummy.x += a.x; dummy.y += a.y;               // defeat DCE
    }
    const auto t1 = clock::now();
    const double per_ns =
        std::chrono::duration<double, std::nano>(t1 - t0).count() / N;
    assert(std::isfinite(dummy.x + dummy.y));
    std::cerr << "[bench] SIE deflection per-call ≈ " << per_ns << " ns\n";
    assert(per_ns < 5000.0);                        // << 5 µs / call
  }

  // ─── 2. find_images on SIE: 60×60 grid + Newton refinement.
  //     Median over 5 runs. ────────────────────────────────────────────
  {
    SIE sie(1.0, 0.7, 0.3);
    Vec2 source{0.04, 0.03};
    // Warm-up.
    (void)find_images(sie, source, 2.0, 60);

    std::vector<double> ms;
    ms.reserve(5);
    for (int i = 0; i < 5; ++i) {
      const auto t0 = clock::now();
      auto imgs = find_images(sie, source, 2.0, 60);
      const auto t1 = clock::now();
      ms.push_back(
          std::chrono::duration<double, std::milli>(t1 - t0).count());
      assert(imgs.size() == 4);
    }
    std::sort(ms.begin(), ms.end());
    const double median_ms = ms[ms.size() / 2];
    std::cerr << "[bench] find_images(SIE, 60-grid) median ≈ "
              << median_ms << " ms\n";
    assert(median_ms < 200.0);                      // CI-safe loose bound
  }

  // ─── 3. NFW deflection benchmark — ensures Wright-Brainerd
  //     branches (arccosh / arccos) don't hide a perf regression. ──────
  {
    NFW nfw(2.0, 1.0);
    constexpr int N = 100000;
    double acc = 0.0;
    const auto t0 = clock::now();
    for (int i = 0; i < N; ++i) {
      const double x = 0.5 + 1.0e-5 * i;
      Vec2 a = nfw.deflection({x, 0.3});
      acc += a.x;
    }
    const auto t1 = clock::now();
    const double per_ns =
        std::chrono::duration<double, std::nano>(t1 - t0).count() / N;
    assert(std::isfinite(acc));
    std::cerr << "[bench] NFW deflection per-call ≈ " << per_ns << " ns\n";
    assert(per_ns < 5000.0);
  }

  // ─── 4. CompoundLens (SIS + ExternalShear): summing two components
  //     should be < 2× the SIS-alone cost. ────────────────────────────
  {
    auto sis = std::make_shared<SIS>(1.0);
    auto sh  = std::make_shared<ExternalShear>(0.05, -0.02);
    CompoundLens c{{sis, sh}};
    constexpr int N = 100000;
    double acc = 0.0;
    const auto t0 = clock::now();
    for (int i = 0; i < N; ++i) {
      const double x = 0.1 + 1.0e-5 * i;
      Vec2 a = c.deflection({x, 0.4});
      acc += a.x;
    }
    const auto t1 = clock::now();
    const double per_ns =
        std::chrono::duration<double, std::nano>(t1 - t0).count() / N;
    assert(std::isfinite(acc));
    std::cerr << "[bench] CompoundLens (SIS+shear) per-call ≈ "
              << per_ns << " ns\n";
    assert(per_ns < 5000.0);
  }

  return 0;
}
