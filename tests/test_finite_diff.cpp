// Validate forward-mode autograd (Dual<double>) end-to-end against
// central finite differences. Verifies the dual chain rule survives
// through the entire Voigt + intensity-weighted line shape.

#include <cassert>
#include <cmath>

#include "argus/argus.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

// Templated single-line cross-section evaluator. With T = double we
// recover the standard floating-point answer; with T = Dual<double>
// we get derivatives wrt whichever parameter is seeded with d=1.
template <typename T>
T single_line_xsec(T x, T sigma_g, T gamma_l, T S) {
  return S * argus::voigt(x, sigma_g, gamma_l);
}

}  // namespace

int main() {
  using namespace argus;
  using D = Dual<double>;

  // Test point.
  const double x_v     = 0.012;
  const double sigma_v = 0.04;
  const double gamma_v = 0.06;
  const double S_v     = 1.0e-20;
  const double h       = 1.0e-7;

  // 1. d/d(sigma)
  {
    D out = single_line_xsec(D{x_v, 0.0}, D{sigma_v, 1.0},
                             D{gamma_v, 0.0}, D{S_v, 0.0});
    const double dp = single_line_xsec<double>(x_v, sigma_v + h, gamma_v, S_v);
    const double dm = single_line_xsec<double>(x_v, sigma_v - h, gamma_v, S_v);
    const double fd = (dp - dm) / (2.0 * h);
    assert(close(out.d, fd, 1.0e-5, 1.0e-30));
  }

  // 2. d/d(gamma)
  {
    D out = single_line_xsec(D{x_v, 0.0}, D{sigma_v, 0.0},
                             D{gamma_v, 1.0}, D{S_v, 0.0});
    const double dp = single_line_xsec<double>(x_v, sigma_v, gamma_v + h, S_v);
    const double dm = single_line_xsec<double>(x_v, sigma_v, gamma_v - h, S_v);
    const double fd = (dp - dm) / (2.0 * h);
    assert(close(out.d, fd, 1.0e-5, 1.0e-30));
  }

  // 3. d/d(x)
  {
    D out = single_line_xsec(D{x_v, 1.0}, D{sigma_v, 0.0},
                             D{gamma_v, 0.0}, D{S_v, 0.0});
    const double dp = single_line_xsec<double>(x_v + h, sigma_v, gamma_v, S_v);
    const double dm = single_line_xsec<double>(x_v - h, sigma_v, gamma_v, S_v);
    const double fd = (dp - dm) / (2.0 * h);
    assert(close(out.d, fd, 1.0e-5, 1.0e-30));
  }

  // 4. d/d(S) — scalar linear, should be exactly voigt(x, sigma, gamma).
  {
    D out = single_line_xsec(D{x_v, 0.0}, D{sigma_v, 0.0},
                             D{gamma_v, 0.0}, D{S_v, 1.0});
    const double phi = voigt<double>(x_v, sigma_v, gamma_v);
    assert(close(out.d, phi, 1.0e-12));
  }

  // 5. Sweep across the line profile and assert FD agreement at every
  //    sample. This catches any numeric instability in the polynomial
  //    evaluation away from the line centre.
  {
    for (double x = -0.4; x <= 0.4; x += 0.05) {
      D out = single_line_xsec(D{x, 0.0}, D{sigma_v, 0.0},
                               D{gamma_v, 1.0}, D{S_v, 0.0});
      const double dp = single_line_xsec<double>(x, sigma_v, gamma_v + h, S_v);
      const double dm = single_line_xsec<double>(x, sigma_v, gamma_v - h, S_v);
      const double fd = (dp - dm) / (2.0 * h);
      assert(close(out.d, fd, 1.0e-4, 1.0e-30));
    }
  }

  return 0;
}
