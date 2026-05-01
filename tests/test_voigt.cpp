#include <cassert>
#include <cmath>

#include "argus/argus.hpp"

int main() {
  using namespace argus;

  // 1. Voigt at line centre with zero pressure broadening should match a
  //    pure Gaussian.
  {
    const double sigma = 0.05;
    const double v0 = voigt<double>(0.0, sigma, 1e-9);
    const double g0 = gaussian<double>(0.0, sigma);
    // pseudo-Voigt is approximate (~1%) so loosen tolerance.
    assert(std::fabs(v0 - g0) / g0 < 0.05);
  }

  // 2. Voigt at line centre with zero Doppler should match a pure Lorentzian.
  {
    const double gamma = 0.05;
    const double v0 = voigt<double>(0.0, 1e-9, gamma);
    const double l0 = lorentz<double>(0.0, gamma);
    assert(std::fabs(v0 - l0) / l0 < 0.05);
  }

  // 3. Voigt is symmetric about the line centre.
  {
    const double s = 0.04, g = 0.06;
    const double left  = voigt<double>(-0.1, s, g);
    const double right = voigt<double>( 0.1, s, g);
    assert(std::fabs(left - right) / left < 1e-6);
  }

  // 4. Voigt approximately integrates to 1 over a wide grid (area-normalised).
  {
    const double s = 0.04, g = 0.05;
    double total = 0.0;
    const double dx = 0.001;
    for (double x = -10.0; x <= 10.0; x += dx) {
      total += voigt<double>(x, s, g) * dx;
    }
    // pseudo-Voigt is good to ~1%, so allow 3%.
    assert(std::fabs(total - 1.0) < 0.03);
  }

  // 5. Voigt with Dual<double> seed: derivative wrt sigma at line centre is
  //    the analytic derivative of a Gaussian (the dominant term in this regime).
  {
    using Dx = Dual<double>;
    const double sigma_v = 0.05;
    Dx x{0.0, 0.0};
    Dx sigma{sigma_v, 1.0};      // seed: d/d_sigma
    Dx gamma{1e-9, 0.0};
    Dx v = voigt(x, sigma, gamma);
    // The derivative should be finite and have the sign of -1/sigma at peak
    // (broader Gaussian -> lower peak).
    assert(std::isfinite(v.d));
    assert(v.d < 0.0);
  }

  return 0;
}
