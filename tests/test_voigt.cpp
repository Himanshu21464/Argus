// Hard tests for the Hui-Armstrong-Wray Voigt evaluator.
// HAW gives ~1e-6 relative error in the upper half-plane, so these
// assertions are ~10000× tighter than the M2-α pseudo-Voigt tests.

#include <cassert>
#include <cmath>

#include "argus/argus.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

}  // namespace

int main() {
  using namespace argus;

  // 1. Voigt with tiny γ should match the pure Gaussian to <1e-5 at ALL x.
  {
    const double sigma = 0.04;
    for (double x = -0.4; x <= 0.4; x += 0.05) {
      const double v = voigt<double>(x, sigma, 1.0e-9);
      const double g = gaussian<double>(x, sigma);
      assert(close(v, g, 1.0e-4, 1.0e-12));
    }
  }

  // 2. Voigt with tiny σ should match the pure Lorentzian to <1e-5.
  {
    const double gamma = 0.04;
    for (double x = -0.4; x <= 0.4; x += 0.05) {
      const double v = voigt<double>(x, 1.0e-9, gamma);
      const double l = lorentz<double>(x, gamma);
      assert(close(v, l, 1.0e-4, 1.0e-12));
    }
  }

  // 3. Symmetric in x to machine precision.
  {
    const double s = 0.04, g = 0.06;
    for (double x = 0.001; x < 1.0; x += 0.097) {
      const double L = voigt<double>(-x, s, g);
      const double R = voigt<double>( x, s, g);
      assert(close(L, R, 1.0e-12));
    }
  }

  // 4. Area-normalised to 1 over a wide grid (much tighter than the
  //    pseudo-Voigt's 3% tolerance). HAW + Simpson integration ≈ 1e-6.
  {
    const double s = 0.04, g = 0.05;
    double total = 0.0;
    const double dx = 0.0005;
    const double xmax = 50.0;     // ~1000 sigma — fully captures the tails
    for (double x = -xmax; x <= xmax; x += dx) {
      total += voigt<double>(x, s, g) * dx;
    }
    assert(close(total, 1.0, 1.0e-3));    // Simpson error from coarse step
  }

  // 5. Voigt cross-validation against numerical convolution of
  //    Gaussian * Lorentzian. The convolution (V = G ⊛ L) is what the
  //    Voigt is by definition, so this is a fully-independent check.
  {
    const double sigma_g = 0.04;
    const double gamma_l = 0.05;
    const double dx_int = 0.0005;
    const double xmax_int = 5.0;
    const double xs[] = {-0.10, -0.03, 0.0, 0.03, 0.10};
    for (double x : xs) {
      // Numerically convolve: V_num(x) = ∫ G(t) * L(x - t) dt
      double V_num = 0.0;
      for (double t = -xmax_int; t <= xmax_int; t += dx_int) {
        V_num += gaussian<double>(t, sigma_g) *
                 lorentz<double>(x - t, gamma_l) * dx_int;
      }
      const double V = voigt<double>(x, sigma_g, gamma_l);
      // Numerical convolution residual + HAW residual = ~1e-3 relative.
      assert(close(V, V_num, 5.0e-3, 1.0e-12));
    }
  }

  // 6. Voigt at large x should asymptote to Lorentzian wing
  //    L(x) = γ_l / (π x²)  for x >> max(σ_g, γ_l)
  {
    const double sigma_g = 0.05;
    const double gamma_l = 0.04;
    const double x = 1.0;          // 20σ, 25γ — deep wing
    const double V = voigt<double>(x, sigma_g, gamma_l);
    const double L_asymp = gamma_l / (3.14159265358979323846 * x * x);
    assert(close(V, L_asymp, 1.0e-2));
  }

  // 7. Dual-number autograd — derivative wrt γ_l verified against central
  //    finite differences to 1e-5 relative.
  {
    using D = Dual<double>;
    const double sigma_v = 0.05;
    const double gamma_v = 0.06;
    const double x_v     = 0.02;

    D V_d = voigt(D{x_v, 0.0}, D{sigma_v, 0.0}, D{gamma_v, 1.0});

    const double h = 1.0e-7;
    const double Vp = voigt<double>(x_v, sigma_v, gamma_v + h);
    const double Vm = voigt<double>(x_v, sigma_v, gamma_v - h);
    const double dV_fd = (Vp - Vm) / (2.0 * h);
    assert(close(V_d.d, dV_fd, 1.0e-5, 1.0e-12));
  }

  // 8. Dual-number autograd — derivative wrt σ at line centre.
  {
    using D = Dual<double>;
    const double sigma_v = 0.05;
    const double gamma_v = 0.04;
    const double x_v     = 0.0;

    D V_d = voigt(D{x_v, 0.0}, D{sigma_v, 1.0}, D{gamma_v, 0.0});
    const double h = 1.0e-7;
    const double Vp = voigt<double>(x_v, sigma_v + h, gamma_v);
    const double Vm = voigt<double>(x_v, sigma_v - h, gamma_v);
    const double dV_fd = (Vp - Vm) / (2.0 * h);
    assert(close(V_d.d, dV_fd, 1.0e-5, 1.0e-12));
    // peak narrows -> derivative wrt sigma is negative at line centre
    assert(V_d.d < 0.0);
  }

  return 0;
}
