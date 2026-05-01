// Tests for self-broadening (γ_self) and pressure-shift (δ_air) effects
// in the line opacity computation.

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

  // 1. With γ_self ≠ γ_air, increasing self-VMR should shift the cross-
  //    section noticeably. Use γ_self = 2 * γ_air to make the effect
  //    obvious.
  {
    Line line{
      /*nu0=*/3000.0,    /*S=*/1.0e-21,
      /*gamma_air=*/0.05, /*gamma_self=*/0.10,
      /*n_air=*/0.5,     /*delta_air=*/0.0,
      /*E_lower=*/0.0
    };
    LineListOpacity lo("X", {line}, 18.0);

    // Sample at a wing point where Lorentzian dominates — shift in γ_l
    // produces clean signal there.
    const std::vector<double> wn{3000.20};   // 0.20 cm⁻¹ off centre
    const std::vector<double> T{500.0};
    const std::vector<double> P{1.0};

    Tensor s_air  = lo.cross_section_with_self(wn, T, P, {0.0});
    Tensor s_half = lo.cross_section_with_self(wn, T, P, {0.5});
    Tensor s_pure = lo.cross_section_with_self(wn, T, P, {1.0});

    // monotone — wider Lorentz means more wing absorption
    assert(s_air[0]  < s_half[0]);
    assert(s_half[0] < s_pure[0]);

    // pure-self has γ_self = 2*γ_air → wing absorption ~2× larger
    assert(close(s_pure[0] / s_air[0], 2.0, 0.20));
  }

  // 2. Default cross_section() = self-VMR=0 (i.e. air-broadened limit).
  {
    Line line{3000.0, 1.0e-21, 0.05, 0.10, 0.5, 0.0, 0.0};
    LineListOpacity lo("X", {line}, 18.0);
    const std::vector<double> wn{3000.20};
    Tensor a = lo.cross_section(wn, {500.0}, {1.0});
    Tensor b = lo.cross_section_with_self(wn, {500.0}, {1.0}, {0.0});
    assert(close(a[0], b[0], 1.0e-12));
  }

  // 3. Pressure shift: with δ_air > 0, the effective line centre moves
  //    to higher wavenumber by δ_air * P_atm. At a point that was the
  //    centre at low P, increasing P should reduce the cross-section.
  {
    Line line{
      /*nu0=*/3000.0,    /*S=*/1.0e-21,
      /*gamma_air=*/0.05, /*gamma_self=*/0.05,
      /*n_air=*/0.5,     /*delta_air=*/0.1,    // 0.1 cm⁻¹/atm shift
      /*E_lower=*/0.0
    };
    LineListOpacity lo("X", {line}, 18.0);

    // Sample exactly at the unshifted centre.
    const std::vector<double> wn{3000.0};
    const std::vector<double> T{500.0};

    Tensor lo_p = lo.cross_section(wn, T, {0.001});  // 1 mbar — negligible shift
    Tensor hi_p = lo.cross_section(wn, T, {2.0});    // 2 bar — line shifted ~0.2 cm⁻¹

    // At higher pressure the line centre is shifted away from our
    // sample point, so the cross-section there drops.
    assert(hi_p[0] < lo_p[0]);
  }

  // 4. With δ_air = 0, pressure does not shift the centre but does
  //    broaden the line. Sample at the centre: cross-section grows
  //    weakly with P (Lorentzian peak height = 1/(π γ) shrinks, but
  //    Voigt at centre with both contributions changes too).
  {
    Line line{3000.0, 1.0e-21, 0.05, 0.05, 0.5, 0.0, 0.0};
    LineListOpacity lo("X", {line}, 18.0);
    const std::vector<double> wn{3000.0};
    const std::vector<double> T{500.0};

    Tensor lo_p = lo.cross_section(wn, T, {0.001});
    Tensor hi_p = lo.cross_section(wn, T, {2.0});

    // Centre point with no shift: increasing pressure broadens γ, so
    // peak height drops (γ grows linearly in P). Verify this.
    assert(hi_p[0] < lo_p[0]);
  }

  return 0;
}
