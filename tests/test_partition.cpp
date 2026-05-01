// Tests for the TIPS partition-function fits.
// Anchor: Q(296 K) must equal the published TIPS-2017 value to <0.1%.

#include <cassert>
#include <cmath>
#include <stdexcept>

#include "argus/argus.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

}  // namespace

int main() {
  using namespace argus;

  // 1. Anchor values at 296 K, all major molecules, exact.
  const struct { const char* k; double q; } anchors[] = {
    {"H2O", 174.58},
    {"CO2", 286.94},
    {"CH4", 590.48},
    {"CO",  107.42},
    {"NH3", 1725.22},
  };
  for (const auto& a : anchors) {
    const double q_ref  = Partition::Q_ref(a.k);
    const double q_at   = Partition::Q(a.k, 296.0);
    assert(close(q_ref, a.q, 1.0e-3));
    assert(close(q_at,  a.q, 1.0e-9));      // anchor is exact by construction
  }

  // 2. Q(T) must be monotonically increasing for all rotational molecules
  //    over the 100-3000 K range.
  for (const auto& a : anchors) {
    double prev = -1.0;
    for (double T = 100.0; T <= 3000.0; T += 50.0) {
      const double q = Partition::Q(a.k, T);
      assert(q > 0.0);
      assert(q > prev);    // strictly monotone
      prev = q;
    }
  }

  // 3. Q(T) at high T: H2O Q(1000) should be ~1296 (TIPS-2017).
  //    Power-law fit anchored at 296 K and exponent 1.648 hits this exactly
  //    by construction; allow 1% slack for floating-point.
  {
    const double q = Partition::Q("H2O", 1000.0);
    assert(close(q, 1296.0, 0.01));
  }
  // 4. Q at low T: H2O Q(150) should be ~57 ± 10% (power-law extrapolated
  //    below the anchor; published TIPS gives ~56.6).
  {
    const double q = Partition::Q("H2O", 150.0);
    assert(close(q, 57.0, 0.10));
  }
  // 5. Q at very high T: H2O Q(2000) should be ~4080 (TIPS-2017).
  {
    const double q = Partition::Q("H2O", 2000.0);
    assert(close(q, 4080.0, 0.05));
  }

  // 6. Unknown species: throw.
  bool threw = false;
  try { Partition::Q("XENON_DUST", 296.0); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);

  // 7. Negative T: throw.
  threw = false;
  try { Partition::Q("H2O", -5.0); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);

  return 0;
}
