#include <cassert>
#include <cmath>

#include "argus/argus.hpp"

int main() {
  using namespace argus;
  using D = Dual<double>;

  // 1. Arithmetic. f(x) = x^2 + 3x at x=2 has value 10 and derivative 7.
  {
    D x{2.0, 1.0};
    D f = x * x + 3.0 * x;
    assert(std::fabs(f.v - 10.0) < 1e-12);
    assert(std::fabs(f.d - 7.0)  < 1e-12);
  }

  // 2. Division. f(x) = 1 / x at x=2 has value 0.5 and derivative -0.25.
  {
    D x{2.0, 1.0};
    D f = D{1.0, 0.0} / x;
    assert(std::fabs(f.v - 0.5)   < 1e-12);
    assert(std::fabs(f.d + 0.25)  < 1e-12);
  }

  // 3. exp/log/sqrt.
  {
    D x{1.0, 1.0};
    D ex = exp(x);
    assert(std::fabs(ex.v - std::exp(1.0)) < 1e-12);
    assert(std::fabs(ex.d - std::exp(1.0)) < 1e-12);

    D lx = log(D{2.0, 1.0});
    assert(std::fabs(lx.v - std::log(2.0)) < 1e-12);
    assert(std::fabs(lx.d - 0.5)           < 1e-12);

    D sx = sqrt(D{4.0, 1.0});
    assert(std::fabs(sx.v - 2.0)   < 1e-12);
    assert(std::fabs(sx.d - 0.25)  < 1e-12);
  }

  // 4. Chain rule: f(x) = exp(-x^2) at x=1 has derivative -2x*exp(-x^2)=-2/e.
  {
    D x{1.0, 1.0};
    D f = exp(-(x * x));
    const double expect_v = std::exp(-1.0);
    const double expect_d = -2.0 * std::exp(-1.0);
    assert(std::fabs(f.v - expect_v) < 1e-12);
    assert(std::fabs(f.d - expect_d) < 1e-12);
  }

  return 0;
}
