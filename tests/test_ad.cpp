// Reverse-mode autograd tests.
//
// Cross-validates against forward-mode Dual<T> (already test_finite_diff
// validates Dual against central FD). Once both autograds agree, we
// have two independent gradient implementations that triangulate
// each other.

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "argus/argus.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

}  // namespace

int main() {
  using namespace argus;
  using ad::Tape;
  using ad::Var;

  // ─── 1. Arithmetic: f(x, y) = x*y + 3*x at x=2, y=4. ────────────────
  {
    Tape t;
    Var x = t.input(2.0);
    Var y = t.input(4.0);
    Var f = x * y + x * 3.0;
    t.backward(f);
    // f = 2*4 + 2*3 = 14
    assert(close(f.val, 14.0, 1.0e-12));
    // df/dx = y + 3 = 7
    // df/dy = x   = 2
    assert(close(t.grad(x), 7.0, 1.0e-12));
    assert(close(t.grad(y), 2.0, 1.0e-12));
  }

  // ─── 2. Division: f(x) = 1/x at x=2. f=0.5; df/dx = -1/x^2 = -0.25. ─
  {
    Tape t;
    Var x = t.input(2.0);
    Var f = 1.0 / x;
    t.backward(f);
    assert(close(f.val, 0.5, 1.0e-12));
    assert(close(t.grad(x), -0.25, 1.0e-12));
  }

  // ─── 3. exp / log / sqrt. ───────────────────────────────────────────
  {
    Tape t;
    Var x = t.input(1.5);
    Var f = ad::exp(ad::log(ad::sqrt(x)));
    t.backward(f);
    // exp(log(sqrt(x))) = sqrt(x); df/dx = 0.5/sqrt(x)
    assert(close(f.val, std::sqrt(1.5), 1.0e-12));
    assert(close(t.grad(x), 0.5 / std::sqrt(1.5), 1.0e-10));
  }

  // ─── 4. pow + chain rule. f(x) = (x^2 + 1)^3. ───────────────────────
  // df/dx = 3 * (x^2+1)^2 * 2x = 6x(x^2+1)^2
  {
    Tape t;
    Var x = t.input(1.5);
    Var f = ad::pow(x * x + 1.0, 3.0);
    t.backward(f);
    const double y = 1.5 * 1.5 + 1.0;     // 3.25
    assert(close(f.val, y * y * y, 1.0e-10));
    assert(close(t.grad(x), 6.0 * 1.5 * y * y, 1.0e-10));
  }

  // ─── 5. sin / cos chain. f(x) = sin(cos(x)) at x = 0.7. ─────────────
  // df/dx = cos(cos(x)) * (-sin(x))
  {
    Tape t;
    Var x = t.input(0.7);
    Var f = ad::sin(ad::cos(x));
    t.backward(f);
    assert(close(f.val, std::sin(std::cos(0.7)), 1.0e-12));
    const double expected =
        std::cos(std::cos(0.7)) * (-std::sin(0.7));
    assert(close(t.grad(x), expected, 1.0e-12));
  }

  // ─── 6. Cross-validation against forward-mode Dual<T>: 5-D function
  //     f(x) = x0^2 + x1*sin(x2) + exp(x3*x4)
  //     for several random points. Both autograds should agree to <1e-10.
  {
    auto f_dual = [](const std::vector<Dual<double>>& x) {
      return x[0] * x[0] + x[1] * sin(x[2]) + exp(x[3] * x[4]);
    };
    auto f_tape = [](const std::vector<Var>& x) {
      return x[0] * x[0] + x[1] * ad::sin(x[2]) + ad::exp(x[3] * x[4]);
    };
    const std::vector<std::vector<double>> points = {
      {1.0, 2.0, 3.0, 0.5, -0.3},
      {-2.0, 0.5, -1.2, 0.7, 0.1},
      {0.1, -0.5, 2.0, -0.4, 0.6},
    };
    for (const auto& pt : points) {
      // Forward-mode gradient.
      auto g_fwd = grad(f_dual, pt);
      // Reverse-mode gradient via Tape.
      Tape t;
      std::vector<Var> xs;
      for (double v : pt) xs.push_back(t.input(v));
      Var y = f_tape(xs);
      t.backward(y);
      for (std::size_t i = 0; i < pt.size(); ++i) {
        assert(close(t.grad(xs[i]), g_fwd[i], 1.0e-10));
      }
    }
  }

  // ─── 7. Many-input scaling: a sum-of-squares of 100 vars.
  //     d(sum x_i^2)/d(x_k) = 2 x_k.
  //     Reverse-mode does this in one forward + one backward; a
  //     forward-mode equivalent would be 100x slower.
  {
    Tape t;
    std::vector<Var> xs;
    for (std::size_t i = 0; i < 100; ++i) {
      xs.push_back(t.input(0.1 * static_cast<double>(i)));
    }
    Var sum = t.input(0.0);
    for (const auto& v : xs) sum = sum + v * v;
    t.backward(sum);
    for (std::size_t i = 0; i < 100; ++i) {
      assert(close(t.grad(xs[i]),
                   2.0 * 0.1 * static_cast<double>(i),
                   1.0e-10, 1.0e-12));
    }
  }

  // ─── 8. Reset clears state; reusing the tape after reset works. ─────
  {
    Tape t;
    {
      Var x = t.input(3.0);
      Var f = x * x;
      t.backward(f);
      assert(close(t.grad(x), 6.0, 1.0e-12));
    }
    t.reset();
    {
      Var x = t.input(5.0);
      Var f = x * x * x;
      t.backward(f);
      assert(close(t.grad(x), 75.0, 1.0e-10));
    }
  }

  // ─── 9. Bad inputs throw. ───────────────────────────────────────────
  {
    Tape t1, t2;
    Var a = t1.input(1.0);
    bool threw = false;
    try { t2.backward(a); }   // wrong tape
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    threw = false;
    try { (void)t2.grad(a); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);
  }

  // ─── 10. tanh derivative: d/dx tanh(x) = 1 - tanh(x)^2.
  {
    Tape t;
    Var x = t.input(0.5);
    Var f = ad::tanh(x);
    t.backward(f);
    const double th = std::tanh(0.5);
    assert(close(f.val, th, 1.0e-12));
    assert(close(t.grad(x), 1.0 - th * th, 1.0e-12));
  }

  // ─── 11. Mixing Vars from two different tapes throws on every binary
  //     operator. Earlier code silently produced wrong gradients
  //     because the second Var's index was resolved against the
  //     first's tape storage. ────────────────────────────────────────
  {
    Tape t1, t2;
    Var a = t1.input(2.0);
    Var b = t2.input(3.0);
    bool threw;
    threw = false;
    try { (void)(a + b); } catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    threw = false;
    try { (void)(a - b); } catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    threw = false;
    try { (void)(a * b); } catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    threw = false;
    try { (void)(a / b); } catch (const std::runtime_error&) { threw = true; }
    assert(threw);
  }

  return 0;
}
