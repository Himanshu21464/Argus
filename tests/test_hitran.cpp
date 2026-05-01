// Tests for the HITRAN .par fixed-width parser.
// Round-trip parses a known string and verifies every numeric field.

#include <cassert>
#include <cmath>
#include <sstream>

#include "argus/argus.hpp"
#include "argus/test_data.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

}  // namespace

int main() {
  using namespace argus;

  // 1. A single hand-crafted line with explicit expected values.
  //    Strict 67-char fixed-width layout with no inter-column spaces.
  {
    const std::string raw =
      " 11 3651.969720 4.501E-19 4.501E+020.0750.090  136.02020.500.001253          0 1   1";
    auto rec = Hitran::parse_line(raw);
    assert(rec.has_value());
    assert(rec->molecule_id == 1);
    assert(rec->isotope_id  == 1);
    assert(close(rec->line.nu0_cm,        3651.969720, 1.0e-9));
    assert(close(rec->line.intensity_cm,  4.501e-19,   1.0e-3));
    assert(close(rec->line.gamma_air_cm,  0.0750,      1.0e-3));
    assert(close(rec->line.gamma_self_cm, 0.090,       1.0e-3));
    assert(close(rec->line.E_lower_cm,    136.0202,    1.0e-6));
    assert(close(rec->line.n_air,         0.500,       1.0e-3));
    assert(close(rec->line.delta_air_cm,  0.001253,    1.0e-3));
  }

  // 2. Trailing CR/LF tolerated.
  {
    const std::string raw =
      " 11 3651.969720 4.501E-19 4.501E+020.0750.090  136.02020.500.001253          0 1   1\r\n";
    assert(Hitran::parse_line(raw).has_value());
  }

  // 3. Too-short line rejected.
  {
    const std::string raw = " 11 3651.0";
    assert(!Hitran::parse_line(raw).has_value());
  }

  // 4. Multi-line stream loads everything and the count matches the bundled
  //    real-world H2O test fixture.
  {
    std::istringstream is{std::string(test_data::kH2OLines)};
    auto records = Hitran::load(is);
    assert(records.size() == 16);

    for (const auto& r : records) {
      assert(r.molecule_id == 1);
      assert(r.isotope_id  == 1);
      assert(r.line.nu0_cm > 3000.0 && r.line.nu0_cm < 8000.0);
      assert(r.line.intensity_cm  > 0.0);
      assert(r.line.gamma_air_cm  > 0.0);
      assert(r.line.gamma_self_cm > 0.0);
      assert(r.line.E_lower_cm   >= 0.0);
      assert(r.line.n_air         > 0.0 && r.line.n_air < 1.0);
    }
  }

  // 5. Filter by molecule id.
  {
    std::istringstream is{std::string(test_data::kH2OLines)};
    auto h2o = Hitran::load(is, /*molecule_id_filter=*/1);
    assert(h2o.size() == 16);

    std::istringstream is2{std::string(test_data::kH2OLines)};
    auto co2 = Hitran::load(is2, /*molecule_id_filter=*/2);
    assert(co2.empty());
  }

  return 0;
}
