// HITRAN file I/O round-trip — write the bundled test data to /tmp,
// load it back via Hitran::load_file, and verify field-by-field
// equality with the in-memory parse.

#include <cassert>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "argus/argus.hpp"
#include "argus/test_data.hpp"

namespace {

bool close(double a, double b, double rtol, double atol = 0.0) {
  return std::fabs(a - b) <= atol + rtol * std::fabs(b);
}

}  // namespace

int main() {
  using namespace argus;

  // 1. Write bundled fixture to a temp file.
  const std::string path = "/tmp/argus_h2o_test.par";
  {
    std::ofstream out(path, std::ios::trunc);
    out << test_data::kH2OLines;
  }

  // 2. Load via load_file and via load(istringstream); compare.
  auto from_file = Hitran::load_file(path, /*filter=*/1);
  std::istringstream is{std::string(test_data::kH2OLines)};
  auto from_mem = Hitran::load(is, /*filter=*/1);

  assert(from_file.size() == from_mem.size());
  assert(from_file.size() == 16);

  for (std::size_t i = 0; i < from_file.size(); ++i) {
    const auto& a = from_file[i];
    const auto& b = from_mem[i];
    assert(a.molecule_id == b.molecule_id);
    assert(a.isotope_id  == b.isotope_id);
    assert(close(a.line.nu0_cm,        b.line.nu0_cm,        1.0e-12));
    assert(close(a.line.intensity_cm,  b.line.intensity_cm,  1.0e-12));
    assert(close(a.line.gamma_air_cm,  b.line.gamma_air_cm,  1.0e-12));
    assert(close(a.line.gamma_self_cm, b.line.gamma_self_cm, 1.0e-12));
    assert(close(a.line.E_lower_cm,    b.line.E_lower_cm,    1.0e-12));
    assert(close(a.line.n_air,         b.line.n_air,         1.0e-12));
    assert(close(a.line.delta_air_cm,  b.line.delta_air_cm,  1.0e-12));
  }

  // 3. Missing file throws.
  bool threw = false;
  try { (void)Hitran::load_file("/tmp/this_file_should_not_exist_argus.par"); }
  catch (const std::runtime_error&) { threw = true; }
  assert(threw);

  // 4. Cleanup.
  std::remove(path.c_str());

  return 0;
}
