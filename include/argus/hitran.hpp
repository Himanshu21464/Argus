#pragma once

#include <iosfwd>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "line_list.hpp"

namespace argus {

// Loader for HITRAN .par files (the fixed-width 160-character format).
//
// The .par format is documented at:
//   https://hitran.org/docs/definitions-and-units/
// Each record is exactly 160 characters with this column layout
// (1-indexed, inclusive):
//
//   1 -  2  molecule_id   (i2)
//   3       isotope_id    (i1)
//   4 - 15  nu0           (f12.6, cm^-1)
//   16- 25  intensity     (e10.3, cm^-1 / (molec cm^-2))
//   26- 35  einstein_A    (e10.3, s^-1)
//   36- 40  gamma_air     (f5.4, cm^-1 / atm)
//   41- 45  gamma_self    (f5.4, cm^-1 / atm)
//   46- 55  E_lower       (f10.4, cm^-1)
//   56- 59  n_air         (f4.2)
//   60- 67  delta_air     (f8.6, cm^-1 / atm)
//   ... (quantum / error / ref columns ignored at M2)
//
// The parser is tolerant of trailing whitespace / line endings but rejects
// lines that are too short for the columns we read.
struct HitranRecord {
  int molecule_id = 0;
  int isotope_id  = 0;
  Line line;
};

class Hitran {
 public:
  // Parse one .par-format record. Returns std::nullopt on malformed input.
  static std::optional<HitranRecord> parse_line(std::string_view raw);

  // Read every line from `is`. If `molecule_id_filter > 0` only records
  // matching that molecule id are kept.
  static std::vector<HitranRecord> load(std::istream& is,
                                        int molecule_id_filter = -1);

  // Same, from a file path.
  static std::vector<HitranRecord> load_file(const std::string& path,
                                             int molecule_id_filter = -1);
};

}  // namespace argus
