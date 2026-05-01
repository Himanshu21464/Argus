#include "argus/hitran.hpp"

#include <charconv>
#include <fstream>
#include <istream>
#include <stdexcept>
#include <string>

namespace argus {

namespace {

// Strip leading whitespace from a string_view.
std::string_view ltrim(std::string_view s) {
  std::size_t i = 0;
  while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
  return s.substr(i);
}

// Parse an integer from a fixed-width slice. Returns std::nullopt on
// failure or an empty / all-whitespace slice.
std::optional<int> parse_int(std::string_view slice) {
  slice = ltrim(slice);
  if (slice.empty()) return std::nullopt;
  int value = 0;
  auto [ptr, ec] = std::from_chars(slice.data(),
                                   slice.data() + slice.size(), value);
  if (ec != std::errc{}) return std::nullopt;
  return value;
}

// Parse a double from a fixed-width slice. Accepts HITRAN's mix of
// fixed-point ("3001.034740") and scientific ("9.231E-25") formats.
//
// Note: std::from_chars<double> is supported in libstdc++ since GCC 11.
// We fall back to std::strtod to be safe across older compilers.
std::optional<double> parse_double(std::string_view slice) {
  slice = ltrim(slice);
  if (slice.empty()) return std::nullopt;
  std::string buf(slice);
  // Strip trailing whitespace.
  while (!buf.empty() && (buf.back() == ' ' || buf.back() == '\t' ||
                          buf.back() == '\r' || buf.back() == '\n')) {
    buf.pop_back();
  }
  if (buf.empty()) return std::nullopt;
  char* endp = nullptr;
  const double v = std::strtod(buf.c_str(), &endp);
  if (endp == buf.c_str()) return std::nullopt;
  return v;
}

}  // namespace

std::optional<HitranRecord> Hitran::parse_line(std::string_view raw) {
  // Trim trailing CR/LF so a Windows-encoded line is still acceptable.
  while (!raw.empty() && (raw.back() == '\r' || raw.back() == '\n')) {
    raw.remove_suffix(1);
  }
  // A complete HITRAN record is 160 chars. We need columns up through
  // delta_air (col 67) so 67 is the minimum length we accept.
  constexpr std::size_t kMinLen = 67;
  if (raw.size() < kMinLen) return std::nullopt;

  auto col = [&](std::size_t a, std::size_t b) {
    // 1-indexed inclusive HITRAN columns: substring [a-1, b)
    return raw.substr(a - 1, b - (a - 1));
  };

  HitranRecord rec;

  if (auto v = parse_int(col(1, 2));   v) rec.molecule_id = *v; else return std::nullopt;
  if (auto v = parse_int(col(3, 3));   v) rec.isotope_id  = *v; else return std::nullopt;
  if (auto v = parse_double(col(4, 15));  v) rec.line.nu0_cm        = *v; else return std::nullopt;
  if (auto v = parse_double(col(16, 25)); v) rec.line.intensity_cm  = *v; else return std::nullopt;
  // einstein_A in cols 26-35 ignored at M2.
  if (auto v = parse_double(col(36, 40)); v) rec.line.gamma_air_cm  = *v; else return std::nullopt;
  if (auto v = parse_double(col(41, 45)); v) rec.line.gamma_self_cm = *v; else return std::nullopt;
  if (auto v = parse_double(col(46, 55)); v) rec.line.E_lower_cm    = *v; else return std::nullopt;
  if (auto v = parse_double(col(56, 59)); v) rec.line.n_air         = *v; else return std::nullopt;
  if (auto v = parse_double(col(60, 67)); v) rec.line.delta_air_cm  = *v; else return std::nullopt;

  return rec;
}

std::vector<HitranRecord> Hitran::load(std::istream& is, int filter) {
  std::vector<HitranRecord> out;
  std::string line;
  while (std::getline(is, line)) {
    auto rec = parse_line(line);
    if (!rec) continue;
    if (filter > 0 && rec->molecule_id != filter) continue;
    out.push_back(*rec);
  }
  return out;
}

std::vector<HitranRecord> Hitran::load_file(const std::string& path,
                                            int filter) {
  std::ifstream is(path);
  if (!is) {
    throw std::runtime_error("Hitran::load_file: cannot open " + path);
  }
  return load(is, filter);
}

}  // namespace argus
