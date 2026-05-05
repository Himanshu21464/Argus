#include "argus/jwst_data.hpp"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace argus {

std::vector<double> JWSTSpectrum::wavenumber_cm() const {
  std::vector<double> out;
  out.reserve(wavelength_um.size());
  for (double w : wavelength_um) {
    if (!(w > 0.0)) {
      throw std::runtime_error(
          "JWSTSpectrum::wavenumber_cm: wavelength must be > 0");
    }
    out.push_back(10000.0 / w);
  }
  return out;
}

namespace JWST {

namespace {

// Trim leading + trailing whitespace.
std::string_view trim(std::string_view s) {
  std::size_t i = 0;
  while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r')) ++i;
  std::size_t j = s.size();
  while (j > i && (s[j - 1] == ' ' || s[j - 1] == '\t' || s[j - 1] == '\r')) --j;
  return s.substr(i, j - i);
}

double parse_double_or_throw(std::string_view s, std::size_t row_idx) {
  std::string buf(s);
  char* endp = nullptr;
  const double v = std::strtod(buf.c_str(), &endp);
  if (endp == buf.c_str()) {
    throw std::runtime_error(
        "JWST::load: failed to parse number on row " +
        std::to_string(row_idx) + ": '" + buf + "'");
  }
  return v;
}

}  // namespace

JWSTSpectrum load(std::istream& is, std::string source) {
  JWSTSpectrum out;
  out.source = std::move(source);

  std::string line;
  std::size_t row = 0;
  while (std::getline(is, line)) {
    ++row;
    std::string_view sv = trim(line);
    if (sv.empty()) continue;
    if (sv.front() == '#') continue;       // comment

    // Split on commas — exactly 3 fields required.
    std::vector<std::string_view> fields;
    std::size_t start = 0;
    for (std::size_t i = 0; i <= sv.size(); ++i) {
      if (i == sv.size() || sv[i] == ',') {
        fields.push_back(trim(sv.substr(start, i - start)));
        start = i + 1;
      }
    }
    if (fields.size() != 3) {
      throw std::runtime_error(
          "JWST::load: row " + std::to_string(row) +
          " has " + std::to_string(fields.size()) +
          " fields (expected 3: wavelength_um, depth, sigma)");
    }
    out.wavelength_um.push_back(parse_double_or_throw(fields[0], row));
    out.transit_depth.push_back(parse_double_or_throw(fields[1], row));
    out.sigma_depth.push_back(parse_double_or_throw(fields[2], row));
  }

  if (out.wavelength_um.empty()) {
    throw std::runtime_error("JWST::load: no usable rows in input");
  }
  return out;
}

JWSTSpectrum load_file(const std::string& path, std::string source) {
  std::ifstream is(path);
  if (!is) {
    throw std::runtime_error("JWST::load_file: cannot open '" + path + "'");
  }
  return load(is, std::move(source));
}

}  // namespace JWST

}  // namespace argus
