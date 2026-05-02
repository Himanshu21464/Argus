#include "argus/chain_io.hpp"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "argus/version.hpp"

namespace argus {

namespace {

// Format a double as a hex-float for bit-exact round-trips.
// std::format would be cleaner but we keep the C++20 baseline portable.
std::string hex_double(double x) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%a", x);
  return buf;
}

double parse_hex_or_dec(const std::string& s) {
  // strtod handles both decimal and hex-float.
  char* end = nullptr;
  const double v = std::strtod(s.c_str(), &end);
  if (end == s.c_str()) {
    throw std::runtime_error("ChainIO::load_csv: failed to parse double '" +
                             s + "'");
  }
  return v;
}

}  // namespace

void ChainIO::save_csv(const std::string& path,
                       const std::vector<Parameter>& params,
                       const std::vector<std::vector<double>>& samples,
                       const std::vector<double>& log_posteriors) {
  if (samples.size() != log_posteriors.size()) {
    throw std::invalid_argument(
        "ChainIO::save_csv: samples.size() must equal log_posteriors.size()");
  }
  if (!samples.empty() && samples.front().size() != params.size()) {
    throw std::invalid_argument(
        "ChainIO::save_csv: sample dimension must match parameter count");
  }

  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("ChainIO::save_csv: cannot open '" + path + "'");
  }

  out << "# Argus chain " << version_string() << "\n";
  out << "# n_samples=" << samples.size()
      << " n_dim=" << params.size() << "\n";
  for (std::size_t i = 0; i < params.size(); ++i) {
    if (i > 0) out << ',';
    out << params[i].name;
  }
  if (!params.empty()) out << ',';
  out << "log_posterior\n";

  for (std::size_t s = 0; s < samples.size(); ++s) {
    for (std::size_t d = 0; d < params.size(); ++d) {
      if (d > 0) out << ',';
      out << hex_double(samples[s][d]);
    }
    if (!params.empty()) out << ',';
    out << hex_double(log_posteriors[s]) << '\n';
  }
}

void ChainIO::save_csv(const std::string& path,
                       const std::vector<Parameter>& params,
                       const Retrieval::Result& result) {
  save_csv(path, params, result.samples, result.log_posteriors);
}

LoadedChain ChainIO::load_csv(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("ChainIO::load_csv: cannot open '" + path + "'");
  }

  LoadedChain out;
  std::string line;
  bool got_header = false;

  while (std::getline(in, line)) {
    // Strip trailing CR for CRLF files
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    if (line.front() == '#') continue;     // skip comments

    if (!got_header) {
      // Parse comma-separated header, drop final 'log_posterior' column.
      std::stringstream ss(line);
      std::string field;
      std::vector<std::string> fields;
      while (std::getline(ss, field, ',')) fields.push_back(field);
      if (fields.empty() || fields.back() != "log_posterior") {
        throw std::runtime_error(
            "ChainIO::load_csv: header must end with 'log_posterior'");
      }
      out.param_names.assign(fields.begin(), fields.end() - 1);
      got_header = true;
      continue;
    }

    // Data row.
    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> fields;
    while (std::getline(ss, field, ',')) fields.push_back(field);
    if (fields.size() != out.param_names.size() + 1) {
      throw std::runtime_error(
          "ChainIO::load_csv: row has wrong number of columns");
    }
    std::vector<double> sample;
    sample.reserve(out.param_names.size());
    for (std::size_t d = 0; d < out.param_names.size(); ++d) {
      sample.push_back(parse_hex_or_dec(fields[d]));
    }
    out.samples.push_back(std::move(sample));
    out.log_posteriors.push_back(parse_hex_or_dec(fields.back()));
  }

  if (!got_header) {
    throw std::runtime_error(
        "ChainIO::load_csv: file contains no usable data");
  }
  return out;
}

}  // namespace argus
