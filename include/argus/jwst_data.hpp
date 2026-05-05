#pragma once

#include <istream>
#include <string>
#include <vector>

// Real JWST exoplanet transmission-spectrum loader.
//
// Wraps a single observed spectrum: wavelength, transit depth, and the
// per-bin uncertainty. CSV format is the lingua franca of published
// exoplanet spectra — each major reduction pipeline (FIREFLy, Tiberius,
// Eureka!, Tshirt) emits the same three columns, occasionally with a
// few comment lines on top. This loader tolerates `#`-prefixed
// comments and trims whitespace.
//
// Wavelength is converted to wavenumber on access since the kernel's
// opacity surface is wavenumber-native.

namespace argus {

struct JWSTSpectrum {
  std::string source;                 // citation tag (paper / DOI / instrument)
  std::vector<double> wavelength_um;  // [n_bin]
  std::vector<double> transit_depth;  // [n_bin] dimensionless (Rp/R*)^2
  std::vector<double> sigma_depth;    // [n_bin] 1-sigma uncertainty

  std::size_t size() const noexcept { return wavelength_um.size(); }

  // Wavenumber-native view (cm^-1) for argus opacity / RT calls.
  // ν[cm^-1] = 10000 / λ[μm]. Returned in same order as wavelength_um.
  std::vector<double> wavenumber_cm() const;
};

namespace JWST {

// Parse a CSV stream of real JWST spectrum data. Format:
//   - Optional `#`-prefixed comment lines (skipped)
//   - One row per wavelength bin, three comma-separated columns:
//       wavelength_um, transit_depth, sigma_depth
//   - Trailing whitespace, CRLF and empty lines are tolerated.
// Throws std::runtime_error on malformed rows.
JWSTSpectrum load(std::istream& is, std::string source = "");

// Convenience wrapper: open a file at `path` and call load().
JWSTSpectrum load_file(const std::string& path, std::string source = "");

}  // namespace JWST

}  // namespace argus
