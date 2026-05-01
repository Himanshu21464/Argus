// Argus example 03 — load real HITRAN H2O records, run a transmission
// spectrum on a hot-Jupiter-like atmosphere, dump the result.
//
// What's new vs example 02:
//   * uses argus::Hitran::load() to parse real .par-format records
//     from the bundled test fixture (test_data::kH2OLines)
//   * uses argus::Partition for proper Q(T) intensity scaling
//   * spans both the 2.7 μm and 1.4 μm H2O bands so the spectrum
//     shape mirrors what JWST/NIRSpec sees on hot Jupiters

#include <iomanip>
#include <iostream>
#include <sstream>

#include "argus/argus.hpp"
#include "argus/test_data.hpp"

int main() {
  using namespace argus;

  std::cout << "Argus " << version_string() << " — example 03 (real HITRAN H2O)\n\n";

  // 1. Load the bundled real-world H2O records from the .par fixture.
  std::istringstream is{std::string(test_data::kH2OLines)};
  auto records = Hitran::load(is, /*molecule_id=*/1);
  std::cout << "loaded " << records.size() << " H2O lines from HITRAN fixture\n";

  std::vector<Line> lines;
  lines.reserve(records.size());
  for (const auto& r : records) lines.push_back(r.line);

  auto opacity = std::make_shared<LineListOpacity>("H2O", lines, 18.015);

  // 2. Atmosphere — a typical hot Jupiter (T=1500 K, log P from -6 to 2 bar).
  Species h2o{"H2O", 18.015};
  Atmosphere atm = isothermal(/*T=*/1500.0,
                              /*P_top=*/1.0e-6,
                              /*P_bot=*/1.0e2,
                              /*n_layers=*/80,
                              h2o,
                              /*VMR=*/1.0e-3);

  // 3. Forward model.
  TransmissionModel model;
  model.add_opacity(opacity);

  // 4. JWST-NIRSpec-PRISM-like wavenumber grid: ~1.3-2.8 μm => 3500-7500 cm^-1
  std::vector<double> wn;
  for (double w = 3550.0; w <= 7400.0; w += 50.0) wn.push_back(w);

  Spectrum s = model.forward(atm, wn);

  std::cout << "\nQ(H2O, 1500 K) = " << Partition::Q("H2O", 1500.0)
            << "  (Q_ref = " << Partition::Q_ref("H2O") << ")\n\n";

  std::cout << std::fixed << std::setprecision(2);
  std::cout << "wavenumber [cm-1]    wavelength [um]    transit depth [ppm]\n";
  std::cout << "-------------------  ---------------    ------------------\n";
  for (std::size_t i = 0; i < s.wavenumber_cm.size(); ++i) {
    const double um = 1.0e4 / s.wavenumber_cm[i];
    std::cout << "  " << std::setw(8) << s.wavenumber_cm[i]
              << "          " << std::setw(7) << um
              << "          " << std::setprecision(1)
              << (s.values[i] * 1.0e6) << " ppm\n"
              << std::setprecision(2);
  }

  return 0;
}
