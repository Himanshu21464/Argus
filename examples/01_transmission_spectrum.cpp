// Argus example 01 — transmission spectrum of an isothermal H2O atmosphere.
//
// This is the M1 smoke test: build an atmosphere, attach a (placeholder)
// opacity kernel, run the forward model, print the spectrum.
//
// M3 swaps GreyOpacity for HITRAN-backed Voigt evaluation; the call site
// here does not change.

#include <iostream>

#include "argus/argus.hpp"

int main() {
  using namespace argus;

  std::cout << "Argus " << version_string() << " — example 01\n";

  // 1. Atmosphere.
  Species h2o{"H2O", 18.015};
  Atmosphere atm = isothermal(/*T_k=*/1200.0,
                              /*P_top_bar=*/1.0e-6,
                              /*P_bot_bar=*/1.0e2,
                              /*n_layers=*/40,
                              h2o,
                              /*mixing_ratio=*/1.0e-3);

  // 2. Opacity kernel (placeholder — real version is HITRAN+CUDA Voigt).
  auto opacity = std::make_shared<GreyOpacity>("H2O", /*sigma_cm2=*/1.0e-22);

  // 3. Forward model.
  TransmissionModel model;
  model.add_opacity(opacity);

  // 4. Wavenumber grid (cm^-1) — coarse JWST-NIRSpec-like band.
  std::vector<double> wavenumber_cm;
  for (int i = 0; i < 21; ++i) {
    wavenumber_cm.push_back(2000.0 + 100.0 * i);
  }

  Spectrum s = model.forward(atm, wavenumber_cm);

  // 5. Print.
  std::cout << "wavenumber[cm-1]  transit_depth\n";
  for (std::size_t i = 0; i < s.wavenumber_cm.size(); ++i) {
    std::cout << "  " << s.wavenumber_cm[i] << "          " << s.values[i] << "\n";
  }

  // 6. Argus IR demo: build a tiny graph and print its content address.
  ir::Graph g;
  auto a_node = g.add(ir::NodeKind::kAtmosphere, "atm");
  auto o_node = g.add(ir::NodeKind::kOpacity, "h2o_grey");
  auto m_node = g.add(ir::NodeKind::kForwardModel, "transmission",
                      {a_node, o_node});
  auto sp_node = g.add(ir::NodeKind::kSpectrum, "out", {m_node});
  (void)sp_node;
  std::cout << "ir.content_address = " << g.content_address() << "\n";
  return 0;
}
