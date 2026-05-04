// Argus example 01 — transmission spectrum of an isothermal H2O atmosphere.
//
// This is the smallest end-to-end demo: build an atmosphere, attach a
// baseline (wavelength-flat) opacity kernel, run the forward model,
// print the spectrum. For the HITRAN-backed Voigt-line variant see
// examples/03_real_hitran.cpp — same call site, different OpacityKernel.

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

  // 2. Opacity kernel — wavelength-flat baseline (see example 03 for HITRAN).
  auto opacity = std::make_shared<GreyOpacity>("H2O", /*sigma_cm2=*/1.0e-22);

  // 3. Forward model.
  TransmissionModel model;
  model.add_opacity(opacity);

  // 4. Wavenumber grid (cm^-1) — coarse JWST-NIRSpec-like band.
  Spectrum s = model.forward(atm, make_grid(2000.0, 4000.0, 21));

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
