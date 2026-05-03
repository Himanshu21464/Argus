// Argus example 02 — Voigt-shaped line opacity for an H2O-like atmosphere.
//
// Builds a 4-line H2O line list, attaches a LineListOpacity to the
// transmission model, and prints the resulting transmission spectrum
// over a coarse JWST-NIRSpec-like wavenumber grid. The spectrum now
// shows the absorption peaks at the line centres rather than the flat
// profile of the M1 GreyOpacity stub.

#include <iomanip>
#include <iostream>

#include "argus/argus.hpp"

int main() {
  using namespace argus;

  std::cout << "Argus " << version_string() << " — example 02 (Voigt H2O)\n\n";

  Species h2o{"H2O", 18.015};
  Atmosphere atm = isothermal(/*T=*/1200.0,
                              /*P_top=*/1.0e-6,
                              /*P_bot=*/1.0e2,
                              /*n_layers=*/60,
                              h2o,
                              /*VMR=*/1.0e-3);

  // Toy H2O-like line list (real units, cartoon strengths/positions).
  // Loading real HITRAN .par files into this same struct is shown in
  // examples/03_real_hitran.cpp via `argus::Hitran::load_file`.
  std::vector<Line> lines = {
    Line{2950.0, 8.0e-21, 0.06, 0.0, 0.50, 0.0, 0.0},
    Line{3050.0, 1.2e-20, 0.05, 0.0, 0.45, 0.0, 0.0},
    Line{3650.0, 6.0e-21, 0.05, 0.0, 0.55, 0.0, 0.0},
    Line{3750.0, 9.0e-21, 0.06, 0.0, 0.50, 0.0, 0.0},
  };
  auto opacity = std::make_shared<LineListOpacity>("H2O", lines, 18.015);

  TransmissionModel model;
  model.add_opacity(opacity);

  std::vector<double> wn;
  const double wn_lo = 2900.0;
  const double wn_hi = 3800.0;
  const int n_pts = 46;
  for (int i = 0; i < n_pts; ++i) {
    wn.push_back(wn_lo + i * (wn_hi - wn_lo) / (n_pts - 1));
  }

  Spectrum s = model.forward(atm, wn);

  std::cout << std::fixed << std::setprecision(2);
  std::cout << "wavenumber[cm-1]   transit_depth (R_p/R_*)^2 * 1e3\n";
  for (std::size_t i = 0; i < s.wavenumber_cm.size(); ++i) {
    std::cout << "  " << std::setw(7) << s.wavenumber_cm[i]
              << "         "
              << std::setprecision(4) << (s.values[i] * 1.0e3)
              << std::setprecision(2) << "\n";
  }

  // Demonstrate Dual<double> through a single Voigt evaluation:
  //   dV/d(gamma_l) at the line centre with sigma_g=0.04, gamma_l=0.05
  using D = Dual<double>;
  D x{0.0, 0.0};
  D sg{0.04, 0.0};
  D gl{0.05, 1.0};                  // seed derivative w.r.t. gamma_l
  D V = voigt(x, sg, gl);
  std::cout << "\n[autograd] V(0; 0.04, 0.05)        = " << V.v << "\n";
  std::cout <<   "[autograd] dV/d(gamma_l) at center = " << V.d << "\n";

  return 0;
}
