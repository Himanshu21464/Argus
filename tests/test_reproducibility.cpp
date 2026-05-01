// Bit-exact reproducibility — running the same forward model on the
// same inputs must produce identical outputs every time. Catches
// nondeterminism from uninitialised memory, floating-point reductions
// in nondeterministic order, hash randomisation, etc.

#include <cassert>
#include <cmath>
#include <sstream>

#include "argus/argus.hpp"
#include "argus/test_data.hpp"

int main() {
  using namespace argus;

  // 1. Build the same atmosphere twice and verify byte-equal results
  //    from the forward model — including the IR content address.
  Species h2o{"H2O", 18.015};

  std::istringstream is{std::string(test_data::kH2OLines)};
  auto records = Hitran::load(is, /*filter=*/1);
  std::vector<Line> lines;
  for (const auto& r : records) lines.push_back(r.line);

  auto run_once = [&]() {
    Atmosphere atm = isothermal(/*T=*/1500.0,
                                /*P_top=*/1.0e-6,
                                /*P_bot=*/1.0e2,
                                /*n_layers=*/40,
                                h2o,
                                /*VMR=*/1.0e-3);
    auto op = std::make_shared<LineListOpacity>("H2O", lines, 18.015);
    TransmissionModel m;
    m.add_opacity(op);
    std::vector<double> wn;
    for (double w = 3500.0; w <= 7400.0; w += 25.0) wn.push_back(w);
    return m.forward(atm, wn);
  };

  Spectrum s1 = run_once();
  Spectrum s2 = run_once();

  assert(s1.values.size() == s2.values.size());
  for (std::size_t i = 0; i < s1.values.size(); ++i) {
    // bit-exact — not "close" — true reproducibility
    assert(s1.values[i] == s2.values[i]);
  }

  // 2. IR content address is deterministic across calls.
  auto build_graph = []() {
    ir::Graph g;
    auto a = g.add(ir::NodeKind::kAtmosphere,    "atm_v0");
    auto o = g.add(ir::NodeKind::kOpacity,        "h2o_hitran");
    auto m = g.add(ir::NodeKind::kForwardModel,   "transmission", {a, o});
    g.add(ir::NodeKind::kSpectrum, "spectrum", {m});
    return g.content_address();
  };
  assert(build_graph() == build_graph());

  // 3. Cross-section at a single (nu, T, P) is deterministic.
  {
    auto op = std::make_shared<LineListOpacity>("H2O", lines, 18.015);
    Tensor a = op->cross_section({3651.97}, {1500.0}, {1.0});
    Tensor b = op->cross_section({3651.97}, {1500.0}, {1.0});
    assert(a[0] == b[0]);
  }

  // 4. Q(T) is deterministic.
  {
    const double q1 = Partition::Q("H2O", 1500.0);
    const double q2 = Partition::Q("H2O", 1500.0);
    assert(q1 == q2);
  }

  return 0;
}
