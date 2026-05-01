#include <cassert>

#include "argus/argus.hpp"

int main() {
  using namespace argus;

  ir::Graph g;
  auto a = g.add(ir::NodeKind::kAtmosphere, "atm");
  auto o = g.add(ir::NodeKind::kOpacity, "h2o");
  auto m = g.add(ir::NodeKind::kForwardModel, "transmission", {a, o});
  auto s = g.add(ir::NodeKind::kSpectrum, "out", {m});
  (void)s;

  assert(g.nodes().size() == 4);

  // Same graph topology + names → same content address. Different name → different.
  ir::Graph g2;
  auto a2 = g2.add(ir::NodeKind::kAtmosphere, "atm");
  auto o2 = g2.add(ir::NodeKind::kOpacity, "h2o");
  auto m2 = g2.add(ir::NodeKind::kForwardModel, "transmission", {a2, o2});
  auto s2 = g2.add(ir::NodeKind::kSpectrum, "out", {m2});
  (void)s2;
  assert(g.content_address() == g2.content_address());

  ir::Graph g3;
  auto a3 = g3.add(ir::NodeKind::kAtmosphere, "different_atm");
  (void)a3;
  assert(g.content_address() != g3.content_address());

  return 0;
}
