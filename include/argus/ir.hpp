#pragma once

#include <memory>
#include <string>
#include <variant>
#include <vector>

namespace argus::ir {

// The Argus IR is a typed graph of physics passes:
//   inputs  : Atmosphere, OpacityTable, InstrumentPSF, ObservationGeometry
//   nodes   : compile-time-checked transformations between these types
//   outputs : Spectrum
//
// M1 ships only the type tags and a graph skeleton; M3-6 fills in the
// pass implementations and the autograd tape.

enum class NodeKind {
  kAtmosphere,
  kOpacity,
  kInstrumentPSF,
  kGeometry,
  kPriors,
  kForwardModel,
  kSpectrum,
  kPosterior,
};

struct Node {
  // First enum value (kAtmosphere) as a defensive default — Graph::add
  // always overwrites this immediately, but a default-constructed Node
  // would otherwise have an indeterminate kind (UB on read). Adding a
  // sentinel kNone would shift the underlying integer values of the
  // others and change every Graph's content_address(), so we keep the
  // existing values stable and just pick the lowest as default.
  NodeKind kind = NodeKind::kAtmosphere;
  std::string name;
  std::vector<std::shared_ptr<Node>> inputs;
};

class Graph {
 public:
  std::shared_ptr<Node> add(NodeKind kind, std::string name,
                            std::vector<std::shared_ptr<Node>> inputs = {});

  const std::vector<std::shared_ptr<Node>>& nodes() const noexcept {
    return nodes_;
  }

  // Returns a deterministic content-address (hex SHA-style) of the graph
  // structure — same graph topology + same node names → same address.
  std::string content_address() const;

 private:
  std::vector<std::shared_ptr<Node>> nodes_;
};

}  // namespace argus::ir
