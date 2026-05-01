#include "argus/ir.hpp"

#include <cstdint>
#include <sstream>

namespace argus::ir {

std::shared_ptr<Node> Graph::add(NodeKind kind, std::string name,
                                 std::vector<std::shared_ptr<Node>> inputs) {
  auto node = std::make_shared<Node>();
  node->kind = kind;
  node->name = std::move(name);
  node->inputs = std::move(inputs);
  nodes_.push_back(node);
  return node;
}

std::string Graph::content_address() const {
  // Deterministic 64-bit FNV-1a over the canonicalized node list.
  // Production version will use SHA-256 over a stable serialization.
  constexpr std::uint64_t kFnvOffset = 1469598103934665603ULL;
  constexpr std::uint64_t kFnvPrime = 1099511628211ULL;

  std::uint64_t h = kFnvOffset;
  auto mix_byte = [&](unsigned char b) {
    h ^= b;
    h *= kFnvPrime;
  };
  auto mix_string = [&](const std::string& s) {
    for (char c : s) mix_byte(static_cast<unsigned char>(c));
    mix_byte(0);
  };

  for (const auto& n : nodes_) {
    mix_byte(static_cast<unsigned char>(n->kind));
    mix_string(n->name);
    for (const auto& in : n->inputs) {
      mix_string(in->name);
    }
    mix_byte(0);
  }

  std::ostringstream oss;
  oss << std::hex << h;
  return oss.str();
}

}  // namespace argus::ir
