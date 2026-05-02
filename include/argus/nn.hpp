#pragma once

#include <cstddef>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "tensor.hpp"

namespace argus::nn {

// Lightweight, header-light neural-net primitives sized for the
// normalizing-flow / amortized-SBI stack. Pure C++20, no BLAS/CUDA
// dependency. M3.5 will swap inner loops for cuBLAS/cuDNN where
// available; the shape API stays stable.
//
// Convention: a "vector input" is a std::vector<double> of length
// in_dim; a "batch input" is a std::vector<std::vector<double>> of
// length batch_size, each entry of length in_dim. Outputs follow the
// same convention with out_dim.

// Activation function tag.
enum class Activation {
  None,
  ReLU,
  LeakyReLU,        // alpha = 0.01
  Tanh,
  Sigmoid,
};

// Apply the activation in-place on a vector.
inline void apply_activation(std::vector<double>& x, Activation a) {
  switch (a) {
    case Activation::None:
      break;
    case Activation::ReLU:
      for (auto& v : x) if (v < 0.0) v = 0.0;
      break;
    case Activation::LeakyReLU:
      for (auto& v : x) if (v < 0.0) v *= 0.01;
      break;
    case Activation::Tanh:
      for (auto& v : x) v = std::tanh(v);
      break;
    case Activation::Sigmoid:
      for (auto& v : x) v = 1.0 / (1.0 + std::exp(-v));
      break;
  }
}

// Linear layer: y = W x + b.
//   W is row-major [out_dim x in_dim], stored flat in `weight_`.
//   b is [out_dim].
class Linear {
 public:
  Linear(std::size_t in_dim, std::size_t out_dim);

  std::size_t in_dim()  const noexcept { return in_dim_; }
  std::size_t out_dim() const noexcept { return out_dim_; }

  // Glorot / Xavier uniform initialisation: U(-sqrt(6/(in+out)), +sqrt(...)).
  // Bias zero. Deterministic given a seed.
  void init_xavier(std::uint64_t seed = 0);

  // Direct setters for round-tripping and testing.
  void set_weights(std::vector<double> w);
  void set_bias(std::vector<double> b);
  const std::vector<double>& weights() const noexcept { return weight_; }
  const std::vector<double>& bias()    const noexcept { return bias_;   }

  // Forward pass on a single vector (returns a new vector).
  std::vector<double> forward(const std::vector<double>& x) const;

 private:
  std::size_t in_dim_, out_dim_;
  std::vector<double> weight_;     // [out_dim x in_dim] row-major
  std::vector<double> bias_;        // [out_dim]
};

// A simple MLP: alternating Linear + Activation layers terminating
// in a final Linear with no activation. The architecture is fixed at
// construction; weights can be initialised with Xavier and accessed
// for round-tripping.
class Sequential {
 public:
  // hidden_dims describes the widths of each hidden layer (excluding
  // input and output). `hidden_act` is applied after every hidden
  // Linear; the final layer has no activation.
  Sequential(std::size_t in_dim,
             std::size_t out_dim,
             std::vector<std::size_t> hidden_dims,
             Activation hidden_act = Activation::Tanh);

  std::size_t in_dim()  const noexcept { return in_dim_; }
  std::size_t out_dim() const noexcept { return out_dim_; }
  std::size_t n_layers() const noexcept { return layers_.size(); }

  // Initialise all layers with Xavier; bias zero. Deterministic seed.
  void init_xavier(std::uint64_t seed = 0);

  std::vector<double> forward(const std::vector<double>& x) const;

  // Mutable access to a specific layer (for tests and pretrained-weight
  // loading). Caller must use indices in [0, n_layers()).
  Linear& layer(std::size_t i);
  const Linear& layer(std::size_t i) const;

 private:
  std::size_t in_dim_;
  std::size_t out_dim_;
  std::vector<Linear> layers_;
  Activation hidden_act_;
};

}  // namespace argus::nn
