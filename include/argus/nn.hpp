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

// Affine coupling layer (Real NVP, Dinh et al. 2017; refined as Glow,
// Kingma & Dhariwal 2018). The fundamental building block of
// normalizing-flow posteriors used in amortized SBI.
//
// For input x ∈ R^D, split into x_a (first `split` dims) and x_b
// (remaining D-split dims). The forward transform is:
//
//     y_a = x_a
//     y_b = x_b · exp(s(x_a)) + t(x_a)
//
// where s and t are arbitrary MLPs of input width `split` and output
// width `D - split`. Both share an `Sequential` of input dim `split`
// and output dim `2 * (D - split)` (concatenated [s, t]).
//
// Jacobian determinant (forward direction) is sum_i exp(s_i(x_a)) since
// only y_b depends on x_b in a separable way → log|det J| = Σ s_i(x_a).
//
// Stacking N AffineCoupling layers (alternating which half is passive)
// gives a Real NVP normalizing flow. Sampling: draw z ~ N(0, I), apply
// the inverse stack to get a sample from the learned distribution.
class AffineCoupling {
 public:
  // dim   — total input/output dimensionality D
  // split — number of dimensions in the "passive" half x_a (must be in
  //         [1, dim - 1])
  // hidden_dims — MLP hidden widths for the conditioner
  // act   — activation between hidden layers (default Tanh, smooth)
  AffineCoupling(std::size_t dim,
                 std::size_t split,
                 std::vector<std::size_t> hidden_dims,
                 Activation act = Activation::Tanh);

  std::size_t dim()   const noexcept { return dim_; }
  std::size_t split() const noexcept { return split_; }

  // Initialise the conditioner with Xavier; deterministic seed.
  void init_xavier(std::uint64_t seed = 0);

  // Direct access to the underlying conditioner network for
  // pretrained-weight loading.
  Sequential& conditioner() noexcept { return conditioner_; }
  const Sequential& conditioner() const noexcept { return conditioner_; }

  // Forward: z = f(x). Also returns the log-determinant of the Jacobian
  // |df/dx| so the flow can be used as a normalizing flow for density
  // computation: log p_z(z) = log p_x(x) - log|det J|.
  struct Output {
    std::vector<double> y;
    double log_det_jacobian;
  };
  Output forward(const std::vector<double>& x) const;

  // Inverse: x = f^{-1}(z). Returns the same struct; log_det_jacobian
  // is the log-det of the INVERSE map, i.e. -forward.log_det_jacobian.
  Output inverse(const std::vector<double>& y) const;

 private:
  std::size_t dim_;
  std::size_t split_;
  Sequential conditioner_;     // (split) -> 2 * (dim - split) for [s, t]
};

// Permutation: deterministically swap halves so successive coupling
// layers can transform the previously-passive dimensions. Trivial
// log-det = 0 since permutations preserve volume.
struct HalfSwap {
  std::size_t dim;
  std::size_t split;
  std::vector<double> apply(const std::vector<double>& x) const;
};

// Stack of AffineCoupling layers separated by half-swap permutations,
// trained against a standard-Gaussian base distribution. The standard
// Real NVP / Glow architecture for amortized SBI in exoplanet
// retrieval pipelines (see e.g. DINGO for gravitational-wave PE).
//
// Forward direction: x (data) -> z (base). z is N(0, I) under a
// well-trained flow.
// Inverse direction: z (base) -> x (sample from learned distribution).
//
// log_density(x) = log p_z(z) + log|det df/dx|  (change-of-variables)
//                = -0.5 * |z|^2 - 0.5 * D * log(2π) + Σ log_det_i
//
// Sampling: z ~ N(0, I), then x = inverse(z). Useful for amortized
// SBI: train the flow on (params, observation) pairs to learn
// p(params | observation), sample posterior conditional on a new
// observation by inverting.
class NormalizingFlow {
 public:
  // dim          — data dimensionality
  // n_couplings  — number of coupling layers; each is followed by a
  //                half-swap so dimensions get equal time as "active"
  // split        — initial split index (default dim/2)
  // hidden_dims  — conditioner MLP hidden widths
  // act          — conditioner activation
  NormalizingFlow(std::size_t dim,
                  std::size_t n_couplings,
                  std::size_t split,
                  std::vector<std::size_t> hidden_dims,
                  Activation act = Activation::Tanh);

  std::size_t dim()         const noexcept { return dim_; }
  std::size_t n_couplings() const noexcept { return couplings_.size(); }

  void init_xavier(std::uint64_t seed = 0);

  // Mutable access for pretrained-weight loading.
  AffineCoupling& coupling(std::size_t i);
  const AffineCoupling& coupling(std::size_t i) const;

  // Forward: x (data) -> z (base). Returns z and the cumulative
  // log_det_jacobian Σ log|det J_i|.
  AffineCoupling::Output forward(const std::vector<double>& x) const;

  // Inverse: z (base) -> x (data). Returns x and the cumulative
  // log_det_jacobian (negative of forward).
  AffineCoupling::Output inverse(const std::vector<double>& z) const;

  // log p_x(x) under the flow (with N(0, I) base):
  //   log_p_x = -0.5 |z|^2 - 0.5 D ln(2π) + log|det dz/dx|
  double log_density(const std::vector<double>& x) const;

  // Draw a sample from the flow given a base sample z ~ N(0, I).
  // Caller controls the RNG (test determinism, parallel sampling).
  std::vector<double> sample(std::mt19937_64& rng) const;

 private:
  std::size_t dim_;
  std::size_t init_split_;
  std::vector<AffineCoupling> couplings_;
  std::vector<HalfSwap> swaps_;
};

}  // namespace argus::nn
