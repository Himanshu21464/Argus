#include "argus/nn.hpp"

#include <cmath>
#include <stdexcept>

namespace argus::nn {

Linear::Linear(std::size_t in_dim, std::size_t out_dim)
    : in_dim_(in_dim), out_dim_(out_dim),
      weight_(in_dim * out_dim, 0.0),
      bias_(out_dim, 0.0) {
  if (in_dim_ == 0 || out_dim_ == 0) {
    throw std::invalid_argument("Linear: dimensions must be positive");
  }
}

void Linear::init_xavier(std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  const double bound = std::sqrt(6.0 / static_cast<double>(in_dim_ + out_dim_));
  std::uniform_real_distribution<double> u(-bound, bound);
  for (auto& w : weight_) w = u(rng);
  std::fill(bias_.begin(), bias_.end(), 0.0);
}

void Linear::set_weights(std::vector<double> w) {
  if (w.size() != weight_.size()) {
    throw std::invalid_argument(
        "Linear::set_weights: size must equal in_dim * out_dim");
  }
  weight_ = std::move(w);
}

void Linear::set_bias(std::vector<double> b) {
  if (b.size() != bias_.size()) {
    throw std::invalid_argument(
        "Linear::set_bias: size must equal out_dim");
  }
  bias_ = std::move(b);
}

std::vector<double> Linear::forward(const std::vector<double>& x) const {
  if (x.size() != in_dim_) {
    throw std::invalid_argument(
        "Linear::forward: input size must equal in_dim");
  }
  std::vector<double> y(out_dim_);
  // y[i] = bias[i] + sum_j weight[i, j] * x[j]
  for (std::size_t i = 0; i < out_dim_; ++i) {
    double acc = bias_[i];
    const double* wrow = weight_.data() + i * in_dim_;
    for (std::size_t j = 0; j < in_dim_; ++j) {
      acc += wrow[j] * x[j];
    }
    y[i] = acc;
  }
  return y;
}

// ─── Sequential ───────────────────────────────────────────────────────

Sequential::Sequential(std::size_t in_dim,
                       std::size_t out_dim,
                       std::vector<std::size_t> hidden_dims,
                       Activation hidden_act)
    : in_dim_(in_dim), out_dim_(out_dim), hidden_act_(hidden_act) {
  if (in_dim_ == 0 || out_dim_ == 0) {
    throw std::invalid_argument("Sequential: dimensions must be positive");
  }
  std::size_t prev = in_dim_;
  for (std::size_t h : hidden_dims) {
    if (h == 0) {
      throw std::invalid_argument("Sequential: hidden dim must be positive");
    }
    layers_.emplace_back(prev, h);
    prev = h;
  }
  layers_.emplace_back(prev, out_dim_);
}

void Sequential::init_xavier(std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  for (std::size_t i = 0; i < layers_.size(); ++i) {
    // Use a derived seed per layer so each Linear sees a unique stream
    // while the whole network remains deterministic given the master seed.
    layers_[i].init_xavier(rng());
  }
}

std::vector<double> Sequential::forward(
    const std::vector<double>& x) const {
  if (x.size() != in_dim_) {
    throw std::invalid_argument(
        "Sequential::forward: input size must equal in_dim");
  }
  std::vector<double> h = x;
  for (std::size_t i = 0; i < layers_.size(); ++i) {
    h = layers_[i].forward(h);
    // Activation on every hidden layer; final Linear is left raw.
    if (i + 1 < layers_.size()) apply_activation(h, hidden_act_);
  }
  return h;
}

Linear& Sequential::layer(std::size_t i) {
  if (i >= layers_.size()) {
    throw std::out_of_range("Sequential::layer: index out of range");
  }
  return layers_[i];
}

const Linear& Sequential::layer(std::size_t i) const {
  if (i >= layers_.size()) {
    throw std::out_of_range("Sequential::layer: index out of range");
  }
  return layers_[i];
}

}  // namespace argus::nn
