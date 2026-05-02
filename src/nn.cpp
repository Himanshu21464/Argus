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

// ─── AffineCoupling ───────────────────────────────────────────────────

AffineCoupling::AffineCoupling(std::size_t dim,
                               std::size_t split,
                               std::vector<std::size_t> hidden_dims,
                               Activation act)
    : dim_(dim),
      split_(split),
      conditioner_(split, 2 * (dim - split),
                   std::move(hidden_dims), act) {
  if (dim_ < 2) {
    throw std::invalid_argument("AffineCoupling: dim must be >= 2");
  }
  if (split_ < 1 || split_ >= dim_) {
    throw std::invalid_argument(
        "AffineCoupling: split must be in [1, dim - 1]");
  }
}

void AffineCoupling::init_xavier(std::uint64_t seed) {
  conditioner_.init_xavier(seed);
}

AffineCoupling::Output AffineCoupling::forward(
    const std::vector<double>& x) const {
  if (x.size() != dim_) {
    throw std::invalid_argument(
        "AffineCoupling::forward: input size must equal dim");
  }
  const std::size_t n_b = dim_ - split_;
  const auto split_off = static_cast<std::ptrdiff_t>(split_);
  std::vector<double> x_a(x.begin(), x.begin() + split_off);
  std::vector<double> x_b(x.begin() + split_off, x.end());

  // Conditioner output: [s_0, s_1, ..., s_{n_b-1}, t_0, ..., t_{n_b-1}]
  std::vector<double> st = conditioner_.forward(x_a);

  Output out;
  out.y.resize(dim_);
  for (std::size_t i = 0; i < split_; ++i) out.y[i] = x_a[i];

  double log_det = 0.0;
  for (std::size_t i = 0; i < n_b; ++i) {
    const double s = st[i];
    const double t = st[n_b + i];
    out.y[split_ + i] = x_b[i] * std::exp(s) + t;
    log_det += s;
  }
  out.log_det_jacobian = log_det;
  return out;
}

AffineCoupling::Output AffineCoupling::inverse(
    const std::vector<double>& y) const {
  if (y.size() != dim_) {
    throw std::invalid_argument(
        "AffineCoupling::inverse: input size must equal dim");
  }
  const std::size_t n_b = dim_ - split_;
  const auto split_off = static_cast<std::ptrdiff_t>(split_);
  // y_a == x_a since the passive half is unchanged.
  std::vector<double> x_a(y.begin(), y.begin() + split_off);
  std::vector<double> st = conditioner_.forward(x_a);

  Output out;
  out.y.resize(dim_);
  for (std::size_t i = 0; i < split_; ++i) out.y[i] = x_a[i];

  double log_det = 0.0;
  for (std::size_t i = 0; i < n_b; ++i) {
    const double s = st[i];
    const double t = st[n_b + i];
    out.y[split_ + i] = (y[split_ + i] - t) * std::exp(-s);
    log_det += -s;       // log-det of inverse map
  }
  out.log_det_jacobian = log_det;
  return out;
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
