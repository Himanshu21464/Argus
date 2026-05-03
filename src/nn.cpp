#include "argus/nn.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "argus/version.hpp"

namespace argus::nn {

namespace {

std::string hex_dbl(double x) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%a", x);
  return buf;
}

double parse_dbl(const std::string& s) {
  char* end = nullptr;
  const double v = std::strtod(s.c_str(), &end);
  if (end == s.c_str()) {
    throw std::runtime_error("NormalizingFlow::load: bad number '" + s + "'");
  }
  return v;
}

}  // namespace

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

// ─── HalfSwap + NormalizingFlow ───────────────────────────────────────

std::vector<double> HalfSwap::apply(const std::vector<double>& x) const {
  // Cycle: [a, b] -> [b, a] when split * 2 == dim;
  // For unequal splits, treat as a left rotation by `split`.
  if (x.size() != dim) {
    throw std::invalid_argument("HalfSwap: input size must equal dim");
  }
  std::vector<double> y(dim);
  // y = x rotated left by `split`
  for (std::size_t i = 0; i < dim; ++i) {
    y[i] = x[(i + split) % dim];
  }
  return y;
}

NormalizingFlow::NormalizingFlow(std::size_t dim,
                                 std::size_t n_couplings,
                                 std::size_t split,
                                 std::vector<std::size_t> hidden_dims,
                                 Activation act)
    : dim_(dim), init_split_(split) {
  if (dim_ < 2) {
    throw std::invalid_argument("NormalizingFlow: dim must be >= 2");
  }
  if (n_couplings == 0) {
    throw std::invalid_argument(
        "NormalizingFlow: n_couplings must be >= 1");
  }
  if (split < 1 || split >= dim_) {
    throw std::invalid_argument(
        "NormalizingFlow: split must be in [1, dim - 1]");
  }
  couplings_.reserve(n_couplings);
  swaps_.reserve(n_couplings);
  for (std::size_t i = 0; i < n_couplings; ++i) {
    couplings_.emplace_back(dim_, split, hidden_dims, act);
    swaps_.push_back({dim_, split});
  }
}

void NormalizingFlow::init_xavier(std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  for (auto& c : couplings_) c.init_xavier(rng());
}

AffineCoupling& NormalizingFlow::coupling(std::size_t i) {
  if (i >= couplings_.size()) {
    throw std::out_of_range(
        "NormalizingFlow::coupling: index out of range");
  }
  return couplings_[i];
}

const AffineCoupling& NormalizingFlow::coupling(std::size_t i) const {
  if (i >= couplings_.size()) {
    throw std::out_of_range(
        "NormalizingFlow::coupling: index out of range");
  }
  return couplings_[i];
}

AffineCoupling::Output NormalizingFlow::forward(
    const std::vector<double>& x) const {
  if (x.size() != dim_) {
    throw std::invalid_argument(
        "NormalizingFlow::forward: input size must equal dim");
  }
  AffineCoupling::Output state;
  state.y = x;
  state.log_det_jacobian = 0.0;
  for (std::size_t i = 0; i < couplings_.size(); ++i) {
    auto step = couplings_[i].forward(state.y);
    state.y = std::move(step.y);
    state.log_det_jacobian += step.log_det_jacobian;
    // Permutation between layers (no log-det contribution).
    if (i + 1 < couplings_.size()) {
      state.y = swaps_[i].apply(state.y);
    }
  }
  return state;
}

AffineCoupling::Output NormalizingFlow::inverse(
    const std::vector<double>& z) const {
  if (z.size() != dim_) {
    throw std::invalid_argument(
        "NormalizingFlow::inverse: input size must equal dim");
  }
  AffineCoupling::Output state;
  state.y = z;
  state.log_det_jacobian = 0.0;
  // Reverse order; un-permute, then invert coupling.
  for (std::size_t r = 0; r < couplings_.size(); ++r) {
    const std::size_t i = couplings_.size() - 1 - r;
    if (i + 1 < couplings_.size()) {
      // Un-rotate: rotating left by `dim - split` is the inverse
      // of rotating left by `split`.
      HalfSwap inv_swap{dim_, dim_ - swaps_[i].split};
      state.y = inv_swap.apply(state.y);
    }
    auto step = couplings_[i].inverse(state.y);
    state.y = std::move(step.y);
    state.log_det_jacobian += step.log_det_jacobian;
  }
  return state;
}

double NormalizingFlow::log_density(const std::vector<double>& x) const {
  auto fw = forward(x);
  double zz = 0.0;
  for (double v : fw.y) zz += v * v;
  constexpr double kLog2Pi = 1.8378770664093453;   // ln(2π)
  return -0.5 * zz
         - 0.5 * static_cast<double>(dim_) * kLog2Pi
         + fw.log_det_jacobian;
}

std::vector<double> NormalizingFlow::sample(std::mt19937_64& rng) const {
  std::normal_distribution<double> nd(0.0, 1.0);
  std::vector<double> z(dim_);
  for (auto& v : z) v = nd(rng);
  return inverse(z).y;
}

void NormalizingFlow::save(const std::string& path) const {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("NormalizingFlow::save: cannot open '" +
                             path + "'");
  }
  out << "# argus.nn.NormalizingFlow " << argus::version_string() << "\n";
  out << "# dim=" << dim_
      << " n_couplings=" << couplings_.size()
      << " init_split=" << init_split_ << "\n";
  for (std::size_t k = 0; k < couplings_.size(); ++k) {
    const auto& cond = couplings_[k].conditioner();
    out << "# coupling " << k
        << " conditioner_layers=" << cond.n_layers() << "\n";
    for (std::size_t l = 0; l < cond.n_layers(); ++l) {
      const auto& lin = cond.layer(l);
      out << "layer " << l << " in=" << lin.in_dim()
          << " out=" << lin.out_dim() << "\n";
      for (double w : lin.weights()) out << hex_dbl(w) << "\n";
      out << "bias\n";
      for (double b : lin.bias()) out << hex_dbl(b) << "\n";
    }
  }
}

void NormalizingFlow::load(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("NormalizingFlow::load: cannot open '" +
                             path + "'");
  }
  std::string line;
  std::size_t k = 0, l = 0;
  bool have_layer_header = false;
  std::size_t expected_in = 0, expected_out = 0;
  std::vector<double> w_buf;
  std::vector<double> b_buf;
  bool reading_bias = false;
  std::size_t expected_w = 0;

  auto commit_layer = [&]() {
    if (k >= couplings_.size() || l >= couplings_[k].conditioner().n_layers()) {
      throw std::runtime_error(
          "NormalizingFlow::load: more layers in file than in flow");
    }
    auto& lin = couplings_[k].conditioner().layer(l);
    if (lin.in_dim() != expected_in || lin.out_dim() != expected_out) {
      throw std::runtime_error(
          "NormalizingFlow::load: layer shape mismatch at coupling " +
          std::to_string(k) + " layer " + std::to_string(l));
    }
    if (w_buf.size() != lin.weights().size()) {
      throw std::runtime_error(
          "NormalizingFlow::load: weight count mismatch at coupling " +
          std::to_string(k) + " layer " + std::to_string(l));
    }
    if (b_buf.size() != lin.bias().size()) {
      throw std::runtime_error(
          "NormalizingFlow::load: bias count mismatch at coupling " +
          std::to_string(k) + " layer " + std::to_string(l));
    }
    lin.set_weights(std::move(w_buf));
    lin.set_bias(std::move(b_buf));
    w_buf.clear();
    b_buf.clear();
    have_layer_header = false;
    ++l;
    if (l >= couplings_[k].conditioner().n_layers()) {
      l = 0;
      ++k;
    }
  };

  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    if (line.front() == '#') continue;

    if (line.rfind("layer ", 0) == 0) {
      if (have_layer_header) commit_layer();
      // Parse "layer L in=N out=M". The L index is discarded — layer
      // ordering is implicit in the file.
      std::stringstream ss(line);
      std::string token;
      ss >> token;       // "layer"
      std::size_t layer_index_unused;
      ss >> layer_index_unused;
      (void)layer_index_unused;
      std::string in_kv, out_kv;
      ss >> in_kv >> out_kv;
      if (in_kv.rfind("in=", 0) != 0 || out_kv.rfind("out=", 0) != 0) {
        throw std::runtime_error(
            "NormalizingFlow::load: malformed layer header");
      }
      expected_in  = std::stoul(in_kv.substr(3));
      expected_out = std::stoul(out_kv.substr(4));
      expected_w   = expected_in * expected_out;
      reading_bias = false;
      have_layer_header = true;
      continue;
    }
    if (line == "bias") { reading_bias = true; continue; }
    if (!have_layer_header) continue;

    const double v = parse_dbl(line);
    if (!reading_bias) {
      if (w_buf.size() < expected_w) {
        w_buf.push_back(v);
      } else {
        throw std::runtime_error(
            "NormalizingFlow::load: extra weight before 'bias' marker");
      }
    } else {
      b_buf.push_back(v);
    }
  }
  if (have_layer_header) commit_layer();

  // Completeness check: a truncated file (fewer layers than the flow
  // expects) was previously silently accepted, leaving the remaining
  // couplings in their pre-load state. Catch that as an error.
  if (k != couplings_.size() || l != 0) {
    throw std::runtime_error(
        "NormalizingFlow::load: file ended before all flow layers were "
        "filled (loaded " + std::to_string(k) + "/" +
        std::to_string(couplings_.size()) + " complete couplings, " +
        "partial layer " + std::to_string(l) + ")");
  }
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

// ─── ConditionalAffineCoupling ───────────────────────────────────────

ConditionalAffineCoupling::ConditionalAffineCoupling(
    std::size_t dim, std::size_t cond_dim, std::size_t split,
    std::vector<std::size_t> hidden_dims, Activation act)
    : dim_(dim),
      cond_dim_(cond_dim),
      split_(split),
      conditioner_(split + cond_dim, 2 * (dim - split),
                   std::move(hidden_dims), act) {
  if (dim_ < 2) {
    throw std::invalid_argument("ConditionalAffineCoupling: dim must be >= 2");
  }
  if (split_ < 1 || split_ >= dim_) {
    throw std::invalid_argument(
        "ConditionalAffineCoupling: split must be in [1, dim - 1]");
  }
}

void ConditionalAffineCoupling::init_xavier(std::uint64_t seed) {
  conditioner_.init_xavier(seed);
}

namespace {

std::vector<double> concat(const std::vector<double>& a,
                           const std::vector<double>& b) {
  std::vector<double> out;
  out.reserve(a.size() + b.size());
  out.insert(out.end(), a.begin(), a.end());
  out.insert(out.end(), b.begin(), b.end());
  return out;
}

}  // namespace

AffineCoupling::Output ConditionalAffineCoupling::forward(
    const std::vector<double>& x, const std::vector<double>& cond) const {
  if (x.size() != dim_) {
    throw std::invalid_argument(
        "ConditionalAffineCoupling::forward: input size must equal dim");
  }
  if (cond.size() != cond_dim_) {
    throw std::invalid_argument(
        "ConditionalAffineCoupling::forward: cond size must equal cond_dim");
  }
  const std::size_t n_b = dim_ - split_;
  const auto split_off = static_cast<std::ptrdiff_t>(split_);
  std::vector<double> x_a(x.begin(), x.begin() + split_off);
  std::vector<double> x_b(x.begin() + split_off, x.end());

  // Conditioner input: [x_a; cond].
  std::vector<double> input = concat(x_a, cond);
  std::vector<double> st = conditioner_.forward(input);

  AffineCoupling::Output out;
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

AffineCoupling::Output ConditionalAffineCoupling::inverse(
    const std::vector<double>& y, const std::vector<double>& cond) const {
  if (y.size() != dim_) {
    throw std::invalid_argument(
        "ConditionalAffineCoupling::inverse: input size must equal dim");
  }
  if (cond.size() != cond_dim_) {
    throw std::invalid_argument(
        "ConditionalAffineCoupling::inverse: cond size must equal cond_dim");
  }
  const std::size_t n_b = dim_ - split_;
  const auto split_off = static_cast<std::ptrdiff_t>(split_);
  std::vector<double> x_a(y.begin(), y.begin() + split_off);
  std::vector<double> input = concat(x_a, cond);
  std::vector<double> st = conditioner_.forward(input);

  AffineCoupling::Output out;
  out.y.resize(dim_);
  for (std::size_t i = 0; i < split_; ++i) out.y[i] = x_a[i];

  double log_det = 0.0;
  for (std::size_t i = 0; i < n_b; ++i) {
    const double s = st[i];
    const double t = st[n_b + i];
    out.y[split_ + i] = (y[split_ + i] - t) * std::exp(-s);
    log_det += -s;
  }
  out.log_det_jacobian = log_det;
  return out;
}

// ─── ConditionalNormalizingFlow ──────────────────────────────────────

ConditionalNormalizingFlow::ConditionalNormalizingFlow(
    std::size_t dim, std::size_t cond_dim, std::size_t n_couplings,
    std::size_t split, std::vector<std::size_t> hidden_dims, Activation act)
    : dim_(dim), cond_dim_(cond_dim), init_split_(split) {
  if (dim_ < 2) {
    throw std::invalid_argument(
        "ConditionalNormalizingFlow: dim must be >= 2");
  }
  if (n_couplings == 0) {
    throw std::invalid_argument(
        "ConditionalNormalizingFlow: n_couplings must be >= 1");
  }
  if (split < 1 || split >= dim_) {
    throw std::invalid_argument(
        "ConditionalNormalizingFlow: split must be in [1, dim - 1]");
  }
  couplings_.reserve(n_couplings);
  swaps_.reserve(n_couplings);
  for (std::size_t i = 0; i < n_couplings; ++i) {
    couplings_.emplace_back(dim_, cond_dim_, split, hidden_dims, act);
    swaps_.push_back({dim_, split});
  }
}

void ConditionalNormalizingFlow::init_xavier(std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  for (auto& c : couplings_) c.init_xavier(rng());
}

ConditionalAffineCoupling& ConditionalNormalizingFlow::coupling(std::size_t i) {
  if (i >= couplings_.size()) {
    throw std::out_of_range(
        "ConditionalNormalizingFlow::coupling: index out of range");
  }
  return couplings_[i];
}

const ConditionalAffineCoupling&
ConditionalNormalizingFlow::coupling(std::size_t i) const {
  if (i >= couplings_.size()) {
    throw std::out_of_range(
        "ConditionalNormalizingFlow::coupling: index out of range");
  }
  return couplings_[i];
}

AffineCoupling::Output ConditionalNormalizingFlow::forward(
    const std::vector<double>& x, const std::vector<double>& cond) const {
  if (x.size() != dim_) {
    throw std::invalid_argument(
        "ConditionalNormalizingFlow::forward: input size must equal dim");
  }
  AffineCoupling::Output state;
  state.y = x;
  state.log_det_jacobian = 0.0;
  for (std::size_t i = 0; i < couplings_.size(); ++i) {
    auto step = couplings_[i].forward(state.y, cond);
    state.y = std::move(step.y);
    state.log_det_jacobian += step.log_det_jacobian;
    if (i + 1 < couplings_.size()) {
      state.y = swaps_[i].apply(state.y);
    }
  }
  return state;
}

AffineCoupling::Output ConditionalNormalizingFlow::inverse(
    const std::vector<double>& z, const std::vector<double>& cond) const {
  if (z.size() != dim_) {
    throw std::invalid_argument(
        "ConditionalNormalizingFlow::inverse: input size must equal dim");
  }
  AffineCoupling::Output state;
  state.y = z;
  state.log_det_jacobian = 0.0;
  for (std::size_t r = 0; r < couplings_.size(); ++r) {
    const std::size_t i = couplings_.size() - 1 - r;
    if (i + 1 < couplings_.size()) {
      HalfSwap inv_swap{dim_, dim_ - swaps_[i].split};
      state.y = inv_swap.apply(state.y);
    }
    auto step = couplings_[i].inverse(state.y, cond);
    state.y = std::move(step.y);
    state.log_det_jacobian += step.log_det_jacobian;
  }
  return state;
}

double ConditionalNormalizingFlow::log_density(
    const std::vector<double>& x, const std::vector<double>& cond) const {
  auto fw = forward(x, cond);
  double zz = 0.0;
  for (double v : fw.y) zz += v * v;
  constexpr double kLog2Pi = 1.8378770664093453;
  return -0.5 * zz
         - 0.5 * static_cast<double>(dim_) * kLog2Pi
         + fw.log_det_jacobian;
}

std::vector<double> ConditionalNormalizingFlow::sample(
    const std::vector<double>& cond, std::mt19937_64& rng) const {
  std::normal_distribution<double> nd(0.0, 1.0);
  std::vector<double> z(dim_);
  for (auto& v : z) v = nd(rng);
  return inverse(z, cond).y;
}

}  // namespace argus::nn
