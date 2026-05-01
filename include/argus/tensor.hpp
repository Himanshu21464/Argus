#pragma once

#include <cstddef>
#include <initializer_list>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace argus {

// Minimal CPU-resident dense tensor. The production kernel will replace this
// with a device-aware tensor type (CUDA / ROCm) backing the same shape API.
// The intent is that callers depend only on this header — implementation
// switches do not propagate.
class Tensor {
 public:
  Tensor() = default;

  explicit Tensor(std::vector<std::size_t> shape)
      : shape_(std::move(shape)), data_(num_elements(shape_), 0.0) {}

  Tensor(std::vector<std::size_t> shape, std::vector<double> data)
      : shape_(std::move(shape)), data_(std::move(data)) {
    if (data_.size() != num_elements(shape_)) {
      throw std::invalid_argument("Tensor: shape/data size mismatch");
    }
  }

  const std::vector<std::size_t>& shape() const noexcept { return shape_; }
  std::size_t size() const noexcept { return data_.size(); }
  std::size_t rank() const noexcept { return shape_.size(); }

  double& operator[](std::size_t i) { return data_[i]; }
  double operator[](std::size_t i) const { return data_[i]; }

  double* data() noexcept { return data_.data(); }
  const double* data() const noexcept { return data_.data(); }

  // 2-D index helper for the common (layer, species) case.
  double& at(std::size_t i, std::size_t j) {
    return data_[i * shape_[1] + j];
  }
  double at(std::size_t i, std::size_t j) const {
    return data_[i * shape_[1] + j];
  }

 private:
  static std::size_t num_elements(const std::vector<std::size_t>& s) {
    return std::accumulate(s.begin(), s.end(), std::size_t{1},
                           std::multiplies<>{});
  }

  std::vector<std::size_t> shape_;
  std::vector<double> data_;
};

}  // namespace argus
