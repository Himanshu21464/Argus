#pragma once

#include <string_view>

namespace argus {

constexpr int version_major = 0;
constexpr int version_minor = 1;
constexpr int version_patch = 0;

std::string_view version_string() noexcept;

}  // namespace argus
