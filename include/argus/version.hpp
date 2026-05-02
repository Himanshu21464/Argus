#pragma once

#include <string_view>

namespace argus {

constexpr int version_major = 0;
constexpr int version_minor = 7;
constexpr int version_patch = 5;

std::string_view version_string() noexcept;

}  // namespace argus
