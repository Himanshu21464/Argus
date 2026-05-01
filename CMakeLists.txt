cmake_minimum_required(VERSION 3.20)
project(argus
  VERSION 0.1.0
  DESCRIPTION "Argus — a differentiable atmospheric retrieval kernel"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "" FORCE)
endif()

add_compile_options(
  -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
)

# ─── L0 — kernel core ──────────────────────────────────────────────────
add_library(argus_core STATIC
  src/version.cpp
  src/atmosphere.cpp
  src/opacity.cpp
  src/radiative_transfer.cpp
  src/ir.cpp
  src/geometry.cpp
  src/line_list.cpp
)
target_include_directories(argus_core PUBLIC include)
target_compile_features(argus_core PUBLIC cxx_std_20)

# ─── examples ──────────────────────────────────────────────────────────
add_executable(example_01_transmission_spectrum
  examples/01_transmission_spectrum.cpp
)
target_link_libraries(example_01_transmission_spectrum PRIVATE argus_core)
set_target_properties(example_01_transmission_spectrum PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples
)

add_executable(example_02_voigt_h2o
  examples/02_voigt_h2o.cpp
)
target_link_libraries(example_02_voigt_h2o PRIVATE argus_core)
set_target_properties(example_02_voigt_h2o PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples
)

# ─── tests ─────────────────────────────────────────────────────────────
enable_testing()
add_subdirectory(tests)
