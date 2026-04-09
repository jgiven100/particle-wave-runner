# Minimum CMake version required for project
cmake_minimum_required(VERSION 3.20)

# Project name and language
project(particle_wave_runner LANGUAGES CXX)

# -----------------------------------------------------------------------------
# Global configuration
# -----------------------------------------------------------------------------

# Use the C++20 standard
set(CMAKE_CXX_STANDARD 20)

# Require the specified C++20 standard (do not silently downgrade)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Disable compiler-specific extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# Point to known VTK location (append list to avoid overwrites)
list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/vtk/install")

# -----------------------------------------------------------------------------
# Default build type
# -----------------------------------------------------------------------------

# Set default build type
#
# Debug          -> -O0 -g
# RelWithDebInfo -> -O2 -g
# Release        -> -O3 -DNDEBUG
#
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# -----------------------------------------------------------------------------
# Build options
# -----------------------------------------------------------------------------

# Option to build unit tests
option(PWR_BUILD_TESTS "Build tests" ON)

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------

# Find MPI (required for parallel execution)
find_package(MPI REQUIRED)

# Find VTK (required for writing output files)
#
# CommonCore      -> basic VTK infrastructure
# CommonDataModel -> data structures
# IOXML           -> writers
#
find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    IOXML
)

# VTK is quiet by default, so pass a message to show the package is found
message(STATUS "VTK version: ${VTK_VERSION}")

# -----------------------------------------------------------------------------
# Core library
# -----------------------------------------------------------------------------

# Create the core project library
add_library(pwr_core
    src/mesh.cpp
    src/output_manager.cpp
    src/solver_explicit.cpp
    src/solver_implicit.cpp
    src/vtk_writer_mesh.cpp
)

# Set compiler options (for gcc)
#
# -Wall      -> common warnings
# -Wextra    -> additional useful warnings
# -Wpedantic -> enforce strict ISO compliance
#
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(pwr_core
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
    )
endif()

# Public include directories
target_include_directories(pwr_core
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

# Link required external libraries
target_link_libraries(pwr_core
    PUBLIC
        MPI::MPI_CXX

        VTK::CommonCore
        VTK::CommonDataModel
        VTK::IOXML
)

# Explicitly require C++20 features for this target
target_compile_features(pwr_core PUBLIC cxx_std_20)

# -----------------------------------------------------------------------------
# Executable
# -----------------------------------------------------------------------------

# Create the main application executable
add_executable(particle_wave_runner
    src/main.cpp
)

# Link the executable against the core library
target_link_libraries(particle_wave_runner 
    PRIVATE
        pwr_core
)

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

# Check if build tests option is "on"
if (PWR_BUILD_TESTS)
    # Enable testing
    enable_testing()

    # Add tests/ directory
    add_subdirectory(tests/)
endif()

