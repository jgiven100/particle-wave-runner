# -----------------------------------------------------------------------------
# Global configuration
# -----------------------------------------------------------------------------

# Point to know Catch2 location (append list to avoid overwrites)
list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/Catch2-3.13.0/install")

# -----------------------------------------------------------------------------
# Test executable
# -----------------------------------------------------------------------------

# Create the main test executable
add_executable(particle_wave_runner_tests
    field_test.cpp
    test_main.cpp
    mesh_test.cpp
    mpi_utilities_test.cpp
    particle_poisson_test.cpp
)

# Include the header files
target_include_directories(particle_wave_runner_tests
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

# Link the executable against the core library
target_link_libraries(particle_wave_runner_tests
    PRIVATE
        pwr_core
)

# Explicitly require C++20 features for this target
target_compile_features(particle_wave_runner_tests PRIVATE cxx_std_20)

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------

# Find Catch2
find_package(Catch2 REQUIRED)

# Link test executable against Catch2
target_link_libraries(particle_wave_runner_tests
    PRIVATE
        Catch2::Catch2
)

# -----------------------------------------------------------------------------
# Register test with ctest
# -----------------------------------------------------------------------------

# Add test in order to use ctest
add_test(NAME particle_wave_runner_tests COMMAND particle_wave_runner_tests)

