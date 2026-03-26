# -----------------------------------------------------------------------------
# Test executable
# -----------------------------------------------------------------------------

# Create the main test executable
add_executable(particle_wave_runner_tests
    test_main.cpp
    mesh_test.cpp

    # Add Catch2 implementation
    ${PROJECT_SOURCE_DIR}/external/catch2/catch_amalgamated.cpp
)

# Include the header files
target_include_directories(particle_wave_runner_tests
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/external/catch2
)

# Link the executable against the core library
target_link_libraries(particle_wave_runner_tests
    PRIVATE
        pwr_core
)

# -----------------------------------------------------------------------------
# Register test with ctest
# -----------------------------------------------------------------------------

# Add test in order to use ctest
add_test(NAME particle_wave_runner_tests COMMAND particle_wave_runner_tests)

