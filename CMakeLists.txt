cmake_minimum_required(VERSION 3.20)
project(particle_wave_runner LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(PWR_BUILD_TESTS "Build tests" ON)
option(PWR_BUILD_DOCS "Build documentation" ON)

# MPI
find_package(MPI REQUIRED)

# Libarary

# Executable
add_executable(particle_wave_runner src/main.cpp)
target_link_libraries(particle_wave_runner PRIVATE MPI::MPI_CXX)

# Tests
if(PWR_BUILD_TESTS)
    # TODO
endif()

# Docs
if(PWR_BUILD_DOCS)
    # TODO
endif()

