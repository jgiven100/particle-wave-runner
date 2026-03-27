#define CATCH_CONFIG_RUNNER

#include <catch2/catch_all.hpp>

#include "mpi_environment.h"

int main(int argc, char **argv) {
    // Create MPI Environment (RAII wrapper for MPI_Init and MPI_Finalize)
    pwr::MPIEnvironment mpi_env(argc, argv);

    // Create Catch2 session
    Catch::Session session;

    return session.run(argc, argv);
}
