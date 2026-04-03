#define CATCH_CONFIG_RUNNER

#include <catch2/catch_all.hpp>
#include <sstream>
#include <string>

#include "mpi_environment.h"
#include "mpi_utilities.h"
#include "utilities.h"

int main(int argc, char **argv) {
    // Create MPI Environment (RAII wrapper for MPI_Init and MPI_Finalize)
    pwr::MPIEnvironment mpi_env(argc, argv);

    // Number of MPI processes must be less than 5
    if (pwr::MPIUtilities::Size() > 4) {
        std::stringstream oss;
        oss << "Too many MPI processes (>4) for `./particle_wave_runner_test`.";
        pwr::Utilities::PrintErrorOnRoot(oss.str());

        return -1;
    }

    // Create Catch2 session
    Catch::Session session;

    // Quite excpet for root (rank 0)
    if (pwr::MPIUtilities::Rank() != 0) {
        freopen("/dev/null", "w", stdout);
        freopen("/dev/null", "w", stderr);
    }

    // Run Catch2 session
    int result = session.run(argc, argv);

    return result;
}
