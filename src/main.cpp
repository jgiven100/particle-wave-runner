#include <mpi.h>

#include <cstddef>
#include <iostream>
#include <memory>

#include "mesh.h"
#include "mesh_base.h"
#include "mpi_environment.h"
#include "mpi_utilities.h"

int main(int argc, char **argv) {
    // Create MPI Enviroment (RAII wrapper for MPI_Init and MPI_Finalize)
    pwr::MPIEnvironment mpi_env(argc, argv);

    // Get MPI rank
    const int rank = pwr::MPIUtilities::Rank();

    // Get MPI size
    const int size = pwr::MPIUtilities::Size();

    std::cout << "this is rank " << rank << " out of " << size << std::endl;

    if (rank == 0) {
        const double x_min = -3.;
        const double y_min = -2.;
        const double z_min = -1.;
        const double x_max = 1.;
        const double y_max = 2.;
        const double z_max = 3.;
        const std::size_t nx = 4;
        const std::size_t ny = 4;
        const std::size_t nz = 4;

        std::cout << "MESH 1" << std::endl;
        std::shared_ptr<pwr::MeshBase> mesh0 = std::make_shared<pwr::Mesh>(
            x_min, y_min, z_min, x_max, y_max, z_max, nx, ny, nz);
    }

    return 0;
}
