#include <mpi.h>

#include <cstddef>
#include <iostream>
#include <memory>

#include "mesh.h"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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

        // std::cout << "\nMESH 2" << std::endl;
        // std::shared_ptr<pwr::MeshBase> mesh1 =
        //     std::make_shared<pwr::Mesh>(x_max, y_max, z_max, nx, ny, nz);

        // std::cout << "\nMESH 3" << std::endl;
        // std::shared_ptr<pwr::MeshBase> mesh2 =
        //     std::make_shared<pwr::Mesh>(nx, ny, nz);
    }

    MPI_Finalize();

    return 0;
}
