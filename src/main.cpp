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
        const double x_min = -10.;
        const double y_min = -9.;
        const double z_min = -8.;
        const double x_max = 8.;
        const double y_max = 9.;
        const double z_max = 10.;
        const std::size_t nx = 100;
        const std::size_t ny = 200;
        const std::size_t nz = 300;

        std::cout << "MESH 1" << std::endl;
        std::shared_ptr<pwr::Mesh> mesh0 = std::make_shared<pwr::Mesh>(
            x_min, y_min, z_min, x_max, y_max, z_max, nx, ny, nz);

        std::cout << "\nMESH 2" << std::endl;
        std::shared_ptr<pwr::Mesh> mesh1 =
            std::make_shared<pwr::Mesh>(x_max, y_max, z_max, nx, ny, nz);

        std::cout << "\nMESH 3" << std::endl;
        std::shared_ptr<pwr::Mesh> mesh2 =
            std::make_shared<pwr::Mesh>(nx, ny, nz);
    }

    MPI_Finalize();

    return 0;
}
