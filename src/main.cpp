#include <mpi.h>

#include <cstddef>
#include <memory>
#include <string>

#include "mesh.h"
#include "mesh_base.h"
#include "mpi_environment.h"
#include "mpi_utilities.h"
#include "output_manager.h"
#include "vtk_writer_base.h"
#include "vtk_writer_mesh.h"

int main(int argc, char **argv) {
    // Create MPI Enviroment (RAII wrapper for MPI_Init and MPI_Finalize)
    pwr::MPIEnvironment mpi_env(argc, argv);

    // Get rank and size
    const int rank = pwr::MPIUtilities::Rank();
    const int size = pwr::MPIUtilities::Size();

    // Generate mesh
    const double x_min = -3.;
    const double y_min = -2.;
    const double z_min = -1.;
    const double x_max = 1.;
    const double y_max = 2.;
    const double z_max = 3.;
    const std::size_t nx = 4;
    const std::size_t ny = 4;
    const std::size_t nz = 4;

    const std::shared_ptr<const pwr::MeshBase> mesh =
        std::make_shared<const pwr::Mesh>(x_min, y_min, z_min, x_max, y_max,
                                          z_max, nx, ny, nz);

    // Create output manager
    const std::size_t max_step = 100;  // To testing padding

    const std::shared_ptr<pwr::OutputManager> output_manager =
        std::make_shared<pwr::OutputManager>(mesh, max_step);

    // Write mesh
    output_manager->WriteMeshFiles();

    return 0;
}
