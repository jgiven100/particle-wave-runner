#include <mpi.h>

#include <cstddef>
#include <memory>
#include <string>

#include "field.h"
#include "mesh.h"
#include "mesh_base.h"
#include "mpi_environment.h"
#include "mpi_utilities.h"
#include "output_manager.h"
#include "solver_base.h"
#include "solver_explicit.h"
#include "solver_implicit.h"
#include "vtk_writer_base.h"
#include "vtk_writer_mesh.h"

int main(int argc, char **argv) {
    // Create MPI Enviroment (RAII wrapper for MPI_Init and MPI_Finalize)
    pwr::MPIEnvironment mpi_env(argc, argv);

    // Get rank and size
    const int rank = pwr::MPIUtilities::Rank();
    const int size = pwr::MPIUtilities::Size();

    // ------------------------------------------------------------------------
    // Generate mesh
    // ------------------------------------------------------------------------
    const double x_min = 0.;
    const double y_min = 0.;
    const double z_min = 0.;
    const double x_max = 1.;
    const double y_max = 1.;
    const double z_max = 1.;
    const std::size_t nx = 3;
    const std::size_t ny = 3;
    const std::size_t nz = 3;

    const std::shared_ptr<const pwr::MeshBase> mesh =
        std::make_shared<const pwr::Mesh>(x_min, y_min, z_min, x_max, y_max,
                                          z_max, nx, ny, nz);

    // ------------------------------------------------------------------------
    // Set field
    // ------------------------------------------------------------------------
    const std::shared_ptr<pwr::Field> field =
        std::make_shared<pwr::Field>(mesh->GetNumNodesActive());

    for (std::size_t i = 0; i < field->GetFieldSize(); ++i) {
        (*field)[i] = rank;
    }

    // ------------------------------------------------------------------------
    // Set operator (physics engine)
    // ------------------------------------------------------------------------

    // ------------------------------------------------------------------------
    // Set solver
    // ------------------------------------------------------------------------

    const std::shared_ptr<pwr::SolverBase> solver_explicit =
        std::make_shared<pwr::SolverExplicit>();

    const std::shared_ptr<pwr::SolverBase> solver_implicit =
        std::make_shared<pwr::SolverImplicit>();

    // ------------------------------------------------------------------------
    // Run
    // ------------------------------------------------------------------------

    solver_explicit->Step();

    solver_implicit->Step();

    // ------------------------------------------------------------------------
    // Create output manager
    // ------------------------------------------------------------------------
    const std::size_t step = 0;
    const std::size_t max_step = 100;  // To testing padding

    const std::shared_ptr<pwr::OutputManager> output_manager =
        std::make_shared<pwr::OutputManager>(mesh, max_step);

    output_manager->AddField(field, "Field_0");
    output_manager->AddField(field, "Field_1");

    output_manager->WriteMeshFiles(step);

    output_manager->WriteFieldFiles(step);

    // ------------------------------------------------------------------------
    // Done
    // ------------------------------------------------------------------------
    return 0;
}
