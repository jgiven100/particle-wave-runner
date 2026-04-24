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
#include "particle_base.h"
#include "particle_poisson.h"
#include "solver_base.h"
#include "solver_explicit.h"
#include "solver_implicit.h"
#include "utilities.h"
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
    const std::size_t nx = 16;
    const std::size_t ny = 16;
    const std::size_t nz = 16;

    const std::shared_ptr<const pwr::MeshBase> mesh =
        std::make_shared<const pwr::Mesh>(x_min, y_min, z_min, x_max, y_max,
                                          z_max, nx, ny, nz);

    // ------------------------------------------------------------------------
    // Set particles
    // ------------------------------------------------------------------------
    const std::shared_ptr<pwr::ParticleBase> particle =
        std::make_shared<pwr::ParticlePoisson>();

    // ------------------------------------------------------------------------
    // Set fields
    // ------------------------------------------------------------------------
    const std::shared_ptr<pwr::Field> field_rank =
        std::make_shared<pwr::Field>(mesh->GetNumNodesActive());

    const std::shared_ptr<pwr::Field> field_owned =
        std::make_shared<pwr::Field>(mesh->GetNumNodesActive());
    const auto &nodal_ownership = mesh->GetNodalOwnership();

    for (std::size_t i = 0; i < field_rank->GetFieldSize(); ++i) {
        (*field_rank)[i] = rank;
        (*field_owned)[i] = nodal_ownership[i];
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

    output_manager->AddField(field_rank, "rank");
    output_manager->AddField(field_owned, "owned");

    output_manager->WriteMeshFiles(step);

    output_manager->WriteFieldFiles(step);

    // ------------------------------------------------------------------------
    // Done
    // ------------------------------------------------------------------------
    return 0;
}
