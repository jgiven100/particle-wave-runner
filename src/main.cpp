#include <mpi.h>

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "field.h"
#include "mesh.h"
#include "mesh_base.h"
#include "mpi_environment.h"
#include "mpi_utilities.h"
#include "particle_base.h"
#include "particle_poisson.h"
#include "solver_base.h"
#include "solver_explicit.h"
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
    const double z_max = 10.;
    const std::size_t nx = 1;
    const std::size_t ny = 1;
    const std::size_t nz = 10;

    const std::shared_ptr<const pwr::MeshBase> mesh =
        std::make_shared<const pwr::Mesh>(x_min, y_min, z_min, x_max, y_max,
                                          z_max, nx, ny, nz);

    // ------------------------------------------------------------------------
    // Generate particles
    // ------------------------------------------------------------------------
    const double nz_double = static_cast<double>(nz);
    const double ny_double = static_cast<double>(ny);
    const double nx_double = static_cast<double>(nx);

    std::vector<std::shared_ptr<pwr::ParticleBase>> particles;
    std::size_t id = 0;

    for (std::size_t nz_i = 0; nz_i < nz; ++nz_i) {
        // Set coordinate z-direction
        const double z_i =
            (0.5 + static_cast<double>(nz_i)) * (z_max - z_min) / nz_double;
        for (std::size_t ny_i = 0; ny_i < ny; ++ny_i) {
            // Set coordinate y-direction
            const double y_i =
                (0.5 + static_cast<double>(ny_i)) * (y_max - y_min) / ny_double;
            for (std::size_t nx_i = 0; nx_i < nx; ++nx_i) {
                // Set coordinate x-direction
                const double x_i = (0.5 + static_cast<double>(nx_i)) *
                                   (x_max - x_min) / nx_double;

                // Set coords array
                std::array<double, 3> coords{x_i, y_i, z_i};

                // Generate particle
                const std::shared_ptr<pwr::ParticleBase> particle =
                    std::make_shared<pwr::ParticlePoisson>(id, coords);

                // Save it to vector
                particles.emplace_back(particle);

                // Update id
                id++;
            }
        }
    }

    // ------------------------------------------------------------------------
    // Set fields
    // ------------------------------------------------------------------------
    const std::shared_ptr<pwr::Field> field_u =
        std::make_shared<pwr::Field>(mesh->GetNumNodesActive());

    for (std::size_t i = 0; i < field_u->GetFieldSize(); ++i) {
        (*field_u)[i] = 0;
    }

    const std::vector<std::shared_ptr<pwr::Field>> fields{field_u};
    const std::vector<std::string> fields_names{"temperature"};

    // ------------------------------------------------------------------------
    // Set solver
    // ------------------------------------------------------------------------
    const std::shared_ptr<pwr::SolverBase> solver_explicit =
        std::make_shared<pwr::SolverExplicit>(mesh, particles, fields,
                                              fields_names);

    // ------------------------------------------------------------------------
    // Run
    // ------------------------------------------------------------------------
    solver_explicit->Run();

    // ------------------------------------------------------------------------
    // Done
    // ------------------------------------------------------------------------
    return 0;
}
