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

    // Set x, y, and z limits
    const double x_min = 0.;
    const double y_min = 0.;
    const double z_min = 0.;
    const double x_max = 1.;
    const double y_max = 1.;
    const double z_max = 1.;

    // Number of elements in each direction
    const std::size_t nx = 4;
    const std::size_t ny = 4;
    const std::size_t nz = 4;

    // Create mesh object for *this* partition
    const std::shared_ptr<const pwr::MeshBase> mesh =
        std::make_shared<const pwr::Mesh>(x_min, y_min, z_min, x_max, y_max,
                                          z_max, nx, ny, nz);

    // ------------------------------------------------------------------------
    // Generate particles
    // ------------------------------------------------------------------------

    // Particles per direction per cell
    const std::size_t n = 4;
    const double n_double = static_cast<double>(n);

    // Particles per element
    const std::size_t num_particles_elem = n * n * n;

    // Set element-wise (local) particle indices (per direction)
    std::vector<std::size_t> particle_index_elem(3 * num_particles_elem, 0);

    // Loop particles per direction
    std::size_t p_local = 0;
    for (std::size_t p_z_i = 0; p_z_i < n; ++p_z_i) {
        for (std::size_t p_y_i = 0; p_y_i < n; ++p_y_i) {
            for (std::size_t p_x_i = 0; p_x_i < n; ++p_x_i) {
                // Save particle index per direction
                particle_index_elem[3 * p_local + 0] = p_x_i;
                particle_index_elem[3 * p_local + 1] = p_y_i;
                particle_index_elem[3 * p_local + 2] = p_z_i;

                // Increment
                p_local++;
            }
        }
    }

    // Compute element size
    const double e_dx = (x_max - x_min) / static_cast<double>(nx);
    const double e_dy = (y_max - y_min) / static_cast<double>(ny);
    const double e_dz = (z_max - z_min) / static_cast<double>(nz);

    // Compute particle size
    const double p_dx = e_dx / n_double;
    const double p_dy = e_dy / n_double;
    const double p_dz = e_dz / n_double;

    // Grab number of partition elements
    const std::size_t num_elem_partition = mesh->GetNumElemPartition();

    // Grab number of total (partition + ghost) elements
    const std::size_t num_elem_total = mesh->GetNumElemTotal();

    // Grab global element id
    const std::vector<std::size_t> elem_id_global = mesh->GetElemIdGlobal();

    // Grab global element indices (in each direction)
    const std::vector<std::size_t> elem_index_global =
        mesh->GetElemIndexGlobal();

    // Place to save list of particle objects
    std::vector<std::shared_ptr<pwr::ParticleBase>> particles;
    particles.reserve(num_elem_total * num_particles_elem);

    // Loop total (partition + ghost) elements
    for (std::size_t e = 0; e < num_elem_total; ++e) {
        // Set ownership to `false` if particle is in ghost element
        const bool owned = (e >= num_elem_partition) ? false : true;

        // Set global element id
        const std::size_t eid_global = elem_id_global[e];

        // Set global element indices
        const std::size_t e_x_i = elem_index_global[3 * e + 0];
        const std::size_t e_y_i = elem_index_global[3 * e + 1];
        const std::size_t e_z_i = elem_index_global[3 * e + 2];

        // Compute element corner
        const double e_x0 = x_min + e_dx * static_cast<double>(e_x_i);
        const double e_y0 = y_min + e_dy * static_cast<double>(e_y_i);
        const double e_z0 = z_min + e_dz * static_cast<double>(e_z_i);

        // Loop particles per element
        for (std::size_t p = 0; p < num_particles_elem; ++p) {
            // Set local particle id      TODO -- local here means
            // partition-wise
            std::size_t pid_local = 0;  // TODO -- does this matter?

            // Set global particle id
            std::size_t pid_global = (eid_global * num_particles_elem + p);

            // // Decompose into local indices in each direction
            // const std::size_t p_x_i = p % n;
            // const std::size_t p_y_i = (p / n) % n;
            // const std::size_t p_z_i = p / (n * n);

            // Grab element-wise (local) indices in each direction
            const std::size_t p_x_i = particle_index_elem[3 * p + 0];
            const std::size_t p_y_i = particle_index_elem[3 * p + 1];
            const std::size_t p_z_i = particle_index_elem[3 * p + 2];

            // Set coords array
            const double x = e_x0 + (static_cast<double>(p_x_i) + 0.5) * p_dx;
            const double y = e_y0 + (static_cast<double>(p_y_i) + 0.5) * p_dy;
            const double z = e_z0 + (static_cast<double>(p_z_i) + 0.5) * p_dz;
            std::array<double, 3> coords{x, y, z};

            // Generate particle
            const std::shared_ptr<pwr::ParticleBase> particle =
                std::make_shared<pwr::ParticlePoisson>(pid_local, pid_global,
                                                       owned, coords);

            // Save it to vector
            particles.emplace_back(particle);
        }
    }

    //
    const std::string particles_name = "temperature";

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
        std::make_shared<pwr::SolverExplicit>(mesh, particles, particles_name,
                                              fields, fields_names);

    // ------------------------------------------------------------------------
    // Run
    // ------------------------------------------------------------------------
    solver_explicit->Run();

    // ------------------------------------------------------------------------
    // Done
    // ------------------------------------------------------------------------
    return 0;
}
