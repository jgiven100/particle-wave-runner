#include "field.h"

#include <catch2/catch_all.hpp>
#include <cstddef>
#include <memory>
#include <vector>

#include "mesh.h"
#include "mpi_utilities.h"

TEST_CASE("Field", "[field]") {
    // Grab rank and size
    const int rank = pwr::MPIUtilities::Rank();
    const int size = pwr::MPIUtilities::Size();

    // Number of elements per direction
    const std::size_t nx = 3;
    const std::size_t ny = 3;
    const std::size_t nz = 3;

    // Create shared pointer to mesh object
    const auto mesh = std::make_shared<const pwr::Mesh>(nx, ny, nz);

    // Get the number of active nodes for this partition
    const std::size_t num_nodes_active = mesh->GetNumNodesActive();

    // Get nodal ownership for this partition
    const auto& nodal_ownership = mesh->GetNodalOwnership();

    // Create shared point to field object
    const auto field = std::make_shared<pwr::Field>(num_nodes_active);

    // Set field values
    for (std::size_t i = 0; i < num_nodes_active; ++i) {
        (*field)[i] = static_cast<double>(nodal_ownership[i]);
    }

    // ------------------------------------------------------------------------
    // Field
    // ------------------------------------------------------------------------
    {
        INFO("Field");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Grab field
        const auto& f = field->GetField();

        // Loop each value in field
        for (std::size_t i = 0; i < f.size(); ++i) {
            // -- Check if getter matches nodal ownership from mesh -----------
            if (f[i] != static_cast<double>(nodal_ownership[i])) {
                num_failed_tests_local++;
            }
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root (rank 0)
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // Field

    // ------------------------------------------------------------------------
    // Field size
    // ------------------------------------------------------------------------
    {
        INFO("Field size");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // -- Check if getter matches number of nodes from mesh ---------------
        if (field->GetFieldSize() != num_nodes_active) {
            num_failed_tests_local++;
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root (rank 0)
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // Field size

}  // TEST_CASE("Field")
