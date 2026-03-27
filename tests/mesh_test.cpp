#include "mesh.h"

#include <catch2/catch_all.hpp>
#include <cstddef>
#include <memory>
#include <vector>

#include "mesh_base.h"

TEST_CASE("MESH") {
    // ------------------------------------------------------------------------
    // Create mesh objects
    // ------------------------------------------------------------------------

    // Place to save meshes created using constructor
    std::vector<std::shared_ptr<const pwr::MeshBase>> meshes(3);

    // Global coordinates
    const double x_min = 0.;
    const double y_min = 0.;
    const double z_min = 0.;
    const double x_max = 1.;
    const double y_max = 1.;
    const double z_max = 1.;

    // Number of elements per direction
    const std::size_t nx = 2;
    const std::size_t ny = 2;
    const std::size_t nz = 2;

    // Create shared pointer to mesh object with full constructor
    meshes[0] = std::make_shared<const pwr::Mesh>(x_min, y_min, z_min, x_max,
                                                  y_max, z_max, nx, ny, nz);

    // Create shared pointer to mesh object with constructor delegation
    meshes[1] =
        std::make_shared<const pwr::Mesh>(x_max, y_max, z_max, nx, ny, nz);

    // Create shared pointer to mesh object with constructor delegation
    meshes[2] = std::make_shared<const pwr::Mesh>(nx, ny, nz);

    // ------------------------------------------------------------------------
    // Test number of partition elements
    // ------------------------------------------------------------------------
    const std::size_t num_elem = 8;

    for (const auto &mesh : meshes) {
        REQUIRE(mesh->GetNumElemPartition() == num_elem);
    }

    // Placeholder
    REQUIRE(true);
}
