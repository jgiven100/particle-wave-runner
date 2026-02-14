#include "mesh.h"

#include <cassert>
#include <iostream>
#include <limits>

namespace pwr {
// ----------------------------------------------------------------------------
// Setup mesh
// ----------------------------------------------------------------------------
void Mesh::SetupMesh_() {
    // Initialized, partition, and check mesh
    InitializeMesh_();
    PartitionMesh_();
    CheckMesh_();
}  // Mesh::SetupMesh_

// ----------------------------------------------------------------------------
// Initialize mesh
// ----------------------------------------------------------------------------
void Mesh::InitializeMesh_() {
    // Sanity check: dimenstions
    assert(x_max_ > x_min_);
    assert(y_max_ > y_min_);
    assert(z_max_ > z_min_);

    // Sanity check: number of elements
    assert(nx_ > 0);
    assert(ny_ > 0);
    assert(nz_ > 0);

    // Total elements
    num_elem_ = nx_ * ny_ * nz_;

    // Number of nodes in each direction
    const std::size_t num_nodes_x = nx_ + 1;
    const std::size_t num_nodes_y = ny_ + 1;
    const std::size_t num_nodes_z = nz_ + 1;

    // Compute element size
    dx_ = (x_max_ - x_min_) / static_cast<double>(nx_);
    dy_ = (y_max_ - y_min_) / static_cast<double>(ny_);
    dz_ = (z_max_ - z_min_) / static_cast<double>(nz_);

    // Total nodes
    const std::size_t num_nodes = num_nodes_x * num_nodes_y * num_nodes_z;

    // Generate nodal coordinates
    nodal_coords_.resize(3 * num_nodes);

    // Loop nodes in each direction
    std::size_t n_index = 0;
    for (std::size_t k = 0; k < num_nodes_z; ++k) {
        // Set z coordinate
        const double z = z_min_ + dz_ * static_cast<double>(k);

        for (std::size_t j = 0; j < num_nodes_y; ++j) {
            // Set y coordinate
            const double y = y_min_ + dy_ * static_cast<double>(j);

            for (std::size_t i = 0; i < num_nodes_x; ++i) {
                // Set x coordinate
                const double x = x_min_ + dx_ * static_cast<double>(i);

                // Nodal coordinates
                nodal_coords_[n_index * 3 + 0] = x;
                nodal_coords_[n_index * 3 + 1] = y;
                nodal_coords_[n_index * 3 + 2] = z;

                // Update node index
                n_index++;

            }  // Loop x-dir nodes

        }  // Loop y-dir nodes

    }  // Loop z-dir nodes

    for (size_t n = 0; n < num_nodes; ++n) {
        const double x = nodal_coords_[n * 3 + 0];
        const double y = nodal_coords_[n * 3 + 1];
        const double z = nodal_coords_[n * 3 + 2];

        std::cout << "x,y,z = " << x << "," << y << "," << z << std::endl;
    }

}  // Mesh::InitializeMesh_

// ----------------------------------------------------------------------------
// Partition mesh
// ----------------------------------------------------------------------------
void Mesh::PartitionMesh_() {
    // TODO
}  // Mesh::PartitionMesh_

// ----------------------------------------------------------------------------
// Check mesh
// ----------------------------------------------------------------------------
void Mesh::CheckMesh_() {
    // Setup is complete
    setup_complete_ = true;
}  // Mesh::CheckMesh_

}  // namespace pwr
