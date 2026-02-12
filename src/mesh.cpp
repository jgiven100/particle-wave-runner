#include "mesh.h"

#include <cassert>

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

    // Compute element size
    dx_ = (x_max_ - x_min_) / static_cast<double>(nx_);
    dy_ = (y_max_ - y_min_) / static_cast<double>(ny_);
    dz_ = (z_max_ - z_min_) / static_cast<double>(nz_);

    // Total elements
    num_elem_ = nx_ * ny_ * nz_;

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
    // TODO
    setup_complete_ = true;
}  // Mesh::CheckMesh_

}  // namespace pwr
