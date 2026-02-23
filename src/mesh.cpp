#include "mesh.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>

#include "mpi_utilities.h"

namespace pwr {
// ----------------------------------------------------------------------------
// Setup mesh
// ----------------------------------------------------------------------------
void Mesh::SetupMesh_() {
    // Initialize, partition, and check mesh
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

    // Set total number of elements (all partitions)
    num_elem_ = nx_ * ny_ * nz_;

    // Compute element size
    dx_ = (x_max_ - x_min_) / static_cast<double>(nx_);
    dy_ = (y_max_ - y_min_) / static_cast<double>(ny_);
    dz_ = (z_max_ - z_min_) / static_cast<double>(nz_);

}  // Mesh::InitializeMesh_

// ----------------------------------------------------------------------------
// Partition mesh
// ----------------------------------------------------------------------------
void Mesh::PartitionMesh_() {
    // Get rank and size
    const int rank = pwr::MPIUtilities::Rank();
    const int size = pwr::MPIUtilities::Size();

    // Set partition
    std::vector<std::size_t> partitions_start;
    std::vector<std::size_t> partitions_size;
    SetPartitionsPerDirection_(partitions_start, partitions_size);

    // Compute rank coordinates

    // Compute owned element ranges

    // Compute ghost element ranges

    // Build connectivity

    // Set number of elements
    SetNumberElements_(partitions_size);

}  // Mesh::PartitionMesh_

// ----------------------------------------------------------------------------
// Check mesh
// ----------------------------------------------------------------------------
void Mesh::CheckMesh_() {
    // Setup is complete
    setup_complete_ = true;

}  // Mesh::CheckMesh_

// ----------------------------------------------------------------------------
// Set partitions per direction
// ----------------------------------------------------------------------------
void Mesh::SetPartitionsPerDirection_(
    std::vector<std::size_t> &partitions_start,
    std::vector<std::size_t> &partitions_size) {
    // Get size
    const int size = pwr::MPIUtilities::Size();
    const int rank = pwr::MPIUtilities::Rank();

    // Resize and set zero
    partitions_start.resize(3 * size, 0);
    partitions_size.resize(3 * size, 0);

    // Start with one partition per direction
    partitions_size[0] = nx_;
    partitions_size[1] = ny_;
    partitions_size[2] = nz_;

    // Next unassigned index
    std::size_t n_index = 1;

    // Split each direction sequentially
    while (n_index < size) {
        // Find block with the largest number of elements
        std::size_t b_index = 0;
        std::size_t max_num_elements = 0;

        // Loop each block
        for (size_t i = 0; i < n_index; ++i) {
            // Set the number of elements per direction in the current block
            const std::size_t nx_i = partitions_size[3 * i + 0];
            const std::size_t ny_i = partitions_size[3 * i + 1];
            const std::size_t nz_i = partitions_size[3 * i + 2];

            // Total elements in the current block
            const std::size_t num_elements_i = nx_i * ny_i * nz_i;

            // Update the index if the largest block is found
            if (num_elements_i > max_num_elements) {
                max_num_elements = num_elements_i;
                b_index = i;
            }
        }

        // Set elements in each direction
        const std::size_t nx = partitions_size[3 * b_index + 0];
        const std::size_t ny = partitions_size[3 * b_index + 1];
        const std::size_t nz = partitions_size[3 * b_index + 2];

        // Split in the x-direction
        if (nx >= ny && nx >= nz && nx > 1) {
            // Split into "left" and "right" blocks
            const std::size_t left = nx / 2;
            const std::size_t right = nx - left;

            // Sanity check: successful x-direction split
            assert(nx == left + right);

            // Update existing block in place with "left"
            partitions_size[3 * b_index + 0] = left;

            // Assign new block in next open position with "right"
            partitions_size[3 * n_index + 0] = right;
            partitions_size[3 * n_index + 1] = partitions_size[3 * b_index + 1];
            partitions_size[3 * n_index + 2] = partitions_size[3 * b_index + 2];

            // Assign starting index for new block
            partitions_start[3 * n_index + 0] =
                partitions_start[3 * b_index + 0] + left;
            partitions_start[3 * n_index + 1] =
                partitions_start[3 * b_index + 1];
            partitions_start[3 * n_index + 2] =
                partitions_start[3 * b_index + 2];
        }

        // Split in the y-direction
        else if (ny >= nx && ny >= nz && ny > 1) {
            // Split into "left" and "right" blocks
            const std::size_t left = ny / 2;
            const std::size_t right = ny - left;

            // Sanity check: successful y-direction split
            assert(ny == left + right);

            // Update existing block in place with "left"
            partitions_size[3 * b_index + 1] = left;

            // Assign new block in next open position with "right"
            partitions_size[3 * n_index + 0] = partitions_size[3 * b_index + 0];
            partitions_size[3 * n_index + 1] = right;
            partitions_size[3 * n_index + 2] = partitions_size[3 * b_index + 2];

            // Assign starting index for new block
            partitions_start[3 * n_index + 0] =
                partitions_start[3 * b_index + 0];
            partitions_start[3 * n_index + 1] =
                partitions_start[3 * b_index + 1] + left;
            partitions_start[3 * n_index + 2] =
                partitions_start[3 * b_index + 2];
        }

        // Split in the z-direction
        else if (nz > 1) {
            // Split into "left" and "right" blocks
            const std::size_t left = nz / 2;
            const std::size_t right = nz - left;

            // Sanity check: successful z-direction split
            assert(nz == left + right);

            // Update existing block in place with "left"
            partitions_size[3 * b_index + 2] = left;

            // Assign new block in next open position with "right"
            partitions_size[3 * n_index + 0] = partitions_size[3 * b_index + 0];
            partitions_size[3 * n_index + 1] = partitions_size[3 * b_index + 1];
            partitions_size[3 * n_index + 2] = right;

            // Assign starting index for new block
            partitions_start[3 * n_index + 0] =
                partitions_start[3 * b_index + 0];
            partitions_start[3 * n_index + 1] =
                partitions_start[3 * b_index + 1];
            partitions_start[3 * n_index + 2] =
                partitions_start[3 * b_index + 2] + left;
        }

        // Cannot split further
        else {
            break;
        }

        // Update index for next unassigned block
        n_index++;
    }

    // Sanity check: sum of partition elements equals total elements
    std::size_t num_elem = 0;
    for (std::size_t i = 0; i < size; ++i) {
        num_elem += partitions_size[3 * i + 0] * partitions_size[3 * i + 1] *
                    partitions_size[3 * i + 2];
    }
    assert(num_elem == num_elem_);

    // Sanity check: this block's range of element indices are less than global
    // elements per direction
    assert(partitions_start[3 * rank + 0] + partitions_size[3 * rank + 0] <=
           nx_);
    assert(partitions_start[3 * rank + 1] + partitions_size[3 * rank + 1] <=
           ny_);
    assert(partitions_start[3 * rank + 2] + partitions_size[3 * rank + 2] <=
           nz_);

}  // Mesh::SetPartitionsPerDirection_

// ----------------------------------------------------------------------------
// Set number of elements
// ----------------------------------------------------------------------------
void Mesh::SetNumberElements_(const std::vector<std::size_t> &partitions_size) {
    // Get rank
    const int rank = pwr::MPIUtilities::Rank();

    // Set number of elements partitions
    num_elem_partition_ = partitions_size[3 * rank + 0] *
                          partitions_size[3 * rank + 1] *
                          partitions_size[3 * rank + 2];

    // Sanity check: no zero-element partitions allowed
    assert(num_elem_partition_ > 0);

    // Set number of elements ghost
    num_elem_ghost_ = 0;  // TODO update

}  // Mesh::SetNumberElements_

}  // namespace pwr
