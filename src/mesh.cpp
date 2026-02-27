#include "mesh.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
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

    // Set partitions
    std::vector<std::size_t> partitions_start;
    std::vector<std::size_t> partitions_size;
    SetPartitionsPerDirection_(partitions_start, partitions_size);

    // Set number of elements (partition + ghost)
    SetNumberElements_(partitions_start, partitions_size);

    // Set elements global id
    SetElementsGlobalId_(partitions_start, partitions_size);

    // Set elements neighborhood
    SetElementsNeighborhood_();

    // Set connectivity
    SetElementsConnectivity_();

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

    // Sanity check: this block's range of element indices are less than or
    // equal to global elements per direction
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
void Mesh::SetNumberElements_(const std::vector<std::size_t> &partitions_start,
                              const std::vector<std::size_t> &partitions_size) {
    // Get rank
    const int rank = pwr::MPIUtilities::Rank();

    // Grab number of elements per direction in the partitions
    const std::size_t nx = partitions_size[3 * rank + 0];
    const std::size_t ny = partitions_size[3 * rank + 1];
    const std::size_t nz = partitions_size[3 * rank + 2];

    // Set number of elements in the partition
    num_elem_partition_ = nx * ny * nz;

    // Sanity check: no zero-element partitions allowed
    assert(num_elem_partition_ > 0);

    // Grab starting element index per direction in the partition
    const std::size_t nx_0 = partitions_start[3 * rank + 0];
    const std::size_t ny_0 = partitions_start[3 * rank + 1];
    const std::size_t nz_0 = partitions_start[3 * rank + 2];

    // Set final element index per direction in the partition
    const std::size_t nx_f = nx_0 + nx;
    const std::size_t ny_f = ny_0 + ny;
    const std::size_t nz_f = nz_0 + nz;

    // Set element neighboorhood width (depending on basis function order)
    const std::size_t w = neighborhood_width_;

    // Number of ghost element rows at the start of this block
    const std::size_t nx_ghost_0 = std::min(w, nx_0);
    const std::size_t ny_ghost_0 = std::min(w, ny_0);
    const std::size_t nz_ghost_0 = std::min(w, nz_0);

    // Number of ghost element rows at the end of this block
    const std::size_t nx_ghost_f = std::min(w, nx_ - nx_f);
    const std::size_t ny_ghost_f = std::min(w, ny_ - ny_f);
    const std::size_t nz_ghost_f = std::min(w, nz_ - nz_f);

    // Set number of elements per direction (partition + ghost)
    const std::size_t nx_total = nx + nx_ghost_0 + nx_ghost_f;
    const std::size_t ny_total = ny + ny_ghost_0 + ny_ghost_f;
    const std::size_t nz_total = nz + nz_ghost_0 + nz_ghost_f;

    // Set number of ghost elements
    num_elem_ghost_ = (nx_total * ny_total * nz_total) - (nx * ny * nz);

    // Set number of partition + ghost elements
    num_elem_total_ = num_elem_partition_ + num_elem_ghost_;

}  // Mesh::SetNumberElements_

// ----------------------------------------------------------------------------
// Set elements global id
// ----------------------------------------------------------------------------
void Mesh::SetElementsGlobalId_(
    const std::vector<std::size_t> &partitions_start,
    const std::vector<std::size_t> &partitions_size) {
    // Get rank
    const int rank = pwr::MPIUtilities::Rank();

    // Grab number of elements per direction in the partitions
    const std::size_t nx = partitions_size[3 * rank + 0];
    const std::size_t ny = partitions_size[3 * rank + 1];
    const std::size_t nz = partitions_size[3 * rank + 2];

    // Grab starting element index per direction in the partition
    const std::size_t nx_0 = partitions_start[3 * rank + 0];
    const std::size_t ny_0 = partitions_start[3 * rank + 1];
    const std::size_t nz_0 = partitions_start[3 * rank + 2];

    // Set final element index per direction in the partition
    const std::size_t nx_f = nx_0 + nx;
    const std::size_t ny_f = ny_0 + ny;
    const std::size_t nz_f = nz_0 + nz;

    // Set element neighboorhood width (depending on basis function order)
    const std::size_t w = neighborhood_width_;

    // Starting index per direction (including ghost elements)
    const std::size_t nx_start = nx_0 - std::min(w, nx_0);
    const std::size_t ny_start = ny_0 - std::min(w, ny_0);
    const std::size_t nz_start = nz_0 - std::min(w, nz_0);

    // Edning index per direction (including ghost elements)
    const std::size_t nx_end = nx_f + std::min(w, nx_ - nx_f);
    const std::size_t ny_end = ny_f + std::min(w, ny_ - ny_f);
    const std::size_t nz_end = nz_f + std::min(w, nz_ - nz_f);

    // NOTE: Global element ids should be zero-index and based on looping
    //       in the x-, y-, and z-directions (in order)

    // Resize
    elem_id_global_.resize(num_elem_partition_ + num_elem_ghost_,
                           std::numeric_limits<std::size_t>::max());

    // Indices for partition and ghost nodes
    std::size_t e_p = 0;
    std::size_t e_g = num_elem_partition_;

    // Loop partition elements
    for (std::size_t nz_i = nz_start; nz_i < nz_end; ++nz_i) {
        for (std::size_t ny_i = ny_start; ny_i < ny_end; ++ny_i) {
            for (std::size_t nx_i = nx_start; nx_i < nx_end; ++nx_i) {
                // Set global element id
                const std::size_t eid =
                    nx_i + (nx_ * ny_i) + (nx_ * ny_ * nz_i);

                // Is the current element in the partition
                const bool in_partition =
                    (nz_0 <= nz_i && nz_i < nz_f && ny_0 <= ny_i &&
                     ny_i < ny_f && nx_0 <= nx_i && nx_i < nx_f)
                        ? true
                        : false;

                // Current element is in the partition
                if (in_partition) {
                    // Sanity check: index is less than number of partition
                    // elements
                    assert(e_p < num_elem_partition_);

                    // Save global element id
                    elem_id_global_[e_p] = eid;

                    // Increment partition index
                    e_p++;

                    // Skip the rest
                    continue;
                }

                // Sanity check: index is more than number of partition elements
                assert(num_elem_partition_ <= e_g);

                // Sanity check: index is less than number of (partition +
                // ghost) elements
                assert(e_g < num_elem_total_);

                // Save global element id
                elem_id_global_[e_g] = eid;

                // Increment ghost index
                e_g++;
            }
        }
    }

    // Sanity check: looped across all partition elements
    assert(e_p == num_elem_partition_);

    // Sanity check: looped across all ghost elements
    assert(e_g == num_elem_total_);

    // Sanity check: global element ids are reasonable
    for (const auto eid : elem_id_global_) {
        assert(eid < std::numeric_limits<std::size_t>::max());
    }

}  // Mesh::SetElementsGlobalId_

// ----------------------------------------------------------------------------
// Set elements neighborhood
// ----------------------------------------------------------------------------
void Mesh::SetElementsNeighborhood_() {
    // Resize based on shape function order
    elem_neighborhood_.resize(6 *
                              num_elem_partition_);  // TODO make this "dynamic"

    // Default neighbor id is local element id
    for (std::size_t e = 0; e < num_elem_partition_; ++e) {
        for (std::size_t n = 0; n < 6; ++n) {
            elem_neighborhood_[6 * e + n] = e;
        }
    }

    // TODO -- keep going

}  // Mesh::SetElementsNeighborhood_

// ----------------------------------------------------------------------------
// Set elements connectivity
// ----------------------------------------------------------------------------
void Mesh::SetElementsConnectivity_() {
    ;  // TODO

}  // Mesh::SetElementsConnectivity_

}  // namespace pwr
