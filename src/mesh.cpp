#include "mesh.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>

#include "mpi_utilities.h"

namespace pwr {
// ----------------------------------------------------------------------------
// Setup mesh
// ----------------------------------------------------------------------------
void Mesh::SetupMesh_() {
    InitializeMesh_();
    PartitionMesh_();
    NumberMesh_();
    ConnectMesh_();
    CheckMesh_();

}  // Mesh::SetupMesh_

// ----------------------------------------------------------------------------
// Initialize mesh
// ----------------------------------------------------------------------------
void Mesh::InitializeMesh_() {
    // Sanity check: number of elements
    assert(nx_ > 0);
    assert(ny_ > 0);
    assert(nz_ > 0);

    // Set total number of elements (all partitions)
    num_elem_ = nx_ * ny_ * nz_;

    // Set number of nodes
    nodes_x_ = nx_ + 1;
    nodes_y_ = ny_ + 1;
    nodes_z_ = nz_ + 1;

    // Set total number of nodes (all partitions)
    num_nodes_ = nodes_x_ * nodes_y_ * nodes_z_;

    // Sanity check: more nodes than elements
    assert(num_nodes_ > num_elem_);

    // Sanity check: dimensions
    assert(x_max_ > x_min_);
    assert(y_max_ > y_min_);
    assert(z_max_ > z_min_);

    // Compute element size
    dx_ = (x_max_ - x_min_) / static_cast<double>(nx_);
    dy_ = (y_max_ - y_min_) / static_cast<double>(ny_);
    dz_ = (z_max_ - z_min_) / static_cast<double>(nz_);

}  // Mesh::InitializeMesh_

// ----------------------------------------------------------------------------
// Partition mesh
// ----------------------------------------------------------------------------
void Mesh::PartitionMesh_() {
    // Set partitions
    std::vector<std::size_t> partitions_start;
    std::vector<std::size_t> partitions_size;
    SetPartitionsPerDirection_(partitions_start, partitions_size);

    // Set number of elements (partition + ghost) for this partition
    SetElementsNumbering_(partitions_start, partitions_size);

}  // Mesh::PartitionMesh_

// ----------------------------------------------------------------------------
// Number mesh
// ----------------------------------------------------------------------------
void Mesh::NumberMesh_() {
    // Set elements global id
    SetElementsGlobalId_();

    // Set elements neighborhood
    SetElementsNeighborhood_();

}  // Mesh::NumberMesh_

// ----------------------------------------------------------------------------
// Connect mesh
// ----------------------------------------------------------------------------
void Mesh::ConnectMesh_() {
    // Set connectivity
    SetElementsConnectivity_();

    // Set nodal coordinates
    SetNodalCoordinates_();

}  // Mesh::ConnectMesh_

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

    // Sanity check: static_cast has defined behavior
    assert(size >= 0);

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
    while (n_index < static_cast<std::size_t>(size)) {
        // Find block with the largest number of elements
        std::size_t b_index = 0;
        std::size_t max_num_elements = 0;

        // Loop blocks
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
    for (int i = 0; i < size; ++i) {
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
// Set elements numbering
// ----------------------------------------------------------------------------
void Mesh::SetElementsNumbering_(
    const std::vector<std::size_t> &partitions_start,
    const std::vector<std::size_t> &partitions_size) {
    // Get rank
    const int rank = pwr::MPIUtilities::Rank();

    // Set number of elements per direction in the partitions
    nx_partition_ = partitions_size[3 * rank + 0];
    ny_partition_ = partitions_size[3 * rank + 1];
    nz_partition_ = partitions_size[3 * rank + 2];

    // Sanity check: partition size per direction is less than or equal
    // to global size per direction
    assert(nx_partition_ <= nx_);
    assert(ny_partition_ <= ny_);
    assert(nz_partition_ <= nz_);

    // Set number of elements in the partition
    num_elem_partition_ = nx_partition_ * ny_partition_ * nz_partition_;

    // Sanity check: no zero-element partitions allowed
    assert(num_elem_partition_ > 0);

    // Set starting element index per direction in the partition
    index_x_partition_0_ = partitions_start[3 * rank + 0];
    index_y_partition_0_ = partitions_start[3 * rank + 1];
    index_z_partition_0_ = partitions_start[3 * rank + 2];

    // Sanity check: partition starting index is less than global
    // size per direction
    assert(index_x_partition_0_ < nx_);
    assert(index_y_partition_0_ < ny_);
    assert(index_z_partition_0_ < nz_);

    // Set final element index per direction in the partition
    index_x_partition_f_ = index_x_partition_0_ + nx_partition_ - 1;
    index_y_partition_f_ = index_y_partition_0_ + ny_partition_ - 1;
    index_z_partition_f_ = index_z_partition_0_ + nz_partition_ - 1;

    // Sanity check: partition final index is less than global
    // size per direction
    assert(index_x_partition_f_ < nx_);
    assert(index_y_partition_f_ < ny_);
    assert(index_z_partition_f_ < nz_);

    // Grab neighborhood width
    const std::size_t w = neighborhood_width_;

    // Set number of ghost elements at the start of this block
    nx_ghost_0_ = std::min(w, index_x_partition_0_);
    ny_ghost_0_ = std::min(w, index_y_partition_0_);
    nz_ghost_0_ = std::min(w, index_z_partition_0_);

    // Set number of ghost element at the end of this block
    nx_ghost_f_ = std::min(w, nx_ - index_x_partition_f_ - 1);
    ny_ghost_f_ = std::min(w, ny_ - index_y_partition_f_ - 1);
    nz_ghost_f_ = std::min(w, nz_ - index_z_partition_f_ - 1);

    // Set total number of elements per direction (partition + ghost)
    nx_total_ = nx_partition_ + nx_ghost_0_ + nx_ghost_f_;
    ny_total_ = ny_partition_ + ny_ghost_0_ + ny_ghost_f_;
    nz_total_ = nz_partition_ + nz_ghost_0_ + nz_ghost_f_;

    // Sanity check: total partition size per direction is less than or equal
    // to global size per direction
    assert(nx_total_ <= nx_);
    assert(ny_total_ <= ny_);
    assert(nz_total_ <= nz_);

    // Set total number of elements (partition + ghost)
    num_elem_total_ = nx_total_ * ny_total_ * nz_total_;

    // Set number of ghost elements
    num_elem_ghost_ = num_elem_total_ - num_elem_partition_;

}  // Mesh::SetElementsNumbering_

// ----------------------------------------------------------------------------
// Set elements global id
// ----------------------------------------------------------------------------
void Mesh::SetElementsGlobalId_() {
    // Starting index per direction (including ghost elements)
    // Underflow safe: nx_ghost_0_ = std::min(w, index_x_partition_0_)
    // Underflow safe: ny_ghost_0_ = std::min(w, index_y_partition_0_)
    // Underflow safe: nz_ghost_0_ = std::min(w, index_z_partition_0_)
    const std::size_t index_x_0 = index_x_partition_0_ - nx_ghost_0_;
    const std::size_t index_y_0 = index_y_partition_0_ - ny_ghost_0_;
    const std::size_t index_z_0 = index_z_partition_0_ - nz_ghost_0_;

    // Sanity check: partition + ghost starting index is less than global
    // size per direction
    assert(index_x_0 < nx_);
    assert(index_y_0 < ny_);
    assert(index_z_0 < nz_);

    // Ending index per direction (including ghost elements)
    const std::size_t index_x_f = index_x_partition_f_ + nx_ghost_f_;
    const std::size_t index_y_f = index_y_partition_f_ + ny_ghost_f_;
    const std::size_t index_z_f = index_z_partition_f_ + nz_ghost_f_;

    // Sanity check: partition + ghost final index is less than global
    // size per direction
    assert(index_x_f < nx_);
    assert(index_y_f < ny_);
    assert(index_z_f < nz_);

    // NOTE: Global element ids should be zero-index and based on looping
    //       in the x-, y-, and z-directions (in order)

    // Resize
    elem_index_global_.resize(3 * num_elem_total_,
                              std::numeric_limits<std::size_t>::max());

    // Resize
    elem_id_global_.resize(num_elem_total_,
                           std::numeric_limits<std::size_t>::max());

    // Indices for partition and ghost elements
    std::size_t e_partition = 0;
    std::size_t e_ghost = num_elem_partition_;

    // Loop elements
    for (std::size_t z_i = index_z_0; z_i <= index_z_f; ++z_i) {
        for (std::size_t y_i = index_y_0; y_i <= index_y_f; ++y_i) {
            for (std::size_t x_i = index_x_0; x_i <= index_x_f; ++x_i) {
                // Set global element id
                const std::size_t gid = x_i + (nx_ * y_i) + (nx_ * ny_ * z_i);

                // Current element is in the partition
                if (index_z_partition_0_ <= z_i &&
                    z_i <= index_z_partition_f_ &&
                    index_y_partition_0_ <= y_i &&
                    y_i <= index_y_partition_f_ &&
                    index_x_partition_0_ <= x_i &&
                    x_i <= index_x_partition_f_) {
                    // Sanity check: index is less than number of partition
                    // elements
                    assert(e_partition < num_elem_partition_);

                    // Save global indices
                    elem_index_global_[3 * e_partition + 0] = x_i;
                    elem_index_global_[3 * e_partition + 1] = y_i;
                    elem_index_global_[3 * e_partition + 2] = z_i;

                    // Save global element id
                    elem_id_global_[e_partition] = gid;

                    // Increment partition index
                    e_partition++;

                    // Skip the rest
                    continue;
                }

                // Sanity check: index is more than number of partition elements
                assert(num_elem_partition_ <= e_ghost);

                // Sanity check: index is less than number of (partition +
                // ghost) elements
                assert(e_ghost < num_elem_total_);

                // Save global indices
                elem_index_global_[3 * e_ghost + 0] = x_i;
                elem_index_global_[3 * e_ghost + 1] = y_i;
                elem_index_global_[3 * e_ghost + 2] = z_i;

                // Save global element id
                elem_id_global_[e_ghost] = gid;

                // Increment ghost index
                e_ghost++;
            }
        }
    }

    // Sanity check: looped across all partition elements
    assert(e_partition == num_elem_partition_);

    // Sanity check: looped across all ghost elements
    assert(e_ghost == num_elem_total_);

    // Loop elements
    for (std::size_t e = 0; e < num_elem_total_; ++e) {
        // Grab each element index
        const std::size_t x_i = elem_index_global_[3 * e + 0];
        const std::size_t y_i = elem_index_global_[3 * e + 1];
        const std::size_t z_i = elem_index_global_[3 * e + 2];

        // Sanity check: global element indices have been updated from
        // the initialized value
        assert(x_i < std::numeric_limits<std::size_t>::max());
        assert(y_i < std::numeric_limits<std::size_t>::max());
        assert(z_i < std::numeric_limits<std::size_t>::max());

        // Sanity check: global element indices are less than global
        // number of elements per direction
        assert(x_i < nx_);
        assert(y_i < ny_);
        assert(z_i < nz_);

        // Grab each element id
        const std::size_t gid = elem_id_global_[e];

        // Sanity check: global element ids have been updated from
        // the initialized value
        assert(gid < std::numeric_limits<std::size_t>::max());

        // Sanity check: global element ids are less than global
        // number of elements
        assert(gid < num_elem_);
    }

    // Resize global to local vector
    elem_id_local_.resize(num_elem_, std::numeric_limits<std::size_t>::max());

    // Fill global to local vector
    for (std::size_t e = 0; e < elem_id_global_.size(); ++e) {
        elem_id_local_[elem_id_global_[e]] = e;
    }

}  // Mesh::SetElementsGlobalId_

// ----------------------------------------------------------------------------
// Set elements neighborhood
// ----------------------------------------------------------------------------
void Mesh::SetElementsNeighborhood_() {
    // Sanity check: static_cast has defined behavior
    const std::size_t max_int =
        static_cast<std::size_t>(std::numeric_limits<int>::max());
    assert(neighborhood_width_ < max_int);

    // Set number of neighbors per element (depending on basis function order)
    const int w = static_cast<int>(neighborhood_width_);
    const std::size_t num_neighbors =
        (2 * w + 1) * (2 * w + 1) * (2 * w + 1) - 1;

    // Resize based on number neighbors per element
    elem_neighborhood_.resize(num_neighbors * num_elem_partition_);

    // Default neighbor id is local element id
    for (std::size_t e = 0; e < num_elem_partition_; ++e) {
        for (std::size_t n = 0; n < num_neighbors; ++n) {
            elem_neighborhood_[num_neighbors * e + n] = e;
        }
    }

    // Loop elements
    for (std::size_t e = 0; e < num_elem_partition_; ++e) {
        // Grab global element id
        const std::size_t gid = elem_id_global_[e];

        // // Decompose into global indices in each direction
        // const std::size_t x_i = gid % nx_;
        // const std::size_t y_i = (gid / nx_) % ny_;
        // const std::size_t z_i = gid / (nx_ * ny_);

        // Grab global indices in each direction
        const std::size_t x_i = elem_index_global_[3 * e + 0];
        const std::size_t y_i = elem_index_global_[3 * e + 1];
        const std::size_t z_i = elem_index_global_[3 * e + 2];

        // Sanity check: static_cast has defined behavior
        assert(x_i < max_int);
        assert(y_i < max_int);
        assert(z_i < max_int);
        assert(nx_ < max_int);
        assert(ny_ < max_int);
        assert(nz_ < max_int);

        // Count
        std::size_t n = 0;

        // Loop directions
        for (int dz_i = -w; dz_i <= w; ++dz_i) {
            // Set neighbor global index z-direction
            const int nz_i = static_cast<int>(z_i) + dz_i;

            for (int dy_i = -w; dy_i <= w; ++dy_i) {
                // Set neighbor global index y-direction
                const int ny_i = static_cast<int>(y_i) + dy_i;

                for (int dx_i = -w; dx_i <= w; ++dx_i) {
                    // Set neighbor global index x-direction
                    const int nx_i = static_cast<int>(x_i) + dx_i;

                    // Check for outside of global bounds
                    if (nx_i < 0 || static_cast<int>(nx_) <= nx_i || ny_i < 0 ||
                        static_cast<int>(ny_) <= ny_i || nz_i < 0 ||
                        static_cast<int>(nz_) <= nz_i) {
                        // Increment and skip the rest
                        n++;
                        continue;
                    }

                    // Compute global element id
                    const std::size_t n_gid =
                        nx_i + (nx_ * ny_i) + (nx_ * ny_ * nz_i);

                    // Check for center
                    if (gid == n_gid) {
                        // Do *not* increment and skip the rest
                        continue;
                    }

                    // Sanity check: global neighbor id is within bounds
                    assert(n_gid < num_elem_);

                    // Neighbor's local element id
                    const std::size_t n_lid = elem_id_local_[n_gid];

                    // Sanity check: local id has been initialized
                    assert(n_lid < std::numeric_limits<std::size_t>::max());

                    // Save neighbor local element id
                    elem_neighborhood_[num_neighbors * e + n] = n_lid;

                    // Increment
                    n++;
                }
            }
        }
    }

}  // Mesh::SetElementsNeighborhood_

// ----------------------------------------------------------------------------
// Set elements connectivity
// ----------------------------------------------------------------------------
void Mesh::SetElementsConnectivity_() {
    // Connectivity ordering matches gmsh convetion for hex1
    /*
    //  7----------6
    //  |\         |\
    //  | \        | \
    //  |  \       |  \
    //  |   3------|---2
    //  |   |      |   |      y
    //  4----------5   |      ^
    //   \  |       \  |   z  |
    //    \ |        \ |    \ |
    //     \|         \|     \|
    //      0----------1      +------> x
    */

    // Resize based on 8 nodes per element
    conn_.resize(8 * num_elem_total_, std::numeric_limits<std::size_t>::max());

    // Set nodes in x-y plane
    const std::size_t slab_xy = nodes_x_ * nodes_y_;

    // Loop elements
    for (std::size_t e = 0; e < num_elem_total_; ++e) {
        // Grab element indices
        const std::size_t x_i = elem_index_global_[3 * e + 0];
        const std::size_t y_i = elem_index_global_[3 * e + 1];
        const std::size_t z_i = elem_index_global_[3 * e + 2];

        // Helpful plus one index for y and z
        const std::size_t x_ip1 = x_i + 1;
        const std::size_t y_ip1 = y_i + 1;
        const std::size_t z_ip1 = z_i + 1;

        // Compute node ids
        const std::size_t n0 = x_i + (nodes_x_ * y_i) + (slab_xy * z_i);
        const std::size_t n1 = n0 + 1;

        const std::size_t n2 = x_ip1 + (nodes_x_ * y_ip1) + (slab_xy * z_i);
        const std::size_t n3 = n2 - 1;

        const std::size_t n4 = x_i + (nodes_x_ * y_i) + (slab_xy * z_ip1);
        const std::size_t n5 = n4 + 1;

        const std::size_t n6 = x_ip1 + (nodes_x_ * y_ip1) + (slab_xy * z_ip1);
        const std::size_t n7 = n6 - 1;

        // Save node ids
        conn_[8 * e + 0] = n0;
        conn_[8 * e + 1] = n1;
        conn_[8 * e + 2] = n2;
        conn_[8 * e + 3] = n3;
        conn_[8 * e + 4] = n4;
        conn_[8 * e + 5] = n5;
        conn_[8 * e + 6] = n6;
        conn_[8 * e + 7] = n7;
    }

    // Loop elements
    for (std::size_t e = 0; e < num_elem_total_; ++e) {
        // Loop nodes
        for (std::size_t n = 0; n < 8; ++n) {
            // Grab each node id
            const std::size_t nid = conn_[8 * e + n];

            // Sanity check: global node ids have been updated from
            // the initialized value
            assert(nid < std::numeric_limits<std::size_t>::max());

            // Sanity check: global node ids are less than global
            // number of nodes
            assert(nid < num_nodes_);
        }
    }

}  // Mesh::SetElementsConnectivity_

// ----------------------------------------------------------------------------
// Set nodal coordinates
// ----------------------------------------------------------------------------
void Mesh::SetNodalCoordinates_() {
    // Resize based on {x,y,z} per node
    nodal_coords_.resize(3 * num_nodes_, std::numeric_limits<double>::max());

    // Keep track of which nodes have been assigned
    std::vector<char> node_assigned(num_nodes_, 0);

    // Create offset map
    const std::vector<double> offset{
        0.,  0.,  0.,   // node 0
        dx_, 0.,  0.,   // node 1
        dx_, dy_, 0.,   // node 2
        0.,  dy_, 0.,   // node 3
        0.,  0.,  dz_,  // node 4
        dx_, 0.,  dz_,  // node 5
        dx_, dy_, dz_,  // node 6
        0.,  dy_, dz_,  // node 7
    };

    // Loop elements
    for (std::size_t e = 0; e < num_elem_total_; ++e) {
        // Grab element indices
        const std::size_t x_i = elem_index_global_[3 * e + 0];
        const std::size_t y_i = elem_index_global_[3 * e + 1];
        const std::size_t z_i = elem_index_global_[3 * e + 2];

        // Compute coordinates at element corner
        const double e_x = static_cast<double>(x_i) * dx_;
        const double e_y = static_cast<double>(y_i) * dy_;
        const double e_z = static_cast<double>(z_i) * dz_;

        // Loop nodes
        for (std::size_t n = 0; n < 8; ++n) {
            // Grab current node
            const std::size_t nid = conn_[8 * e + n];

            // Sanity check: node id is reasonable
            assert(nid < num_nodes_);

            // Go to next node if this node has alreday been assigned
            if (node_assigned[nid] != 0) {
                continue;
            }

            // Save nodal coordinates
            nodal_coords_[3 * nid + 0] = e_x + offset[3 * n + 0];
            nodal_coords_[3 * nid + 1] = e_y + offset[3 * n + 1];
            nodal_coords_[3 * nid + 2] = e_z + offset[3 * n + 2];

            // Update assigned nodes
            node_assigned[nid] = 1;
        }
    }

    // Loop nodes
    for (std::size_t n = 0; n < num_nodes_; ++n) {
        // Go to next node if this node is not assigned
        if (node_assigned[n] == 0) {
            continue;
        }

        // Sanity check: nodal coordinates in each direction are assigned
        assert(nodal_coords_[3 * n + 0] < std::numeric_limits<double>::max());
        assert(nodal_coords_[3 * n + 1] < std::numeric_limits<double>::max());
        assert(nodal_coords_[3 * n + 2] < std::numeric_limits<double>::max());
    }

}  // Mesh::SetNodalCoordinates_

}  // namespace pwr
