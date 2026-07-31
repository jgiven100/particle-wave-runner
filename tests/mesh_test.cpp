#include "mesh.h"

#include <algorithm>
#include <catch2/catch_all.hpp>
#include <cmath>
#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <vector>

#include "mesh_base.h"
#include "mpi_utilities.h"
#include "numbers.h"

using pwr::numbers::TOL;

TEST_CASE("Mesh", "[mesh]") {
    // Grab rank and size
    const int rank = pwr::MPIUtilities::Rank();
    const int size = pwr::MPIUtilities::Size();

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
    const std::size_t nx = 3;
    const std::size_t ny = 3;
    const std::size_t nz = 3;

    // Create shared pointer to mesh object with full constructor
    meshes[0] = std::make_shared<const pwr::Mesh>(x_min, y_min, z_min, x_max,
                                                  y_max, z_max, nx, ny, nz);

    // Create shared pointer to mesh object with constructor delegation
    meshes[1] =
        std::make_shared<const pwr::Mesh>(x_max, y_max, z_max, nx, ny, nz);

    // Create shared pointer to mesh object with constructor delegation
    meshes[2] = std::make_shared<const pwr::Mesh>(nx, ny, nz);

    // Set solution (number of partition elements)
    // Computed by hand [02 April 2026]
    const std::vector<std::size_t> num_elem_partition = {
        27, 0,  0,  0,   // MPI size = 1
        9,  18, 0,  0,   // MPI size = 2
        6,  9,  12, 0,   // MPI size = 3
        4,  6,  9,  8};  // MPI size = 4

    // Set solution (number of known ghost elements)
    // Computed by hand [02 April 2026]
    const std::vector<std::size_t> num_elem_ghost = {
        0,  0,  0,  0,    // MPI size = 1
        9,  9,  0,  0,    // MPI size = 2
        12, 9,  15, 0,    // MPI size = 3
        14, 12, 9,  19};  // MPI size = 4

    // Node indices for each face
    const std::vector<std::size_t> face_indices = {0, 1, 2, 3,   // Bottom face
                                                   4, 5, 6, 7,   // Top face
                                                   0, 1, 5, 4,   // South face
                                                   3, 0, 4, 7,   // West face
                                                   3, 2, 6, 7,   // North face
                                                   2, 1, 5, 6};  // East face

    // ------------------------------------------------------------------------
    // Number of partition elements
    // ------------------------------------------------------------------------
    {
        INFO("Number of partition elements");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto& mesh : meshes) {
            // Set index for flattened vector
            const int index = 4 * (size - 1) + rank;

            // -- Check if getter matches hand calc'd solution ----------------
            if (mesh->GetNumElemPartition() != num_elem_partition[index]) {
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

    }  // Number of partition elements

    // ------------------------------------------------------------------------
    // Number of ghost elements
    // ------------------------------------------------------------------------
    {
        INFO("Number of ghost elements");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto& mesh : meshes) {
            // Set index for flattened vector
            const int index = 4 * (size - 1) + rank;

            // -- Check if getter matches hand calc'd solution ----------------
            if (mesh->GetNumElemGhost() != num_elem_ghost[index]) {
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

    }  // Number of ghost elements

    // ------------------------------------------------------------------------
    // Number of partition + ghost elements
    // ------------------------------------------------------------------------
    {
        INFO("Number of partition + ghost elements");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto& mesh : meshes) {
            // Set index for flattened vector
            const int index = 4 * (size - 1) + rank;

            // Number of partition + ghost elements
            const std::size_t num_elem_total =
                num_elem_partition[index] + num_elem_ghost[index];

            // -- Check if getter matches hand calc'd solution ----------------
            if (mesh->GetNumElemTotal() != num_elem_total) {
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

    }  // Number of partition + ghost elements

    // ------------------------------------------------------------------------
    // Number of active nodes
    // ------------------------------------------------------------------------
    {
        INFO("Number of active nodes");

        // Set solution (number of active nodes)
        // Computed by hand [02 April 2026]
        const std::vector<std::size_t> num_nodes_active = {
            64, 0,  0,  0,    // MPI size = 1
            48, 64, 0,  0,    // MPI size = 2
            48, 48, 64, 0,    // MPI size = 3
            48, 48, 48, 64};  // MPI size = 4

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto& mesh : meshes) {
            // Set index for flattened vector
            const int index = 4 * (size - 1) + rank;

            // -- Check if getter matches hand calc'd solution ----------------
            if (mesh->GetNumNodesActive() != num_nodes_active[index]) {
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

    }  // Number of active nodes

    // ------------------------------------------------------------------------
    // Global element index
    // ------------------------------------------------------------------------
    {
        INFO("Global element index");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto& mesh : meshes) {
            // Grab global element indices
            const auto& elem_index_global = mesh->GetElemIndexGlobal();

            // Grab global element id
            const auto& elem_id_global = mesh->GetElemIdGlobal();

            // Loop partition + ghost elements
            for (std::size_t e = 0; e < mesh->GetNumElemTotal(); ++e) {
                // Set index in each direction
                const std::size_t x_i = elem_index_global[3 * e + 0];
                const std::size_t y_i = elem_index_global[3 * e + 1];
                const std::size_t z_i = elem_index_global[3 * e + 2];

                // -- Check if index is greater than global max ---------------
                if (x_i >= nx) {
                    num_failed_tests_local++;
                }
                if (y_i >= ny) {
                    num_failed_tests_local++;
                }
                if (z_i >= nz) {
                    num_failed_tests_local++;
                }

                // Compute expected global id
                const std::size_t gid = x_i + nx * y_i + nx * ny * z_i;

                // -- Check if calc'd id matches the saved valued -------------
                if (gid != elem_id_global[e]) {
                    num_failed_tests_local++;
                }
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

    }  // Global element index

    // ------------------------------------------------------------------------
    // Global element id
    // ------------------------------------------------------------------------
    {
        INFO("Global element id");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto& mesh : meshes) {
            // Grab global element id
            const auto& elem_id_global = mesh->GetElemIdGlobal();

            // Loop each id
            for (const auto id : elem_id_global) {
                // -- Check if id is greater than global max ------------------
                if (id >= nx * ny * nz) {
                    num_failed_tests_local++;
                }
            }

            // Grab number of partition elements
            const std::size_t num_elem_partition = mesh->GetNumElemPartition();

            // Vector with just owned element ids
            const std::vector<std::size_t> elem_id_owned(
                elem_id_global.begin(),
                elem_id_global.begin() + num_elem_partition);

            // Collect all rank's owned element ids
            std::vector<std::size_t> elem_id_all;
            pwr::MPIUtilities::AllGatherV(elem_id_owned, elem_id_all);

            // -- Check if combined element ids is correct size ---------------
            if (elem_id_all.size() != (nx * ny * nz)) {
                num_failed_tests_local++;
            }

            // Create std::set (no repeats) for combined owned element ids
            const std::set elem_id_all_set(elem_id_all.begin(),
                                           elem_id_all.end());

            // -- Check for duplicates in combined element ids ----------------
            if (elem_id_all.size() != elem_id_all_set.size()) {
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

    }  // Global element id

    // ------------------------------------------------------------------------
    // Local element connectivity
    // ------------------------------------------------------------------------
    {
        INFO("Local element connectivity");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto& mesh : meshes) {
            // Grab connectivity (using local node ids)
            const auto& conn = mesh->GetElemConnLocal();

            // Grab number of partition + ghost elements
            const std::size_t num_elem_total = mesh->GetNumElemTotal();

            // Grab number of active nodes
            const std::size_t num_nodes_active = mesh->GetNumNodesActive();

            // -- Check if each element has 8 nodes ---------------------------
            if (conn.size() / 8 != num_elem_total) {
                num_failed_tests_local++;
            }

            // Create set of sets to store unique nodes for each element
            std::set<std::set<std::size_t>> elems;

            // Create map to count occurances for each face
            std::map<std::vector<std::size_t>, int> faces;

            // Loop elements
            for (std::size_t e = 0; e < num_elem_total; ++e) {
                // Create set with nodes for each element
                std::set<std::size_t> elem(conn.begin() + (e * 8),
                                           conn.begin() + ((e + 1) * 8));

                // -- Check if each element has 8 unique nodes ----------------
                if (elem.size() != 8) {
                    num_failed_tests_local++;
                }

                // Loop nodes
                for (const auto& n : elem) {
                    // -- Check if each node id is reasonable -----------------
                    if (n >= num_nodes_active) {
                        num_failed_tests_local++;
                    }
                }

                // -- Check if there are duplicate elements -------------------
                const auto it = elems.find(elem);
                if (it != elems.end()) {
                    num_failed_tests_local++;
                }
                elems.insert(elem);

                // Loop faces and count shared faces
                for (int f = 0; f < 6; ++f) {
                    // Set nodes and sort
                    const std::size_t n0 =
                        conn[8 * e + face_indices[4 * f + 0]];
                    const std::size_t n1 =
                        conn[8 * e + face_indices[4 * f + 1]];
                    const std::size_t n2 =
                        conn[8 * e + face_indices[4 * f + 2]];
                    const std::size_t n3 =
                        conn[8 * e + face_indices[4 * f + 3]];
                    std::vector<std::size_t> face = {n0, n1, n2, n3};
                    std::sort(face.begin(), face.end());

                    // Find current face
                    const auto it = faces.find(face);
                    if (it != faces.end()) {
                        // If face already exists, increment and move to
                        // next face
                        it->second++;
                        continue;
                    }
                    faces.insert({face, 1});
                }
            }

            // -- Check if each face is found at most 2 times -----------------
            for (const auto& [_, value] : faces) {
                if (value > 2) {
                    num_failed_tests_local++;
                }
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

    }  // Local element connectivity

    // ------------------------------------------------------------------------
    // Global element connectivity
    // ------------------------------------------------------------------------
    {
        INFO("Global element connectivity");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto& mesh : meshes) {
            // Grab connectivity (using global node ids)
            const auto& conn = mesh->GetElemConnGlobal();

            // Grab number of partition + ghost elements
            const std::size_t num_elem_total = mesh->GetNumElemTotal();

            // Set number of global nodes
            const std::size_t num_nodes_global = (nx + 1) * (ny + 1) * (nz + 1);

            // -- Check if each element has 8 nodes ---------------------------
            if (conn.size() / 8 != num_elem_total) {
                num_failed_tests_local++;
            }

            // Create set of sets to store unique nodes for each element
            std::set<std::set<std::size_t>> elems;

            // Create map to count occurances for each face
            std::map<std::vector<std::size_t>, int> faces;

            // Loop elements
            for (std::size_t e = 0; e < num_elem_total; ++e) {
                // Create set with nodes for each element
                std::set<std::size_t> elem(conn.begin() + (e * 8),
                                           conn.begin() + ((e + 1) * 8));

                // -- Check if each element has 8 unique nodes ----------------
                if (elem.size() != 8) {
                    num_failed_tests_local++;
                }

                // Loop nodes
                for (const auto& n : elem) {
                    // -- Check if each node id is reasonable -----------------
                    if (n >= num_nodes_global) {
                        num_failed_tests_local++;
                    }
                }

                // -- Check if there are duplicate elements -------------------
                const auto it = elems.find(elem);
                if (it != elems.end()) {
                    num_failed_tests_local++;
                }
                elems.insert(elem);

                // Loop faces and count shared faces
                for (int f = 0; f < 6; ++f) {
                    // Set nodes and sort
                    const std::size_t n0 =
                        conn[8 * e + face_indices[4 * f + 0]];
                    const std::size_t n1 =
                        conn[8 * e + face_indices[4 * f + 1]];
                    const std::size_t n2 =
                        conn[8 * e + face_indices[4 * f + 2]];
                    const std::size_t n3 =
                        conn[8 * e + face_indices[4 * f + 3]];
                    std::vector<std::size_t> face = {n0, n1, n2, n3};
                    std::sort(face.begin(), face.end());

                    // Find current face
                    const auto it = faces.find(face);
                    if (it != faces.end()) {
                        // If face already exists, increment and move to
                        // next face
                        it->second++;
                        continue;
                    }
                    faces.insert({face, 1});
                }
            }

            // -- Check if each face is found at most 2 times -----------------
            for (const auto& [_, value] : faces) {
                if (value > 2) {
                    num_failed_tests_local++;
                }
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

    }  // Global element connectivity

    // ------------------------------------------------------------------------
    // Nodal ownership
    // ------------------------------------------------------------------------
    {
        INFO("Nodal ownership");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Set total nodes (all partitions)
        const std::size_t total_nodes = (nx + 1) * (ny + 1) * (nz + 1);

        // Loop each 3x3x3 mesh
        for (const auto& mesh : meshes) {
            // Grab nodal ownership
            const auto& owned = mesh->GetNodalOwnership();

            // Count owned on this partition
            const std::size_t owned_local =
                std::count(owned.begin(), owned.end(), 1);

            // Reduce number of owned
            std::size_t owned_global;
            pwr::MPIUtilities::AllReduceSum(owned_local, owned_global);

            // -- Check if sum of owned matches number active nodes -----------
            if (owned_global != total_nodes) {
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

    }  // Nodal ownership

    // ------------------------------------------------------------------------
    // Nodal coordinates
    // ------------------------------------------------------------------------
    {
        INFO("Nodal coordinates");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto& mesh : meshes) {
            // Grab nodal coordinates
            const auto& coords = mesh->GetNodalCoordinates();

            // Grab number of active nodes
            const std::size_t num_nodes_active = mesh->GetNumNodesActive();

            // -- Check if each node has 3 components {x,y,z} -----------------
            if (coords.size() / 3 != num_nodes_active) {
                num_failed_tests_local++;
            }

            // Loop nodes
            for (std::size_t n = 0; n < num_nodes_active; ++n) {
                // -- Check if x-coordinate is within bounds ------------------
                const double x = coords[n * 3 + 0];
                if ((x < x_min) || (x_max < x)) {
                    num_failed_tests_local++;
                }

                // -- Check if y-coordinate is within bounds ------------------
                const double y = coords[n * 3 + 1];
                if ((y < y_min) || (y_max < y)) {
                    num_failed_tests_local++;
                }

                // -- Check if z-coordinate is within bounds ------------------
                const double z = coords[n * 3 + 2];
                if ((z < z_min) || (z_max < z)) {
                    num_failed_tests_local++;
                }
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

    }  // Nodal coordinates

    // ------------------------------------------------------------------------
    // Element volume
    // ------------------------------------------------------------------------
    {
        INFO("Element volume");

        // Compute element volume
        const double elem_vol = (x_max - x_min) * (y_max - y_min) *
                                (z_max - z_min) / (nx * ny * nz);

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto& mesh : meshes) {
            // Grab connectivity (using local node ids)
            const auto& conn = mesh->GetElemConnLocal();

            // Grab nodal coordinates
            const auto& coords = mesh->GetNodalCoordinates();

            // Grab number of partition elements
            const std::size_t num_elem_partition = mesh->GetNumElemPartition();

            // Loop elements
            for (std::size_t e = 0; e < num_elem_partition; ++e) {
                // Grab near corner
                const std::size_t n0 = conn[8 * e + 0];
                const double x0 = coords[3 * n0 + 0];
                const double y0 = coords[3 * n0 + 1];
                const double z0 = coords[3 * n0 + 2];

                // Grab far corner
                const std::size_t n7 = conn[8 * e + 7];
                const double x7 = coords[3 * n7 + 0];
                const double y7 = coords[3 * n7 + 1];
                const double z7 = coords[3 * n7 + 2];

                // Compute volume
                const double dvol = (x7 - x0) * (y7 - y0) * (z7 - z0);

                // -- Check if volume is close --------------------------------
                if (std::fabs(dvol - elem_vol) > TOL) {
                    num_failed_tests_local++;
                }
            }

            // Grab number of active nodes
            const std::size_t num_nodes_active = mesh->GetNumNodesActive();

            // Create set to store unique nodes
            std::set<std::vector<double>> nodes;

            // Loop nodes
            for (std::size_t n = 0; n < num_nodes_active; ++n) {
                nodes.insert(
                    {coords[3 * n + 0], coords[3 * n + 1], coords[3 * n + 2]});
            }

            // -- Check if each active node is found --------------------------
            if (nodes.size() != num_nodes_active) {
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

    }  // Element volume

    // ------------------------------------------------------------------------
    // Find containing element global id
    // ------------------------------------------------------------------------
    {
        INFO("Find containing element global id");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto& mesh : meshes) {
            // Grab connectivity (using local node ids)
            const auto& conn = mesh->GetElemConnLocal();

            // Grab nodal coordinates
            const auto& coords = mesh->GetNodalCoordinates();

            // Grab number of partition + ghost elements
            const std::size_t num_elem_total = mesh->GetNumElemTotal();

            // Grab global element id
            const auto& gid = mesh->GetElemIdGlobal();

            // Loop elements
            for (std::size_t e = 0; e < num_elem_total; ++e) {
                // Grab near corner
                const std::size_t n0 = conn[8 * e + 0];
                const double x0 = coords[3 * n0 + 0];
                const double y0 = coords[3 * n0 + 1];
                const double z0 = coords[3 * n0 + 2];

                // Grab far corner
                const std::size_t n7 = conn[8 * e + 7];
                const double x7 = coords[3 * n7 + 0];
                const double y7 = coords[3 * n7 + 1];
                const double z7 = coords[3 * n7 + 2];

                // Set center
                const std::vector<double> center = {
                    (x0 + x7) * 0.5, (y0 + y7) * 0.5, (z0 + z7) * 0.5};

                // Find containing element global id
                const std::size_t gid_itr =
                    mesh->FindContainingElemIdGlobal(center);

                // -- Check if global element ids matches ----------------------
                if (gid[e] != gid_itr) {
                    num_failed_tests_local++;
                }
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

    }  // Find containing element global id

    // ------------------------------------------------------------------------
    // Find containing element global id -- serial
    // ------------------------------------------------------------------------
    {
        INFO("Find containing element global id -- serial");

        if (size == 1) {
            // Keep track of failed tests
            int num_failed_tests = 0;

            // Set solution (global element id)
            // Computed by hand [24 July 2026]
            const std::map<std::size_t, std::vector<double>> test_cases = {
                {0, {1. / 6., 1. / 6., 1. / 6.}},    // Element indices: {0,0,0}
                {1, {3. / 6., 1. / 6., 1. / 6.}},    // Element indices: {1,0,0}
                {4, {3. / 6., 3. / 6., 1. / 6.}},    // Element indices: {1,1,0}
                {3, {1. / 6., 3. / 6., 1. / 6.}},    // Element indices: {0,1,0}
                {9, {1. / 6., 1. / 6., 3. / 6.}},    // Element indices: {0,0,1}
                {10, {3. / 6., 1. / 6., 3. / 6.}},   // Element indices: {1,0,1}
                {13, {3. / 6., 3. / 6., 3. / 6.}},   // Element indices: {1,1,1}
                {12, {1. / 6., 3. / 6., 3. / 6.}},   // Element indices: {0,1,1}
                {26, {5. / 6., 5. / 6., 5. / 6.}}};  // Element indices: {2,2,2}

            // Loop each 3x3x3 mesh
            for (const auto& mesh : meshes) {
                // Loop test cases
                for (const auto& [gid, coords] : test_cases) {
                    // Find containing element global id
                    const std::size_t gid_itr =
                        mesh->FindContainingElemIdGlobal(coords);

                    // -- Check if global element ids matches -----------------
                    if (gid != gid_itr) {
                        num_failed_tests++;
                    }
                }
            }

            // Only check on root (rank 0)
            REQUIRE(num_failed_tests == 0);
        }

    }  // Find containing element global id -- serial

    // ------------------------------------------------------------------------
    // Compute local coordinates
    // ------------------------------------------------------------------------
    {
        INFO("Compute local coordinates");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto& mesh : meshes) {
            // Grab connectivity (using local node ids)
            const auto& conn = mesh->GetElemConnLocal();

            // Grab nodal coordinates
            const auto& coords = mesh->GetNodalCoordinates();

            // Grab number of partition + ghost elements
            const std::size_t num_elem_total = mesh->GetNumElemTotal();

            // Grab global element id
            const auto& gid = mesh->GetElemIdGlobal();

            // Loop elements
            for (std::size_t e = 0; e < num_elem_total; ++e) {
                // Grab near corner
                const std::size_t n0 = conn[8 * e + 0];
                const double x0 = coords[3 * n0 + 0];
                const double y0 = coords[3 * n0 + 1];
                const double z0 = coords[3 * n0 + 2];

                // Grab far corner
                const std::size_t n7 = conn[8 * e + 7];
                const double x7 = coords[3 * n7 + 0];
                const double y7 = coords[3 * n7 + 1];
                const double z7 = coords[3 * n7 + 2];

                // Set global coordinates
                const std::vector<double> coords_global = {
                    (x0 + x7) * 0.5 - (x7 - x0) * 0.4,
                    (y0 + y7) * 0.5 + (y7 - y0) * 0.1,
                    (z0 + z7) * 0.5 + (z7 - z0) * 0.5};

                // Compute local coordinates
                std::vector<double> coords_local;
                mesh->ComputeLocalCoordinates(gid[e], coords_global,
                                              coords_local);

                // Compute error
                const double error = std::fabs(coords_local[0] + 0.8) +
                                     std::fabs(coords_local[1] - 0.2) +
                                     std::fabs(coords_local[2] - 1.0);

                // -- Check if local coordinates matches known values ---------
                if (error > TOL) {
                    num_failed_tests_local++;
                }
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

    }  // Compute local coordinates

}  // TEST_CASE("Mesh")
