#include "mesh.h"

#include <algorithm>
#include <array>
#include <catch2/catch_all.hpp>
#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <vector>

#include "mesh_base.h"
#include "mpi_utilities.h"

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
    const std::array<std::size_t, 16> num_elem_partition = {
        27, 0,  0,  0,   // MPI size = 1
        9,  18, 0,  0,   // MPI size = 2
        6,  9,  12, 0,   // MPI size = 3
        4,  6,  9,  8};  // MPI size = 4

    // Set solution (number of known ghost elements)
    // Computed by hand [02 April 2026]
    const std::array<std::size_t, 16> num_elem_ghost = {
        0,  0,  0,  0,    // MPI size = 1
        9,  9,  0,  0,    // MPI size = 2
        12, 9,  15, 0,    // MPI size = 3
        14, 12, 9,  19};  // MPI size = 4

    // Node indices for each face
    const std::array<std::size_t, 24> face_indices = {
        0, 1, 2, 3,   // Bottom face
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
        for (const auto &mesh : meshes) {
            // Set index for flattened array
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
        for (const auto &mesh : meshes) {
            // Set index for flattened array
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
        for (const auto &mesh : meshes) {
            // Set index for flattened array
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
        const std::array<std::size_t, 16> num_nodes_active = {
            64, 0,  0,  0,    // MPI size = 1
            48, 64, 0,  0,    // MPI size = 2
            48, 48, 64, 0,    // MPI size = 3
            48, 48, 48, 64};  // MPI size = 4

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto &mesh : meshes) {
            // Set index for flattened array
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
    // Local element connectivity
    // ------------------------------------------------------------------------
    {
        INFO("Local element connectivity");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Loop each 3x3x3 mesh
        for (const auto &mesh : meshes) {
            // Grab connectivity (using local node ids)
            const auto &conn = mesh->GetElemConnLocal();

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
            std::map<std::array<std::size_t, 4>, int> faces;

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
                for (const auto &n : elem) {
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
                    std::array<std::size_t, 4> face = {n0, n1, n2, n3};
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
            for (const auto &[_, value] : faces) {
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
        for (const auto &mesh : meshes) {
            // Grab connectivity (using global node ids)
            const auto &conn = mesh->GetElemConnGlobal();

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
            std::map<std::array<std::size_t, 4>, int> faces;

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
                for (const auto &n : elem) {
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
                    std::array<std::size_t, 4> face = {n0, n1, n2, n3};
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
            for (const auto &[_, value] : faces) {
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
        for (const auto &mesh : meshes) {
            // Grab nodal ownership
            const auto &owned = mesh->GetNodalOwnership();

            // Count owned on this partition
            const std::size_t owned_local =
                std::count(owned.begin(), owned.end(), 1);

            // Reduce number of owned
            const std::size_t owned_global =
                pwr::MPIUtilities::AllReduceSum(owned_local);

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
        for (const auto &mesh : meshes) {
            // Grab nodal coordinates
            const auto &coords = mesh->GetNodalCoordinates();

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
        for (const auto &mesh : meshes) {
            // Grab connectivity (using local node ids)
            const auto &conn = mesh->GetElemConnLocal();

            // Grab nodal coordinates
            const auto &coords = mesh->GetNodalCoordinates();

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
                if (std::fabs(dvol - elem_vol) > 1.e-12) {
                    num_failed_tests_local++;
                }
            }

            // Grab number of active nodes
            const std::size_t num_nodes_active = mesh->GetNumNodesActive();

            // Create set to store unique nodes
            std::set<std::array<double, 3>> nodes;

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

}  // TEST_CASE("Mesh")
