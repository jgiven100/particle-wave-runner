#ifndef PWR_MESH_H
#define PWR_MESH_H

#include <cassert>
#include <cstddef>
#include <unordered_map>
#include <vector>

#include "mesh_base.h"

namespace pwr {

// Mesh class
class Mesh : public MeshBase {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    Mesh(double x_min, double y_min, double z_min, double x_max, double y_max,
         double z_max, std::size_t nx, std::size_t ny, std::size_t nz)
        : x_min_(x_min),
          y_min_(y_min),
          z_min_(z_min),
          x_max_(x_max),
          y_max_(y_max),
          z_max_(z_max),
          nx_(nx),
          ny_(ny),
          nz_(nz) {
        Setup_();
    }

    // ------------------------------------------------------------------------
    // Constructor delegation
    // ------------------------------------------------------------------------
    Mesh(double x_max, double y_max, double z_max, std::size_t nx,
         std::size_t ny, std::size_t nz)
        : Mesh(0., 0., 0., x_max, y_max, z_max, nx, ny, nz) {}

    // ------------------------------------------------------------------------
    // Constructor delegation
    // ------------------------------------------------------------------------
    Mesh(std::size_t nx, std::size_t ny, std::size_t nz)
        : Mesh(0., 0., 0., 1., 1., 1., nx, ny, nz) {}

    // ------------------------------------------------------------------------
    // Destructor
    // ------------------------------------------------------------------------
    ~Mesh() override = default;

    // ------------------------------------------------------------------------
    // Delete copy constructor
    // ------------------------------------------------------------------------
    Mesh(const Mesh&) = delete;

    // ------------------------------------------------------------------------
    // Delete copy assignment
    // ------------------------------------------------------------------------
    Mesh& operator=(const Mesh&) = delete;

    // ------------------------------------------------------------------------
    // Getters
    // ------------------------------------------------------------------------

    // Get number of elements partition
    std::size_t GetNumElemPartition() const override {
        assert(setup_complete_);
        return num_elem_partition_;
    }

    // Get number of elements ghost
    std::size_t GetNumElemGhost() const override {
        assert(setup_complete_);
        return num_elem_ghost_;
    }

    // Get number of elements partition + ghost
    std::size_t GetNumElemTotal() const override {
        assert(setup_complete_);
        return num_elem_total_;
    }

    // Get number of active nodes
    std::size_t GetNumNodesActive() const override {
        assert(setup_complete_);
        return num_nodes_active_;
    }

    // Get global element-wise connectivity
    const std::vector<std::size_t>& GetElemConnGlobal() const override {
        assert(setup_complete_);
        return conn_global_;
    }

    // Get local element-wise connectivity
    const std::vector<std::size_t>& GetElemConnLocal() const override {
        assert(setup_complete_);
        return conn_local_;
    }

    // Get nodal ownership
    const std::vector<char>& GetNodalOwnership() const override {
        assert(setup_complete_);
        return nodal_ownership_;
    }

    // Get nodal coordinates
    const std::vector<double>& GetNodalCoordinates() const override {
        assert(setup_complete_);
        return nodal_coords_;
    }

   private:
    // ------------------------------------------------------------------------
    // Setup mesh
    // Calls:
    //   InitializeMesh_()
    //   PartitionMesh_()
    //   NumberMesh_()
    //   ConnectMesh_()
    //   CheckSetup_()
    // ------------------------------------------------------------------------
    void Setup_();

    // ------------------------------------------------------------------------
    // Initialize mesh
    // Checks the provided dimensions and number of elements; computes
    // element size
    // ------------------------------------------------------------------------
    void InitializeMesh_();

    // ------------------------------------------------------------------------
    // Partition mesh
    // Calls:
    //   SetPartitionsPerDirection_()
    //   SetPartitionsOrdering_()
    //   SetElementsNumbering_()
    // ------------------------------------------------------------------------
    void PartitionMesh_();

    // ------------------------------------------------------------------------
    // Number mesh
    // Calls:
    //   SetElementsGlobalId_()
    //   SetElementsNeighborhood_()
    // ------------------------------------------------------------------------
    void NumberMesh_();

    // ------------------------------------------------------------------------
    // Connect mesh
    // Calls:
    //   SetElementsConnectivity_()
    //   SetNodalOwnership_()
    //   SetNodalOwnershipRank_()
    //   SetNodalCoordinates_()
    // ------------------------------------------------------------------------
    void ConnectMesh_();

    // ------------------------------------------------------------------------
    // Check mesh setup
    // Making it this far means that `setup_complete_` should be set to `true`
    // ------------------------------------------------------------------------
    void CheckSetup_();

    // ------------------------------------------------------------------------
    // Set partitions per direction
    // Decompose global domain into blocks based on the number of elements
    // per direction
    // ------------------------------------------------------------------------
    void SetPartitionsPerDirection_();

    // ------------------------------------------------------------------------
    // Set partitions ordering
    // Update partition ordering from split history order to spatial order
    // ------------------------------------------------------------------------
    void SetPartitionsOrdering_();

    // ------------------------------------------------------------------------
    // Set elements numbering
    // For partition and ghost elements correspond to this proc's block, set
    // the number of elements, starting indices, ending indices, etc.
    // ------------------------------------------------------------------------
    void SetElementsNumbering_();

    // ------------------------------------------------------------------------
    // Set elements global id
    // For partition and ghost elements, set
    //   `elem_id_global_` (index is local id, value is global id)
    //   `elem_id_local_`  (key is global id, value is local id)
    // ------------------------------------------------------------------------
    void SetElementsGlobalId_();

    // ------------------------------------------------------------------------
    // Set elements neighborhood
    // The element neighborhood size depends on basis function order; determine
    // shell of surrounding elements for partition element
    // ------------------------------------------------------------------------
    void SetElementsNeighborhood_();

    // ------------------------------------------------------------------------
    // Set elements connectivity
    // Connectivity matches gmsh convention for hex1; connectivity set for
    // partition and ghost elements
    // ------------------------------------------------------------------------
    void SetElementsConnectivity_();

    // ------------------------------------------------------------------------
    // Set nodal ownership
    // Determines which active nodes are owned by this process
    // ------------------------------------------------------------------------
    void SetNodalOwnership_();

    // ------------------------------------------------------------------------
    // Set nodal ownership rank
    // Determines which rank owns each active node
    // ------------------------------------------------------------------------
    void SetNodalOwnershipRank_();

    // ------------------------------------------------------------------------
    // Set nodal coordinates
    // Computes and saves the {x,y,z} coordinates for each node
    // ------------------------------------------------------------------------
    void SetNodalCoordinates_();

    // CheckSetup_() has been successfully called
    bool setup_complete_ = false;

    // Minimum global x-direction
    double x_min_;

    // Minimum global y-direction
    double y_min_;

    // Minimum global z-direction
    double z_min_;

    // Maximum global x-direction
    double x_max_;

    // Maximum global y-direction
    double y_max_;

    // Maximum global z-direction
    double z_max_;

    // Number of global elements x-direction
    std::size_t nx_;

    // Number of global elements y-direction
    std::size_t ny_;

    // Number of global elements z-direction
    std::size_t nz_;

    // Total number of elements in the mesh (all partitions)
    std::size_t num_elem_;

    // Starting index of partition elements x-direction
    std::size_t index_x_partition_0_;

    // Starting index of partition elements y-direction
    std::size_t index_y_partition_0_;

    // Starting index of partition elements z-direction
    std::size_t index_z_partition_0_;

    // Ending index of partition elements x-direction
    std::size_t index_x_partition_f_;

    // Ending index of partition elements y-direction
    std::size_t index_y_partition_f_;

    // Ending index of partition elements z-direction
    std::size_t index_z_partition_f_;

    // Number of partition elements x-direction
    std::size_t nx_partition_;

    // Number of partition elements y-direction
    std::size_t ny_partition_;

    // Number of partition elements z-direction
    std::size_t nz_partition_;

    // Number of elements owned by this partition
    std::size_t num_elem_partition_;

    // Number of ghost elements at start x-direction
    std::size_t nx_ghost_0_;

    // Number of ghost elements at start y-direction
    std::size_t ny_ghost_0_;

    // Number of ghost elements at start z-direction
    std::size_t nz_ghost_0_;

    // Number of ghost elements at end x-direction
    std::size_t nx_ghost_f_;

    // Number of ghost elements at end y-direction
    std::size_t ny_ghost_f_;

    // Number of ghost elements at end z-direction
    std::size_t nz_ghost_f_;

    // Number of ghost elements for this partition
    std::size_t num_elem_ghost_;

    // Number of partition + ghost elements x-direction
    std::size_t nx_total_;

    // Number of partition + ghost elements y-direction
    std::size_t ny_total_;

    // Number of partition + ghost elements z-direction
    std::size_t nz_total_;

    // Number of partition + ghost elements
    std::size_t num_elem_total_;

    // Number of global nodes x-direction
    std::size_t nodes_x_;

    // Number of global nodes y-direction
    std::size_t nodes_y_;

    // Number of global nodes z-direction
    std::size_t nodes_z_;

    // Number of global nodes
    std::size_t num_nodes_;

    // Number of active nodes
    std::size_t num_nodes_active_;

    // Element size x-direction
    double dx_;

    // Element size y-direction
    double dy_;

    // Element size z-direction
    double dz_;

    // Neighborhood width
    std::size_t neighborhood_width_ = 1;  // p=1 B-Spline (linear)
    // std::size_t neighborhood_width_ = 2; // p=2 B-Spline (quadratic)
    // std::size_t neighborhood_width_ = 2; // p=3 B-Spline (cubic)
    // std::size_t neighborhood_width_ = 3; // p=4 B-Spline (quartic)

    // Partitions start vector
    // Size: 3 * (total number of MPI ranks)
    // Flat convention:
    //     [start_x0,start_y0,start_z0,...,start_xN,start_yN,start_zN]
    std::vector<std::size_t> partitions_start_;

    // Partitions size vector
    // Size: 3 * (total number of MPI ranks)
    // Flat convention: [size_x0,size_y0,size_z0,...,size_xN,size_yN,size_zN]
    std::vector<std::size_t> partitions_size_;

    // Partitions order vector
    // Size: total number of MPI ranks
    // Index: rank
    // Value: block index
    std::vector<std::size_t> partitions_order_;

    // Element local-to-global id vector
    // Size: number of partition + ghost elements (num_elem_total_)
    // Index: local element id
    // Value: global element id
    std::vector<std::size_t> elem_id_global_;

    // Element global-to-local id map
    // Size: number of partition + ghost elements (num_elem_total_)
    // Key: global element id
    // Value: local element id
    std::unordered_map<std::size_t, std::size_t> elem_id_local_;

    // Element global index
    // Size: (3 * num_elem_total_)
    // Flat convention: [e0x_i,e0y_i,e0z_i,e1x_i,e1y_i,e1z_i,...]
    std::vector<std::size_t> elem_index_global_;

    // Element neighborhood
    // Size: (num_neighbors * num_elem_partition_),
    // where num_neighbors = (2 * neighborhood_width_ + 1) ^ 3 - 1
    // Flat convention: [e0n0,e0n1,...,e0nN,e1n0,e1n1,...,e1nN,...]
    std::vector<std::size_t> elem_neighborhood_;

    // Element connectivity using global node ids
    // Size: (8 * num_elem_total_)
    // Flat convention: [e0n0,e0n1,...,e0n7,e1n0,e1n1,...,e1n7,...]
    std::vector<std::size_t> conn_global_;

    // Element connectivity using local node ids
    // Size: (8 * num_elem_total_)
    // Flat convention: [e0n0,e0n1,...,e0n7,e1n0,e1n1,...,e1n7,...]
    std::vector<std::size_t> conn_local_;

    // Nodal local-to-global id vector
    // Size: number of active nodes (num_nodes_active_)
    // Index: local node id
    // Value: global node id
    std::vector<std::size_t> nodal_id_global_;

    // Nodal global-to-local id map
    // Size: number of active nodes (num_nodes_active_)
    // Key: global node id
    // Value: local node id
    std::unordered_map<std::size_t, std::size_t> nodal_id_local_;

    // Nodal ownership boolean
    // Size: number of active nodes (num_nodes_active_)
    std::vector<char> nodal_ownership_;

    // Nodal ownership global element id
    // Size: number of active nodes (num_nodes_active_)
    std::vector<std::size_t> nodal_ownership_elem_;

    // Nodal ownership rank
    // Size: number of active nodes (num_nodes_active_)
    std::vector<int> nodal_ownership_rank_;

    // Nodal coordinates using local node indexing
    // Size: (3 * num_nodes_active_)
    // Flat convention: [n0x,n0y,n0z,n1x,n1y,n1z,...]
    std::vector<double> nodal_coords_;
};

}  // namespace pwr

#endif  // PWR_MESH_H
