#ifndef PWR_MESH_H
#define PWR_MESH_H

#include <cassert>
#include <cstddef>
#include <vector>

#include "mesh_base.h"

namespace pwr {

// Mesh class
// Store information about particles, nodes, cells, and neighbors
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
        SetupMesh_();
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
    ~Mesh() = default;

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

    // Get element-wise connectivity
    const std::vector<std::size_t>& GetElemConnectivity() const override {
        assert(setup_complete_);
        return conn_;
    }

    // Get nodal coordinates
    const std::vector<double>& GetNodalCoordinates() const override {
        assert(setup_complete_);
        return nodal_coords_;
    }

   private:
    // ------------------------------------------------------------------------
    // Setup mesh
    // ------------------------------------------------------------------------
    void SetupMesh_();

    // ------------------------------------------------------------------------
    // Initialize mesh
    // ------------------------------------------------------------------------
    void InitializeMesh_();

    // ------------------------------------------------------------------------
    // Partition mesh
    // ------------------------------------------------------------------------
    void PartitionMesh_();

    // ------------------------------------------------------------------------
    // Check mesh
    // ------------------------------------------------------------------------
    void CheckMesh_();

    // ------------------------------------------------------------------------
    // Set partitions per direction
    // ------------------------------------------------------------------------
    void SetPartitionsPerDirection_(std::vector<std::size_t>& partitions_start,
                                    std::vector<std::size_t>& partitions_size);

    // ------------------------------------------------------------------------
    // Set number of elements
    // ------------------------------------------------------------------------
    void SetNumberElements_(const std::vector<std::size_t>& partitions_start,
                            const std::vector<std::size_t>& partitions_size);

    // ------------------------------------------------------------------------
    // Set elements global id
    // ------------------------------------------------------------------------
    void SetElementsGlobalId_();

    // ------------------------------------------------------------------------
    // Set elements neighborhood
    // ------------------------------------------------------------------------
    void SetElementsNeighborhood_();

    // ------------------------------------------------------------------------
    // Set elements connectivity
    // ------------------------------------------------------------------------
    void SetElementsConnectivity_();

    // MeshCheck_ has been successfully called
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

    // Element global id with size (num_elem_partition_ + num_elem_ghost_)
    // Index is local element id
    // Value is global element id
    std::vector<std::size_t> elem_id_global_;

    // Element local id with size num_elem_
    // Index is gloabl element id
    // Value is local element id
    std::vector<std::size_t> elem_id_local_;

    // Element neighborhood with size (TODO)
    std::vector<std::size_t> elem_neighborhood_;

    // Element connectivity with size (TODO)
    std::vector<std::size_t> conn_;

    // Nodal coordinates with size 3 * (num_elem_partition_ + num_elem_ghost_)
    std::vector<double> nodal_coords_;
};

}  // namespace pwr
#endif  // PWR_MESH_H
