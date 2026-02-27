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
    void SetElementsGlobalId_(const std::vector<std::size_t>& partitions_start,
                              const std::vector<std::size_t>& partitions_size);

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

    // Minimum global x-dir
    double x_min_;

    // Minimum global y-dir
    double y_min_;

    // Minimum global z-dir
    double z_min_;

    // Maximum global x-dir
    double x_max_;

    // Maximum global y-dir
    double y_max_;

    // Maximum global z-dir
    double z_max_;

    // Number of global elements x-dir
    std::size_t nx_;

    // Number of global elements y-dir
    std::size_t ny_;

    // Number of global elements z-dir
    std::size_t nz_;

    // Total number of elements in the mesh (all partitions)
    std::size_t num_elem_;

    // Number of elements owned by this partition
    std::size_t num_elem_partition_;

    // Number of ghost elements for this partition
    std::size_t num_elem_ghost_;

    // Number of partition + ghost elements
    std::size_t num_elem_total_;

    // Element size x-dir
    double dx_;

    // Element size y-dir
    double dy_;

    // Element size z-dir
    double dz_;

    // Neighborhood width
    std::size_t neighborhood_width_ = 1;  // p=1 B-Spline (linear)
    // std::size_t neighborhood_width_ = 2; // p=2 B-Spline (quadratic)
    // std::size_t neighborhood_width_ = 2; // p=3 B-Spline (cubic)
    // std::size_t neighborhood_width_ = 3; // p=4 B-Spline (quartic)

    // Element global id with size (num_elem_partition_ + num_elem_ghost_)
    std::vector<std::size_t> elem_id_global_;

    // Element neighborhood with size (TODO)
    std::vector<std::size_t> elem_neighborhood_;

    // Element connectivity with size (TODO)
    std::vector<std::size_t> conn_;

    // Nodal coordinates with size 3 * (num_elem_partition_ + num_elem_ghost_)
    std::vector<double> nodal_coords_;
};

}  // namespace pwr
#endif  // PWR_MESH_H
