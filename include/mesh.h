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

    // MeshCheck_ has been successfully called
    bool setup_complete_ = false;

    // Minimum x-dir
    double x_min_;

    // Minimum y-dir
    double y_min_;

    // Minimum z-dir
    double z_min_;

    // Maximum x-dir
    double x_max_;

    // Maximum y-dir
    double y_max_;

    // Maximum z-dir
    double z_max_;

    // Number of elements x-dir
    std::size_t nx_;

    // Number of elements y-dir
    std::size_t ny_;

    // Number of elements z-dir
    std::size_t nz_;

    // Total number of elements in the mesh
    std::size_t num_elem_;

    // Number of elements owned by this partition
    std::size_t num_elem_partition_;

    // Number of ghost elements for this partition
    std::size_t num_elem_ghost_;

    // Element size x-dir
    double dx_;

    // Element size y-dir
    double dy_;

    // Element size z-dir
    double dz_;

    // Element-wise connectivity with size 8*TODO with
    // convention TODO
    std::vector<std::size_t> conn_;

    // Nodal coordinates with size 3*TODO with
    // convetion [x0, y0, z0, x1, y1, z1, ..., xN, yN, zN]
    std::vector<double> nodal_coords_;
};

}  // namespace pwr
#endif  // PWR_MESH_H
