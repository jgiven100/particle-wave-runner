#ifndef VTK_WRITER_MESH_H
#define VTK_WRITER_MESH_H

#include <iostream>  // TODO remove
#include <memory>
#include <string>

#include "vtk_writer_base.h"

namespace pwr {

// Forward declare MeshBase
// This is sufficient since the header only stores a pointer-like object
// and not the complete type
class MeshBase;

// VTK writer class for mesh
class VTKWriterMesh : public VTKWriterBase {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    VTKWriterMesh(const std::shared_ptr<const pwr::MeshBase>& mesh)
        : mesh_(mesh) {
        std::cout << "inside VTKWriterMesh() ctr" << std::endl;
    }

    // ------------------------------------------------------------------------
    // Destructor
    // ------------------------------------------------------------------------
    ~VTKWriterMesh() override = default;

    // ------------------------------------------------------------------------
    // Delete copy constructor
    // ------------------------------------------------------------------------
    VTKWriterMesh(const VTKWriterMesh&) = delete;

    // ------------------------------------------------------------------------
    // Delete copy assignment
    // ------------------------------------------------------------------------
    VTKWriterMesh& operator=(const VTKWriterMesh&) = delete;

    // ------------------------------------------------------------------------
    // Write output file
    // ------------------------------------------------------------------------
    void Write(const std::string& filename) const override;

   private:
    // Shared pointer to mesh object
    const std::shared_ptr<const pwr::MeshBase> mesh_;
};

}  // namespace pwr

#endif  // VTK_WRITER_MESH_H
