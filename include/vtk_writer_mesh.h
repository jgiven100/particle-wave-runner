#ifndef VTK_WRITER_MESH_H
#define VTK_WRITER_MESH_H

#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <cassert>
#include <memory>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>

#include "utilities.h"
#include "vtk_writer_base.h"

namespace pwr {

// Forward declare MeshBase
// This is sufficient since the header only stores a pointer-like object
// and not the complete type
class MeshBase;

// vtk writer class for mesh
class VTKWriterMesh : public VTKWriterBase {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    VTKWriterMesh(const std::shared_ptr<const pwr::MeshBase>& mesh)
        : mesh_(mesh) {
        // Sanity check: pointer to mesh is not null
        assert(mesh_);

        Setup_();
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
    // Add field
    // ------------------------------------------------------------------------
    void AddField(const std::shared_ptr<pwr::Field>& /*field*/,
                  const std::string& /*field_name*/) override {
        assert(setup_complete_);
        std::stringstream error_message;
        error_message << "`AddField` not supported in "
                      << std::source_location::current().function_name();
        pwr::Utilities::PrintErrorOnRoot(error_message.str());
    }

    // ------------------------------------------------------------------------
    // Write output file
    // ------------------------------------------------------------------------
    void Write(const std::string& filename) const override {
        assert(setup_complete_);
        Write_(filename);
    }

    // ------------------------------------------------------------------------
    // Write parallel output file
    // ------------------------------------------------------------------------
    void WriteParallel(
        const std::string& filename,
        const std::vector<std::string>& piece_filenames) const override {
        assert(setup_complete_);
        WriteParallel_(filename, piece_filenames);
    }

   private:
    // ------------------------------------------------------------------------
    // Setup vtk writer for mesh
    // Calls:
    //   InitializeWriter_()
    //   CheckSetup_()
    // ------------------------------------------------------------------------
    void Setup_();

    // ------------------------------------------------------------------------
    // Initialize vtk writer for mesh
    // ------------------------------------------------------------------------
    void InitializeWriter_();

    // ------------------------------------------------------------------------
    // Check setup
    // Making it this far means that `setup_complete_` should be set to `true`
    // ------------------------------------------------------------------------
    void CheckSetup_();

    // ------------------------------------------------------------------------
    // Write output file
    // ------------------------------------------------------------------------
    void Write_(const std::string& filename) const;

    // ------------------------------------------------------------------------
    // Write parallel output file
    // ------------------------------------------------------------------------
    void WriteParallel_(const std::string& filename,
                        const std::vector<std::string>& piece_filenames) const;

    // CheckSetup_() has been successfully called
    bool setup_complete_ = false;

    // Shared pointer to mesh object
    const std::shared_ptr<const pwr::MeshBase> mesh_;

    // Smart vtk pointer to points
    vtkSmartPointer<vtkPoints> points_;

    // Smart vtk pointer to grid
    vtkSmartPointer<vtkUnstructuredGrid> grid_;
};

}  // namespace pwr

#endif  // VTK_WRITER_MESH_H
