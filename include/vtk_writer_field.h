#ifndef VTK_WRITER_FIELD_H
#define VTK_WRITER_FIELD_H

#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <cassert>
#include <memory>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>

#include "field.h"
#include "mesh.h"
#include "utilities.h"
#include "vtk_writer_base.h"

namespace pwr {

// Forward declare MeshBase
// This is sufficient since the header only stores a pointer-like object
// and not the complete type
class MeshBase;

// vtk writer class for field
class VTKWriterField : public VTKWriterBase {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    VTKWriterField(const std::shared_ptr<const pwr::MeshBase>& mesh)
        : mesh_(mesh) {
        // Sanity check: pointer to mesh is not null
        assert(mesh_);

        Setup_();
    }

    // ------------------------------------------------------------------------
    // Destructor
    // ------------------------------------------------------------------------
    ~VTKWriterField() override = default;

    // ------------------------------------------------------------------------
    // Delete copy constructor
    // ------------------------------------------------------------------------
    VTKWriterField(const VTKWriterField&) = delete;

    // ------------------------------------------------------------------------
    // Delete copy assignment
    // ------------------------------------------------------------------------
    VTKWriterField& operator=(const VTKWriterField&) = delete;

    // ------------------------------------------------------------------------
    // Finalize setup
    // ------------------------------------------------------------------------
    void FinalizeSetup() override {
        // Nothing to do here... `Setup_` can be called by constructor
    }

    // ------------------------------------------------------------------------
    // Add field
    // ------------------------------------------------------------------------
    void AddField(const std::shared_ptr<pwr::Field>& field,
                  const std::string& field_name) override {
        assert(setup_complete_);
        fields_.emplace_back(field);
        fields_names_.emplace_back(field_name);
    }

    // ------------------------------------------------------------------------
    // Add particles
    // ------------------------------------------------------------------------
    void AddParticles(
        const std::vector<std::shared_ptr<pwr::ParticleBase>>& /*particles*/,
        const std::string& /*particles_name*/) override {
        assert(setup_complete_);
        std::stringstream error_message;
        error_message << "`AddParticles` not supported in "
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

    // ------------------------------------------------------------------------
    // Write time output file
    // ------------------------------------------------------------------------
    void WriteTime(
        const std::string& filename,
        const std::vector<std::string>& piece_filenames) const override {
        assert(setup_complete_);
        WriteTime_(filename, piece_filenames);
    }

   private:
    // ------------------------------------------------------------------------
    // Setup vtk writer for field
    // Calls:
    //   InitializeWriter_()
    //   CheckSetup_()
    // ------------------------------------------------------------------------
    void Setup_();

    // ------------------------------------------------------------------------
    // Initialize vtk writer for field
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

    // ------------------------------------------------------------------------
    // Write time output file
    // ------------------------------------------------------------------------
    void WriteTime_(const std::string& filename,
                    const std::vector<std::string>& piece_filenames) const;

    // CheckSetup_() has been successfully called
    bool setup_complete_ = false;

    // Shared pointer to mesh object
    const std::shared_ptr<const pwr::MeshBase> mesh_;

    // Container to hold shared pointer to field objects
    std::vector<std::shared_ptr<pwr::Field>> fields_;

    // Container to hold name for each field
    std::vector<std::string> fields_names_;

    // Smart vtk pointer to points
    vtkSmartPointer<vtkPoints> points_;

    // Smart vtk pointer to grid
    vtkSmartPointer<vtkUnstructuredGrid> grid_;
};

}  // namespace pwr

#endif  // VTK_WRITER_FIELD_H
