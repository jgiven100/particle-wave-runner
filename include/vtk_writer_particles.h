#ifndef VTK_WRITER_PARTICLES_H
#define VTK_WRITER_PARTICLES_H

#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <cassert>
#include <memory>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>

#include "particle_base.h"
#include "utilities.h"
#include "vtk_writer_base.h"

namespace pwr {

// Forward declare MeshBase
// This is sufficient since the header only stores a pointer-like object
// and not the complete type
class MeshBase;

// vtk writer class for particles
class VTKWriterParticles : public VTKWriterBase {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    VTKWriterParticles(const std::shared_ptr<const pwr::MeshBase>& mesh)
        : mesh_(mesh) {
        // Sanity check: point to mesh is not null
        assert(mesh_);

        // Need to call `AddParticles()` before calling `Setup_`!
        // Setup_();
    }

    // ------------------------------------------------------------------------
    // Destructor
    // ------------------------------------------------------------------------
    ~VTKWriterParticles() override = default;

    // ------------------------------------------------------------------------
    // Delete copy constructor
    // ------------------------------------------------------------------------
    VTKWriterParticles(const VTKWriterParticles&) = delete;

    // ------------------------------------------------------------------------
    // Delete copy assignment
    // ------------------------------------------------------------------------
    VTKWriterParticles& operator=(const VTKWriterParticles&) = delete;

    // ------------------------------------------------------------------------
    // Finalize setup
    // ------------------------------------------------------------------------
    void FinalizeSetup() override {
        // Call after adding particles
        Setup_();
    }

    // ------------------------------------------------------------------------
    // Add field
    // ------------------------------------------------------------------------
    void AddField(const std::shared_ptr<pwr::Field>& /*field*/,
                  const std::string& /*field_name*/) override {
        // `Setup_` is *not* called by consturctor for particle writer
        // assert(setup_complete_);
        std::stringstream error_message;
        error_message << "`AddField` not supported in "
                      << std::source_location::current().function_name();
        pwr::Utilities::PrintErrorOnRoot(error_message.str());
    }

    // ------------------------------------------------------------------------
    // Add particles
    // ------------------------------------------------------------------------
    void AddParticles(
        const std::vector<std::shared_ptr<pwr::ParticleBase>>& particles,
        const std::string& particles_name) override {
        // `Setup_` is *not* called by consturctor for particle writer
        assert(!setup_complete_);
        particles_ = particles;
        particles_name_ = particles_name;
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
    // Setup vtk writer for particles
    // Calls:
    //   InitializeWriter_()
    //   CheckSetup_()
    // ------------------------------------------------------------------------
    void Setup_();

    // ------------------------------------------------------------------------
    // Initialize vtk writer for particles
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

    // Container to shared pointers to particle objects
    std::vector<std::shared_ptr<pwr::ParticleBase>> particles_;

    // Container to hold name for shared pointers to particle objects
    std::string particles_name_;

    // Smart vtk pointer to points
    vtkSmartPointer<vtkPoints> points_;

    // Smart vtk point to grid
    vtkSmartPointer<vtkUnstructuredGrid> grid_;
};

}  // namespace pwr

#endif  // VTK_WRITER_PARTICLE_H
