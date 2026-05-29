#ifndef PWR_OUTPUT_MANAGER_H
#define PWR_OUTPUT_MANAGER_H

#include <cassert>
#include <cstddef>
#include <memory>
#include <string>

#include "field.h"
#include "mesh_base.h"
#include "particle_base.h"
#include "vtk_writer_base.h"

namespace pwr {

// Forward declare MeshBase
// This is sufficient since the header only stores a pointer-like object
// and not the complete type
class MeshBase;

// Output manager class
class OutputManager {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    OutputManager(const std::shared_ptr<const pwr::MeshBase>& mesh,
                  std::size_t max_step, std::size_t output_step)
        : mesh_(mesh), max_step_(max_step), output_step_(output_step) {
        // Sanity check: pointer to mesh is not null
        assert(mesh_);

        Setup_();
    }

    // ------------------------------------------------------------------------
    // Destructor
    // ------------------------------------------------------------------------
    ~OutputManager() = default;

    // ------------------------------------------------------------------------
    // Delete copy constructor
    // ------------------------------------------------------------------------
    OutputManager(const OutputManager&) = delete;

    // ------------------------------------------------------------------------
    // Delete copy assignment
    // ------------------------------------------------------------------------
    OutputManager& operator=(const OutputManager&) = delete;

    // ------------------------------------------------------------------------
    // Add field
    // ------------------------------------------------------------------------
    void AddField(const std::shared_ptr<pwr::Field>& field,
                  const std::string& field_name) {
        assert(setup_complete_);
        vtk_writer_field_->AddField(field, field_name);
    }

    // ------------------------------------------------------------------------
    // Add particles
    // ------------------------------------------------------------------------
    void AddParticles(
        const std::vector<std::shared_ptr<pwr::ParticleBase>>& particles,
        const std::string& particles_name) {
        // assert(setup_complete_);  // Setup_ called after adding particles
        vtk_writer_particles_->AddParticles(particles, particles_name);
    }

    // ------------------------------------------------------------------------
    // Finalize writers setup
    // ------------------------------------------------------------------------
    void FinalizeWritersSetup() {
        assert(setup_complete_);
        FinalizeWritersSetup_();
    }

    // ------------------------------------------------------------------------
    // Write field files
    // ------------------------------------------------------------------------
    void WriteFieldFiles(const std::size_t step) {
        assert(setup_complete_);
        WriteFieldFiles_(step);
    }

    // ------------------------------------------------------------------------
    // Write field files time
    // ------------------------------------------------------------------------
    void WriteFieldFilesTime() {
        assert(setup_complete_);
        WriteFieldFilesTime_();
    }

    // ------------------------------------------------------------------------
    // Write mesh files
    // ------------------------------------------------------------------------
    void WriteMeshFiles(const std::size_t step) {
        assert(setup_complete_);
        WriteMeshFiles_(step);
    }

    // ------------------------------------------------------------------------
    // Write mesh files time
    // ------------------------------------------------------------------------
    void WriteMeshFilesTime() {
        assert(setup_complete_);
        WriteMeshFilesTime_();
    }

    // ------------------------------------------------------------------------
    // Write particles files
    // ------------------------------------------------------------------------
    void WriteParticlesFiles(const std::size_t step) {
        assert(setup_complete_);
        WriteParticlesFiles_(step);
    }

    // ------------------------------------------------------------------------
    // Write particles files time
    // ------------------------------------------------------------------------
    void WriteParticlesFilesTime() {
        assert(setup_complete_);
        WriteParticlesFilesTime_();
    }

   private:
    // ------------------------------------------------------------------------
    // Setup output manager
    // Calls:
    //   InitializeOutputManager_();
    //   InitializeFieldWriter_();
    //   InitializeMeshWriter_()
    //   InitializeParticleWriter_()
    //   SetTimeSteps_()
    //   CheckSetup_()
    // ------------------------------------------------------------------------
    void Setup_();

    // ------------------------------------------------------------------------
    // Initialize output manager
    // ------------------------------------------------------------------------
    void InitializeOutputManager_();

    // ------------------------------------------------------------------------
    // Initialize field writer
    // ------------------------------------------------------------------------
    void InitializeFieldWriter_();

    // ------------------------------------------------------------------------
    // Initialize mesh writer
    // ------------------------------------------------------------------------
    void InitializeMeshWriter_();

    // ------------------------------------------------------------------------
    // Initialize particle writer
    // ------------------------------------------------------------------------
    void InitializeParticleWriter_();

    // ------------------------------------------------------------------------
    // Set time and steps
    // ------------------------------------------------------------------------
    void SetTimeSteps_();

    // ------------------------------------------------------------------------
    // Check output manager setup
    // Making it this far means that `setup_complete_` should be set to `true`
    // ------------------------------------------------------------------------
    void CheckSetup_();

    // ------------------------------------------------------------------------
    // Finalize writers setup
    // ------------------------------------------------------------------------
    void FinalizeWritersSetup_();

    // ------------------------------------------------------------------------
    // Write field files
    // ------------------------------------------------------------------------
    void WriteFieldFiles_(const std::size_t step);

    // ------------------------------------------------------------------------
    // Write field files time
    // ------------------------------------------------------------------------
    void WriteFieldFilesTime_();

    // ------------------------------------------------------------------------
    // Write mesh files
    // ------------------------------------------------------------------------
    void WriteMeshFiles_(const std::size_t step);

    // ------------------------------------------------------------------------
    // Write mesh files time
    // ------------------------------------------------------------------------
    void WriteMeshFilesTime_();

    // ------------------------------------------------------------------------
    // Write particles files
    // ------------------------------------------------------------------------
    void WriteParticlesFiles_(const std::size_t step);

    // ------------------------------------------------------------------------
    // Write particles files time
    // ------------------------------------------------------------------------
    void WriteParticlesFilesTime_();

    // CheckSetup_() has been successfully called
    bool setup_complete_ = false;

    // Shared pointer to mesh object
    const std::shared_ptr<const pwr::MeshBase> mesh_;

    // Shared pointer to vtk writer object for field
    std::shared_ptr<pwr::VTKWriterBase> vtk_writer_field_;

    // Shared pointer to vtk writer object for mesh
    std::shared_ptr<const pwr::VTKWriterBase> vtk_writer_mesh_;

    // Shared pointer to vtk writer object for particles
    std::shared_ptr<pwr::VTKWriterBase> vtk_writer_particles_;

    // Output directory name
    std::string output_dir_name_;

    // Maximum step
    std::size_t max_step_;

    // Output step
    std::size_t output_step_;

    // Number of files
    std::size_t steps_;

    // Step padding size
    int step_padding_;

    // Rank padding size
    int rank_padding_;

    // Rank
    int rank_;

    // Size
    int size_;
};

}  // namespace pwr

#endif  // PWR_OUTPUT_MANAGER_H
