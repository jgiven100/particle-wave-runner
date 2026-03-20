#ifndef PWR_OUTPUT_MANAGER_H
#define PWR_OUTPUT_MANAGER_H

#include <vtk_writer_base.h>

#include <cassert>
#include <cstddef>
#include <memory>
#include <string>

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
                  std::size_t max_step)
        : mesh_(mesh), max_step_(max_step) {
        // Sanity check: pointer to mesh is not null
        assert(mesh_);

        Setup_();
    }

    // ------------------------------------------------------------------------
    // Constructor delegation
    // ------------------------------------------------------------------------
    OutputManager(const std::shared_ptr<const pwr::MeshBase>& mesh)
        : OutputManager(mesh, 0) {}

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
    // Write mesh files
    // ------------------------------------------------------------------------
    void WriteMeshFiles() {
        assert(setup_complete_);
        WriteMeshFiles_();
    }

    // ------------------------------------------------------------------------
    // Write particle files
    // ------------------------------------------------------------------------
    void WriteParticleFiles() {
        assert(setup_complete_);
        WriteParticleFiles_();
    }

   private:
    // ------------------------------------------------------------------------
    // Setup output manager
    // Calls:
    //   InitializeMeshWriter_()
    //   InitializeParticleWriter_()
    //   CheckSetup_()
    // ------------------------------------------------------------------------
    void Setup_();

    // ------------------------------------------------------------------------
    // Initialize output manager
    // ------------------------------------------------------------------------
    void InitializeOutputManager_();

    // ------------------------------------------------------------------------
    // Initialize mesh writer
    // ------------------------------------------------------------------------
    void InitializeMeshWriter_();

    // ------------------------------------------------------------------------
    // Initialize particle writer
    // ------------------------------------------------------------------------
    void InitializeParticleWriter_();

    // ------------------------------------------------------------------------
    // Check output manager setup
    // Making it this far means that `setup_complete_` should be set to `true`
    // ------------------------------------------------------------------------
    void CheckSetup_();

    // ------------------------------------------------------------------------
    // Write mesh files
    // ------------------------------------------------------------------------
    void WriteMeshFiles_();

    // ------------------------------------------------------------------------
    // Write particle files
    // -----------------------------------------------------------------------
    void WriteParticleFiles_();

    // CheckSetup_() has been successfully called
    bool setup_complete_ = false;

    // Shared pointer to mesh object
    const std::shared_ptr<const pwr::MeshBase> mesh_;

    // Shared pointer to vtk writer object for mesh
    std::shared_ptr<const pwr::VTKWriterBase> vtk_writer_mesh_;

    // Shared pointer to vtk writer object for particles
    std::shared_ptr<const pwr::VTKWriterBase> vtk_writer_particles_;

    // Output directory name
    std::string output_dir_name_;

    // Maximum step
    std::size_t max_step_;

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
