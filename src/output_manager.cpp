#include "output_manager.h"

#include <algorithm>
#include <filesystem>
#include <iostream>  // TODO remove after refactor
#include <memory>
#include <string>

#include "mpi_utilities.h"
#include "utilities.h"
#include "vtk_writer_mesh.h"

namespace pwr {

// ----------------------------------------------------------------------------
// Setup output manager
// ----------------------------------------------------------------------------
void OutputManager::Setup_() {
    InitializeOutputManager_();
    InitializeMeshWriter_();
    InitializeParticleWriter_();
    CheckSetup_();

}  // OutputManager::Setup_

// ----------------------------------------------------------------------------
// Initialize output manager
// ----------------------------------------------------------------------------
void OutputManager::InitializeOutputManager_() {
    // Set rank and size
    rank_ = MPIUtilities::Rank();
    size_ = MPIUtilities::Size();

    // Set output directory name
    output_dir_name_ = "output/";

    // Only worry about file system on root
    if (rank_ == 0) {
        // Get current directory
        std::filesystem::path current_dir = std::filesystem::current_path();
        std::filesystem::path output_dir = current_dir / output_dir_name_;

        // Create the output directory (if needed)
        if (!std::filesystem::exists(output_dir)) {
            const bool new_dir = std::filesystem::create_directory(output_dir);

            // Sanity check: output directory is successfully created
            assert(new_dir);
        }

        // Say where the output directory is located
        // TODO move CLI message to utilities.h later
        std::cout << "Output files are saved at: " << output_dir << std::endl;
    }

    // Ensure `output/` directory exists on all ranks before proceeding
    MPIUtilities::Barrier();

    // Set step padding
    step_padding_ = std::max(1, Utilities::CountDigits(max_step_));

    // Set rank padding
    rank_padding_ = std::max(1, Utilities::CountDigits(size_));

}  // OutputManager::InitializeOutputManager_

// ----------------------------------------------------------------------------
// Initialize VTK writer for mesh
// ----------------------------------------------------------------------------
void OutputManager::InitializeMeshWriter_() {
    // Shared pointer to vtk writer object for mesh
    vtk_writer_mesh_ = std::make_shared<VTKWriterMesh>(mesh_);

    // Sanity check: pointer to vtk writer is not null
    assert(vtk_writer_mesh_);

}  // OutputManager::InitializeMeshWriter_

// ----------------------------------------------------------------------------
// Initialize VTK writer for particles
// ----------------------------------------------------------------------------
void OutputManager::InitializeParticleWriter_() {
    ;  // TODO

}  // OutputManager::InitializeParticleWriter_

// ----------------------------------------------------------------------------
// Check output manager
// ----------------------------------------------------------------------------
void OutputManager::CheckSetup_() {
    // Setup is complete
    setup_complete_ = true;

}  // OutputManager::CheckSetup_

// ----------------------------------------------------------------------------
// Write output files for mesh
// ----------------------------------------------------------------------------
void OutputManager::WriteMeshFiles_() {
    // TODO pass `step` as function input later
    int step = 0;

    // String padding for step
    std::string step_str = std::to_string(step);
    step_str.insert(0, step_padding_ - step_str.size(), '0');

    // String padding for rank
    std::string rank_str = std::to_string(rank_);
    rank_str.insert(0, rank_padding_ - rank_str.size(), '0');

    // Output filename
    const std::string filename_vtu =
        output_dir_name_ + "mesh_step" + step_str + "_" + rank_str + ".vtu";

    // Write file
    vtk_writer_mesh_->Write(filename_vtu);

    // Ensure each `.vtu` is written before `.pvtu` is written
    MPIUtilities::Barrier();

    // Only write ptvu on root
    if (rank_ == 0) {
        // Output filename
        const std::string filename_pvtu =
            output_dir_name_ + "mesh_step" + step_str + ".pvtu";

        // Place to hold each rank's filename
        std::vector<std::string> piece_filenames;

        // Loop each rank and save file name
        for (int r = 0; r < size_; ++r) {
            // String padding for rank
            std::string r_str = std::to_string(r);
            r_str.insert(0, rank_padding_ - r_str.size(), '0');

            // Save name
            piece_filenames.push_back("mesh_step" + step_str + "_" + r_str +
                                      ".vtu");
        }

        // Write parallel file
        vtk_writer_mesh_->WriteParallel(filename_pvtu, piece_filenames);
    }

}  // OutputManager::WriteMeshFiles_

// ----------------------------------------------------------------------------
// Write output files for particles
// ----------------------------------------------------------------------------
void OutputManager::WriteParticleFiles_() {
    ;  // TODO

}  // OutputManager::WriteParticleFiles_

};  // namespace pwr
