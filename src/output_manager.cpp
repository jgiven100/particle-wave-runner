#include "output_manager.h"

#include <algorithm>
#include <filesystem>
#include <memory>
#include <sstream>
#include <string>

#include "mpi_utilities.h"
#include "utilities.h"
#include "vtk_writer_field.h"
#include "vtk_writer_mesh.h"
#include "vtk_writer_particles.h"

namespace pwr {

// ----------------------------------------------------------------------------
// Setup output manager
// ----------------------------------------------------------------------------
void OutputManager::Setup_() {
    InitializeOutputManager_();
    InitializeFieldWriter_();
    InitializeMeshWriter_();
    InitializeParticleWriter_();
    SetTimeSteps_();
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
        std::stringstream info_message;
        info_message << "Output files are save at: " << output_dir;
        Utilities::PrintInfoOnRoot(info_message.str());
    }

    // Ensure `output/` directory exists on all ranks before proceeding
    MPIUtilities::Barrier();

    // Set step padding
    step_padding_ = std::max(1, Utilities::CountDigits(max_step_));

    // Set rank padding
    rank_padding_ = std::max(1, Utilities::CountDigits(size_));

}  // OutputManager::InitializeOutputManager_

// ----------------------------------------------------------------------------
// Initialize VTK writer for field
// ----------------------------------------------------------------------------
void OutputManager::InitializeFieldWriter_() {
    // Shared pointer to vtk writer object for field
    vtk_writer_field_ = std::make_shared<VTKWriterField>(mesh_);

    // Sanity check: pointer to vtk writer is not null
    assert(vtk_writer_field_);

}  // OutputManager::InitializeFieldWriter_

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
    // Shared pointer to vtk writer object for particles
    vtk_writer_particles_ = std::make_shared<VTKWriterParticles>(mesh_);

    // Sanity check: pointer to vtk writer is not null
    assert(vtk_writer_particles_);

}  // OutputManager::InitializeParticleWriter_

// ----------------------------------------------------------------------------
// Set time and steps
// ----------------------------------------------------------------------------
void OutputManager::SetTimeSteps_() {
    // Sanity check: max_step_ and output_step_ are consistent
    assert(max_step_ % output_step_ == 0);

    // Set number of output steps (files) generated; `+1` to count 0th step
    steps_ = (max_step_ / output_step_) + 1;

}  // OutputManager::SetTimeSteps_()

// ----------------------------------------------------------------------------
// Check output manager
// ----------------------------------------------------------------------------
void OutputManager::CheckSetup_() {
    // Setup is complete
    setup_complete_ = true;

}  // OutputManager::CheckSetup_

// ----------------------------------------------------------------------------
// Setup writers
// ----------------------------------------------------------------------------
void OutputManager::FinalizeWritersSetup_() {
    // Must call after adding particles
    vtk_writer_particles_->FinalizeSetup();

}  // OutputManager::FinalizeWritersSetup_

// ----------------------------------------------------------------------------
// Write output files for field
// ----------------------------------------------------------------------------
void OutputManager::WriteFieldFiles_(const std::size_t step) {
    // String padding for step
    std::string step_str = std::to_string(step);
    Utilities::PadString(step_str, step_padding_);

    // String padding for rank
    std::string rank_str = std::to_string(rank_);
    Utilities::PadString(rank_str, rank_padding_);

    // Output filename
    const std::string filename_vtu =
        output_dir_name_ + "field_step" + step_str + "_" + rank_str + ".vtu";

    // Write file
    vtk_writer_field_->Write(filename_vtu);

    // Ensure each `.vtu` is written before `.pvtu` is written
    MPIUtilities::Barrier();

    // Only write pvtu on root
    if (rank_ == 0) {
        // Output filename
        const std::string filename_pvtu =
            output_dir_name_ + "field_step" + step_str + ".pvtu";

        // Place to hold each rank's filename
        std::vector<std::string> piece_filenames(size_);

        // Loop each rank and save file name
        for (int r = 0; r < size_; ++r) {
            // String padding for rank
            std::string r_str = std::to_string(r);
            Utilities::PadString(r_str, rank_padding_);

            // Save name
            piece_filenames[r] = "field_step" + step_str + "_" + r_str + ".vtu";
        }

        // Write parallel file
        vtk_writer_field_->WriteParallel(filename_pvtu, piece_filenames);
    }

}  // OutputManager::WriteFieldFiles_

// ----------------------------------------------------------------------------
// Write output time files for field
// ----------------------------------------------------------------------------
void OutputManager::WriteFieldFilesTime_() {
    // Only write pvd on root
    if (rank_ == 0) {
        // Output filename
        const std::string filename_pvd = output_dir_name_ + "field.pvd";

        // Place to hold each timestep's filename
        std::vector<std::string> piece_filenames(steps_);

        // Loop each output step and save file name
        for (std::size_t s = 0; s < steps_; ++s) {
            // String padding for output step
            std::string s_str = std::to_string(s * output_step_);
            Utilities::PadString(s_str, step_padding_);

            // Save name
            piece_filenames[s] = "field_step" + s_str + ".pvtu";
        }

        // Write time file
        vtk_writer_field_->WriteTime(filename_pvd, piece_filenames);
    }

}  // OutputManager::WriteFieldFilesTime_

// ----------------------------------------------------------------------------
// Write output files for mesh
// ----------------------------------------------------------------------------
void OutputManager::WriteMeshFiles_(const std::size_t step) {
    // String padding for step
    std::string step_str = std::to_string(step);
    Utilities::PadString(step_str, step_padding_);

    // String padding for rank
    std::string rank_str = std::to_string(rank_);
    Utilities::PadString(rank_str, rank_padding_);

    // Output filename
    const std::string filename_vtu =
        output_dir_name_ + "mesh_step" + step_str + "_" + rank_str + ".vtu";

    // Write file
    vtk_writer_mesh_->Write(filename_vtu);

    // Ensure each `.vtu` is written before `.pvtu` is written
    MPIUtilities::Barrier();

    // Only write pvtu on root
    if (rank_ == 0) {
        // Output filename
        const std::string filename_pvtu =
            output_dir_name_ + "mesh_step" + step_str + ".pvtu";

        // Place to hold each rank's filename
        std::vector<std::string> piece_filenames(size_);

        // Loop each rank and save file name
        for (int r = 0; r < size_; ++r) {
            // String padding for rank
            std::string r_str = std::to_string(r);
            Utilities::PadString(r_str, rank_padding_);

            // Save name
            piece_filenames[r] = "mesh_step" + step_str + "_" + r_str + ".vtu";
        }

        // Write parallel file
        vtk_writer_mesh_->WriteParallel(filename_pvtu, piece_filenames);
    }

}  // OutputManager::WriteMeshFiles_

// ----------------------------------------------------------------------------
// Write output time files for mesh
// ----------------------------------------------------------------------------
void OutputManager::WriteMeshFilesTime_() {
    // Only write pvd on root
    if (rank_ == 0) {
        // Output filename
        const std::string filename_pvd = output_dir_name_ + "mesh.pvd";

        // Place to hold each timestep's filename
        std::vector<std::string> piece_filenames(steps_);

        // Loop each output step and save file name
        for (std::size_t s = 0; s < steps_; ++s) {
            // String padding for output step
            std::string s_str = std::to_string(s * output_step_);
            Utilities::PadString(s_str, step_padding_);

            // Save name
            piece_filenames[s] = "mesh_step" + s_str + ".pvtu";
        }

        // Write time file
        vtk_writer_mesh_->WriteTime(filename_pvd, piece_filenames);
    }

}  // OutputManager::WriteMeshFilesTime_

// ----------------------------------------------------------------------------
// Write output files for particles
// ----------------------------------------------------------------------------
void OutputManager::WriteParticlesFiles_(const std::size_t step) {
    // String padding for step
    std::string step_str = std::to_string(step);
    Utilities::PadString(step_str, step_padding_);

    // String padding for rank
    std::string rank_str = std::to_string(rank_);
    Utilities::PadString(rank_str, rank_padding_);

    // Output filename
    const std::string filename_vtu = output_dir_name_ + "particles_step" +
                                     step_str + "_" + rank_str + ".vtu";

    // Write file
    vtk_writer_particles_->Write(filename_vtu);

    // Ensure each `.vtu` is written before `.pvtu` is written
    MPIUtilities::Barrier();

    // Only write pvtu on root
    if (rank_ == 0) {
        // Output filename
        const std::string filename_pvtu =
            output_dir_name_ + "particles_step" + step_str + ".pvtu";

        // Place to hold each rank's filename
        std::vector<std::string> piece_filenames(size_);

        // Loop each rank and save file name
        for (int r = 0; r < size_; ++r) {
            // String padding for rank
            std::string r_str = std::to_string(r);
            Utilities::PadString(r_str, rank_padding_);

            // Save name
            piece_filenames[r] =
                "particles_step" + step_str + "_" + r_str + ".vtu";
        }

        // Write parallel file
        vtk_writer_particles_->WriteParallel(filename_pvtu, piece_filenames);
    }

}  // OutputManager::WriteParticlesFiles_

// ----------------------------------------------------------------------------
// Write output time files for particles
// ----------------------------------------------------------------------------
void OutputManager::WriteParticlesFilesTime_() {
    // Only write pvd on root
    if (rank_ == 0) {
        // Output filename
        const std::string filename_pvd = output_dir_name_ + "particles.pvd";

        // Place to hold each timestep's filename
        std::vector<std::string> piece_filenames(steps_);

        // Loop each output step and save file name
        for (std::size_t s = 0; s < steps_; ++s) {
            // String padding for output step
            std::string s_str = std::to_string(s * output_step_);
            Utilities::PadString(s_str, step_padding_);

            // Save name
            piece_filenames[s] = "particles_step" + s_str + ".pvtu";
        }

        // Write time file
        vtk_writer_particles_->WriteTime(filename_pvd, piece_filenames);
    }

}  // OutputManager::WriteParticlesFilesTime_

};  // namespace pwr
