#include "vtk_writer_particles.h"

#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertex.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <cassert>
#include <cstddef>
#include <fstream>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>

#include "mpi_utilities.h"
#include "utilities.h"

namespace pwr {

// ----------------------------------------------------------------------------
// Setup vtk wrtier for particles
// ----------------------------------------------------------------------------
void VTKWriterParticles::Setup_() {
    InitializeWriter_();
    CheckSetup_();

}  // VTKWriterParticles::Setup_

// ----------------------------------------------------------------------------
// Initialize vtk wrtier for particles
// ----------------------------------------------------------------------------
void VTKWriterParticles::InitializeWriter_() {
    // ------------------------------------------------------------------------
    // Assign points
    // ------------------------------------------------------------------------

    // Grab the number of particles owned by this partition
    const std::size_t num_particles = particles_.size();

    // Create vtk points object
    points_ = vtkSmartPointer<vtkPoints>::New();
    points_->SetNumberOfPoints(num_particles);

    // Loop particles  TODO -- loop partition particles
    for (std::size_t p = 0; p < num_particles; ++p) {
        // Grab particle
        const auto& particle = particles_[p];

        // Pointer to start of particle coordinates
        const auto& coords = particle->GetCoordsGlobal();

        // Set coordinates for this point
        points_->SetPoint(p, coords.data());
    }

    // Sanity check: point to points object is not null
    assert(points_);

    // ------------------------------------------------------------------------
    // Create grid and add vertex elements
    // ------------------------------------------------------------------------

    // Create vtk grid object
    grid_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid_->Allocate(num_particles);
    grid_->SetPoints(points_);

    // Loop particles  TODO -- partition particles
    for (std::size_t p = 0; p < num_particles; ++p) {
        // Create vtk vertex object
        auto vertex = vtkSmartPointer<vtkVertex>::New();

        // Set point id
        vertex->GetPointIds()->SetId(0, p);

        // Add vertex element to grid
        grid_->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
    }

    // Sanity check: pointer to grid object is not null
    assert(grid_);

}  // VTKWriterParticles::InitializeWriter_

// ----------------------------------------------------------------------------
// Check setup
// ----------------------------------------------------------------------------
void VTKWriterParticles::CheckSetup_() {
    // Setup is complete
    setup_complete_ = true;

}  // VTKWriterParticles::CheckSetup_

// ----------------------------------------------------------------------------
// Write output file
// ----------------------------------------------------------------------------
void VTKWriterParticles::Write_(const std::string& filename) const {
    // Clear previous arrays
    auto point_data = grid_->GetPointData();
    point_data->Initialize();

    // Set the number of particles  TODO -- partition particles
    const std::size_t num_particles = particles_.size();

    // ------------------------------------------------------------------------
    // Save particles
    // ------------------------------------------------------------------------

    // Set solution name
    const auto& soln_name = particles_name_;

    // Create array to save particles
    auto particles_array = vtkSmartPointer<vtkDoubleArray>::New();
    particles_array->SetName(soln_name.c_str());
    particles_array->SetNumberOfTuples(num_particles);

    // Loop particles  TODO -- partition particles
    for (std::size_t p = 0; p < num_particles; ++p) {
        // Grab particle and solution
        const auto& particle = particles_[p];
        const auto soln = particle->GetSolution();

        // Save it
        particles_array->SetValue(p, soln);
    }

    // Attach array to grid
    point_data->AddArray(particles_array);

    // ------------------------------------------------------------------------
    // Write
    // ------------------------------------------------------------------------

    // Create writer object
    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetDataModeToBinary();

    // Add grid to writer
    writer->SetInputData(grid_);

    // Write
    const int error = writer->Write();

    // Print error if write failed
    if (error == 0) {
        std::stringstream error_message;
        error_message << "`writer->Write()` failed in "
                      << std::source_location::current().function_name();
        Utilities::PrintErrorOnRoot(error_message.str());
    }

}  // VTKWriterParticles::Write_

// ----------------------------------------------------------------------------
// Write parallel output file
// ----------------------------------------------------------------------------
void VTKWriterParticles::WriteParallel_(
    const std::string& filename,
    const std::vector<std::string>& piece_filenames) const {
    // Open file
    std::ofstream out(filename);

    // Print error if open file failed
    if (!out) {
        std::stringstream error_message;
        error_message << "`std::ofstream out(filename)` failed in "
                      << std::source_location::current().function_name();
        Utilities::PrintErrorOnRoot(error_message.str());
    }

    // Header
    out << R"(<?xml version="1.0"?>)" << "\n";
    out << R"(<VTKFile type="PUnstructuredGrid" )";
    out << R"(version="0.1" byte_order="LittleEndian">)" << "\n";

    // Start PUnstructuredGrid
    out << R"(  <PUnstructuredGrid>)" << "\n";

    // Default type for vtkSmartPointer<vtkPoints>::New() is VTK_FLOAT; can
    // check using vtk::Points::GetDataType() and definitions in vtkSetGet.h
    const char* pointType = "Float32";

    // PPoints declaration (must match the Point in each .vtu)
    out << R"(    <PPoints>)" << "\n";
    out << "      <PDataArray type=\"" << pointType
        << "\" NumberOfComponents=\"3\" Name=\"Points\"/>\n";
    out << R"(    </PPoints>)" << "\n";

    // "Float64" for vtkDoubleArray
    const char* particlesType = "Float64";

    // PPointData: declare the field array written in Write_
    out << R"(    <PPointData>)" << "\n";
    out << "      <PDataArray type=\"" << particlesType
        << "\" NumberOfComponents=\"1\" Name=\"" << particles_name_ << "\"/>\n";
    out << R"(    </PPointData>)" << "\n";

    // PCellData is empty (for now)
    out << R"(    <PCellData>)" << "\n";
    out << R"(    </PCellData>)" << "\n";

    // List pieces
    for (const auto& pname : piece_filenames) {
        out << "    <Piece Source=\"" << pname << "\"/>\n";
    }

    // Footer
    out << R"(  </PUnstructuredGrid>)" << "\n";
    out << R"(</VTKFile>)" << "\n";

    // Close file
    out.close();

}  // VTKWriterParticles::WriteParallel_

// ----------------------------------------------------------------------------
// Write time output file
// ----------------------------------------------------------------------------
void VTKWriterParticles::WriteTime_(
    const std::string& filename,
    const std::vector<std::string>& piece_filenames) const {
    // Open file
    std::ofstream out(filename);

    // Print error if open file failed
    if (!out) {
        std::stringstream error_message;
        error_message << "`std::ofstream out(filename)` failed in "
                      << std::source_location::current().function_name();
        Utilities::PrintErrorOnRoot(error_message.str());
    }

    // Header
    out << R"(<?xml version="1.0"?>)" << "\n";
    out << R"(<VTKFile type="Collection" )";
    out << R"(version="0.1" byte_order="LittleEndian">)" << "\n";

    // Start collection
    out << R"(  <Collection>)" << "\n";

    // List pieces
    for (std::size_t p = 0; p < piece_filenames.size(); ++p) {
        const auto& pname = piece_filenames[p];
        out << "    <DataSet timestep=\"" << p
            << "\" group=\"\" part=\"0\" file=\"" << pname << "\"/>\n";
    }

    // Footer
    out << R"(  </Collection>)" << "\n";
    out << R"(</VTKFile>)" << "\n";

    // Close file
    out.close();

}  // VTKWriterParticles::WriteTime_

}  // namespace pwr
