#include "vtk_writer_field.h"

#include <vtkDoubleArray.h>
#include <vtkHexahedron.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <array>
#include <cassert>
#include <cstddef>
#include <fstream>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>

#include "mesh_base.h"
#include "mpi_utilities.h"
#include "utilities.h"

namespace pwr {

// ----------------------------------------------------------------------------
// Setup vtk writer for field
// ----------------------------------------------------------------------------
void VTKWriterField::Setup_() {
    InitializeWriter_();
    CheckSetup_();

}  // VTKWriterField::Setup_

// ----------------------------------------------------------------------------
// Initialize vtk wrtier for field
// ----------------------------------------------------------------------------
void VTKWriterField::InitializeWriter_() {
    // ------------------------------------------------------------------------
    // Assign points
    // ------------------------------------------------------------------------

    // Grab the number of active nodes for this pratition
    const std::size_t num_nodes_active = mesh_->GetNumNodesActive();

    // Grab nodal coordinates
    const auto &nodal_coords = mesh_->GetNodalCoordinates();

    // Create vtk points object
    points_ = vtkSmartPointer<vtkPoints>::New();
    points_->SetNumberOfPoints(num_nodes_active);

    // Loop nodes
    for (std::size_t n = 0; n < num_nodes_active; ++n) {
        // Pointer to start of nodal coordinates for current node
        const double *coords = &nodal_coords[3 * n];

        // Set coordinates for this point
        points_->SetPoint(n, coords);
    }

    // Sanity check: pointer to points object is not null
    assert(points_);

    // ------------------------------------------------------------------------
    // Create grid and add hex elements
    // ------------------------------------------------------------------------

    // Grab mesh connectivity
    const auto &conn = mesh_->GetElemConnLocal();

    // Sanity check: connectivity array has 8 nodes per element
    assert(conn.size() % 8 == 0);

    // Set the number of elements
    const std::size_t num_elem_partition = mesh_->GetNumElemPartition();

    // Map for right-hand rule ordering (gmsh and vtk)
    std::array<int, 8> vtk_ordering{0, 1, 3, 2, 4, 5, 7, 6};

    // Sanity check: connectivity array has enough nodes for the partition
    // elements (could be more due to ghost elements)
    assert(conn.size() >= 8 * num_elem_partition);

    // Create vtk grid object
    grid_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid_->Allocate(num_elem_partition);
    grid_->SetPoints(points_);

    // Loop partition elements
    for (std::size_t e = 0; e < num_elem_partition; ++e) {
        // Create vtk hexahedron object
        auto hex = vtkSmartPointer<vtkHexahedron>::New();

        // Loop nodes per element
        for (int n = 0; n < 8; ++n) {
            const int index = vtk_ordering[n];
            hex->GetPointIds()->SetId(n, conn[8 * e + index]);
        }

        // Add hex element to grid
        grid_->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
    }

    // Sanity check: pointer to grid object is not null
    assert(grid_);

}  // VTKWriterField::InitializeWriter_

// ----------------------------------------------------------------------------
// Check setup
// ----------------------------------------------------------------------------
void VTKWriterField::CheckSetup_() {
    // Setup is complete
    setup_complete_ = true;

}  // VTKWriterField::CheckSetup_

// ----------------------------------------------------------------------------
// Write output files
// ----------------------------------------------------------------------------
void VTKWriterField::Write_(const std::string &filename) const {
    // Clear previous arrays
    auto point_data = grid_->GetPointData();
    point_data->Initialize();

    // Set the number of active nodes
    const std::size_t num_nodes_active = mesh_->GetNumNodesActive();

    // ------------------------------------------------------------------------
    // Save field
    // ------------------------------------------------------------------------

    // Sanity check: number of fields and names are the same
    assert(fields_.size() == fields_names_.size());

    for (std::size_t i = 0; i < fields_.size(); ++i) {
        // Set field and field name
        const auto &field = fields_[i];
        const auto &field_name = fields_names_[i];

        // Sanity check: field is the same size as number of active nodes
        assert(num_nodes_active == field->GetFieldSize());

        // Create array to save field
        auto field_array = vtkSmartPointer<vtkDoubleArray>::New();
        field_array->SetName(field_name.c_str());
        field_array->SetNumberOfTuples(num_nodes_active);

        // Grab vector
        const auto &f = field->GetField();

        // Loop active nodes
        for (std::size_t n = 0; n < num_nodes_active; ++n) {
            field_array->SetValue(n, f[n]);
        }

        // Attach array to grid
        point_data->AddArray(field_array);
    }

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

}  // VTKWriterField::Write_

// ----------------------------------------------------------------------------
// Writer parallel output files
// ----------------------------------------------------------------------------
void VTKWriterField::WriteParallel_(
    const std::string &filename,
    const std::vector<std::string> &piece_filenames) const {
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
    const char *pointType = "Float32";

    // PPoints declaration (must match the Points in each .vtu)
    out << R"(    <PPoints>)" << "\n";
    out << "      <PDataArray type=\"" << pointType
        << "\" NumberOfComponents=\"3\" Name=\"Points\"/>\n";
    out << R"(    </PPoints>)" << "\n";

    // "Float64" for vtkDoubleArray
    const char *fieldType = "Float64";

    // PPointData: declare the field arrays written in Wrtie_
    out << R"(    <PPointData>)" << "\n";

    for (const auto &field_name : fields_names_) {
        out << "      <PDataArray type=\"" << fieldType
            << "\" NumberOfComponents=\"1\" Name=\"" << field_name << "\"/>\n";
    }

    out << R"(    </PPointData>)" << "\n";

    // PCellData is empty (for now)
    out << R"(    <PCellData>)" << "\n";
    out << R"(    </PCellData>)" << "\n";

    // List pieces
    for (const auto &pname : piece_filenames) {
        out << "    <Piece Source=\"" << pname << "\"/>\n";
    }

    // Footer
    out << R"(  </PUnstructuredGrid>)" << "\n";
    out << R"(</VTKFile>)" << "\n";

    // Close file
    out.close();

}  // VTKWriterField::WriteParallel_

}  // namespace pwr
