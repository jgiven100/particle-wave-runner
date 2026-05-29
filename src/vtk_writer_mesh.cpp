#include "vtk_writer_mesh.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataSetAttributes.h>
#include <vtkHexahedron.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
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
// Setup vtk writer for mesh
// ----------------------------------------------------------------------------
void VTKWriterMesh::Setup_() {
    InitializeWriter_();
    CheckSetup_();

}  // VTKWriterMesh::Setup_

// ----------------------------------------------------------------------------
// Initialize vtk writer for mesh
// ----------------------------------------------------------------------------
void VTKWriterMesh::InitializeWriter_() {
    // ------------------------------------------------------------------------
    // Assign points
    // ------------------------------------------------------------------------

    // Grab the number of active nodes for this partition
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

    // Set the number of elements
    const std::size_t num_elem_total = mesh_->GetNumElemTotal();

    // Map for right-hand rule ordring (gmsh and vtk)
    std::array<int, 8> vtk_ordering{0, 1, 3, 2, 4, 5, 7, 6};

    // Sanity check: connectivity array has 8 nodes per element
    assert(conn.size() == 8 * num_elem_total);

    // Create vtk grid object
    grid_ = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid_->Allocate(num_elem_total);
    grid_->SetPoints(points_);

    // Loop partition and ghost elements
    for (std::size_t e = 0; e < num_elem_total; ++e) {
        // Create vtk hexahedron object
        auto hex = vtkSmartPointer<vtkHexahedron>::New();

        // Loop nodes per elements
        for (int n = 0; n < 8; ++n) {
            const int index = vtk_ordering[n];
            hex->GetPointIds()->SetId(n, conn[8 * e + index]);
        }

        // Add hex element to grid
        grid_->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
    }

    // Sanity check: pointer to grid object is not null
    assert(grid_);

}  // VTKWriterMesh::InitializeWriter_

// ----------------------------------------------------------------------------
// Check setup
// ----------------------------------------------------------------------------
void VTKWriterMesh::CheckSetup_() {
    // Setup is complete
    setup_complete_ = true;

}  // VTKWriterMesh::CheckSetup_

// ----------------------------------------------------------------------------
// Write output file
// ----------------------------------------------------------------------------
void VTKWriterMesh::Write_(const std::string &filename) const {
    // Clear previous arrays
    auto cell_data = grid_->GetCellData();
    cell_data->Initialize();

    // Set the number of elements
    const std::size_t num_elem_total = mesh_->GetNumElemTotal();

    // ------------------------------------------------------------------------
    // Save rank
    // ------------------------------------------------------------------------

    // Grab rank
    const int rank = MPIUtilities::Rank();

    // Create array to save rank info
    auto rank_array = vtkSmartPointer<vtkIntArray>::New();
    rank_array->SetName("rank");
    rank_array->SetNumberOfTuples(num_elem_total);

    // Loop partition and ghost elements
    for (std::size_t e = 0; e < num_elem_total; ++e) {
        rank_array->SetValue(e, rank);
    }

    // Attach array to grid
    cell_data->AddArray(rank_array);

    // ------------------------------------------------------------------------
    // Save ghost
    // ------------------------------------------------------------------------

    // Set the number of partition elements and ghost elements
    const std::size_t num_elem_partition = mesh_->GetNumElemPartition();
    const std::size_t num_elem_ghost = mesh_->GetNumElemGhost();

    // Sanity check: partition plus ghost elements equals total elements
    assert(num_elem_total == num_elem_partition + num_elem_ghost);

    // Create array to save ghost info
    auto ghost_array = vtkSmartPointer<vtkUnsignedCharArray>::New();
    ghost_array->SetName(vtkDataSetAttributes::GhostArrayName());
    ghost_array->SetNumberOfTuples(num_elem_total);

    // Loop partition elements
    for (std::size_t e = 0; e < num_elem_partition; ++e) {
        ghost_array->SetValue(e, 0);
    }

    // Loop ghost elements
    for (std::size_t e = num_elem_partition; e < num_elem_total; ++e) {
        ghost_array->SetValue(e, vtkDataSetAttributes::DUPLICATECELL);
    }

    // Attach array to grid
    cell_data->AddArray(ghost_array);

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

}  // VTKWriterMesh::Write_

// ----------------------------------------------------------------------------
// Write parallel output file
// ----------------------------------------------------------------------------
void VTKWriterMesh::WriteParallel_(
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
    // check using vtkPoints::GetDataType() and definitions in vtkSetGet.h
    const char *pointType = "Float32";

    // PPoints declaration (must match the Points in each .vtu)
    out << R"(    <PPoints>)" << "\n";
    out << "      <PDataArray type=\"" << pointType
        << "\" NumberOfComponents=\"3\" Name=\"Points\"/>\n";
    out << R"(    </PPoints>)" << "\n";

    // PPointData is empty (for now)
    out << R"(    <PPointData>)" << "\n";
    out << R"(    </PPointData>)" << "\n";

    // "Int32" for vtkIntArray
    const char *rankType = "Int32";

    // "UInt8" for vtkUnsignedCharArray (vtkGhostType)
    const char *ghostType = "UInt8";

    // vtkDataSetAttributes::GhostArrayName()
    const char *ghostName = "vtkGhostType";

    // PCellData: declare the cell arrays written in Write_
    out << R"(    <PCellData>)" << "\n";
    out << "      <PDataArray type=\"" << rankType
        << "\" Name=\"rank\" NumberOfComponents=\"1\"/>\n";
    out << "      <PDataArray type=\"" << ghostType << "\" Name=\"" << ghostName
        << "\" NumberOfComponents=\"1\"/>\n";
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

}  // VTKWriterMesh::WriteParallel_

// ----------------------------------------------------------------------------
// Write time output file
// ----------------------------------------------------------------------------
void VTKWriterMesh::WriteTime_(
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
    out << R"(<VTKFile type="Collection" )";
    out << R"(version="0.1" byte_order="LittleEndian">)" << "\n";

    // Start Collection
    out << R"(  <Collection>)" << "\n";

    // List pieces
    for (std::size_t p = 0; p < piece_filenames.size(); ++p) {
        const auto &pname = piece_filenames[p];
        out << "    <DataSet timestep=\"" << p
            << "\" group=\"\" part=\"0\" file=\"" << pname << "\"/>\n";
    }

    // Footer
    out << R"(  </Collection>)" << "\n";
    out << R"(</VTKFile>)" << "\n";

    // Close file
    out.close();

}  // VTKWriterMesh::WriteTime_

}  // namespace pwr
