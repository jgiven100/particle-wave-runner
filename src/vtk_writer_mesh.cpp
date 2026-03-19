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

#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

#include "mesh_base.h"
#include "mpi_utilities.h"

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

    // Grab active nodes
    const auto &active_nodes = mesh_->GetActiveNodes();

    // Set the number of nodes in global mesh
    const std::size_t num_nodes = active_nodes.size();

    // Grab nodal coordinates
    const auto &nodal_coords = mesh_->GetNodalCoordinates();

    // Create vtk points object
    points_ = vtkSmartPointer<vtkPoints>::New();
    points_->SetNumberOfPoints(num_nodes);

    // Loop nodes
    for (std::size_t n = 0; n < num_nodes; ++n) {
        // Go to next node if this node is not active
        if (active_nodes[n] == 0) {
            continue;
        }

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
    const auto &conn = mesh_->GetElemConnectivity();

    // Set the number of elements
    const std::size_t num_elem_total = mesh_->GetNumElemTotal();

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
            hex->GetPointIds()->SetId(n, conn[8 * e + n]);
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
// Write output files
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
    const int rank = pwr::MPIUtilities::Rank();

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
    writer->Write();  // TODO add check: if (writer->Write() == 0) { ERROR; }

}  // VTKWriterMesh::Write_

}  // namespace pwr
