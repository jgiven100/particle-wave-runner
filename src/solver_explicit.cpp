#include "solver_explicit.h"

#include <cassert>
#include <cstddef>
#include <memory>
#include <vector>

namespace pwr {

// ----------------------------------------------------------------------------
// Setup explicit solver
// ----------------------------------------------------------------------------
void SolverExplicit::Setup_() {
    SetOutputManager_();
    CheckSetup_();

}  // SolverExplicit::Setup_

// ----------------------------------------------------------------------------
// Set output manager
// ----------------------------------------------------------------------------
void SolverExplicit::SetOutputManager_() {
    // Create shared pointer to output manager object
    output_manager_ = std::make_shared<pwr::OutputManager>(mesh_, max_step_);

    // Sanity check: fields and fields names are the same size
    assert(fields_.size() == fields_names_.size());

    // Add fields
    for (std::size_t i = 0; i < fields_.size(); ++i) {
        // Grab field
        const auto& field = fields_[i];

        // Grab field name
        const auto& field_name = fields_names_[i];

        // Add field to output manager
        output_manager_->AddField(field, field_name);
    }

}  // SolverExplicit::SetOutputManager_

// ----------------------------------------------------------------------------
// Check explicit solver
// ----------------------------------------------------------------------------
void SolverExplicit::CheckSetup_() {
    // Setup is complete
    setup_complete_ = true;

}  // SolverExplicit::CheckSetup_

// ----------------------------------------------------------------------------
// Run
// ----------------------------------------------------------------------------
void SolverExplicit::Run_() {
    ;  // TODO

    const std::size_t step = 0;
    output_manager_->WriteMeshFiles(step);
    output_manager_->WriteFieldFiles(step);

}  // SolverExplicit::Run_

}  // namespace pwr
