#include "solver_explicit.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <sstream>
#include <vector>

#include "utilities.h"

namespace pwr {

// ----------------------------------------------------------------------------
// Setup explicit solver
// ----------------------------------------------------------------------------
void SolverExplicit::Setup_() {
    SetTimeStepper_();
    SetOutputManager_();
    CheckSetup_();

}  // SolverExplicit::Setup_

// ----------------------------------------------------------------------------
// Set time stepper
// ----------------------------------------------------------------------------
void SolverExplicit::SetTimeStepper_() {
    // Start time
    const double t0 = 0.;  // TODO : read from file

    // Final time
    const double tf = 1.;  // TODO : read from file

    // Sanity check: final time is after start time
    assert(tf > t0);

    // Time step
    const double dt = 1.e-2;  // TODO : read from file

    // Sanity check: time step is less than overall time
    assert(dt < (tf - t0));

    // Compute max steps (round up to nearest integer)
    max_step_ = std::ceil((tf - t0) / dt);

    // Output step
    output_step_ = 10;  // TODO : read from file

    // Sanity check:
    assert(output_step_ <= max_step_);

}  // SolverExplicit::SetTimeStepper_

// ----------------------------------------------------------------------------
// Set output manager
// ----------------------------------------------------------------------------
void SolverExplicit::SetOutputManager_() {
    // Create shared pointer to output manager object
    output_manager_ =
        std::make_shared<pwr::OutputManager>(mesh_, max_step_, output_step_);

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

    // Add particles
    output_manager_->AddParticles(particles_, particles_name_);

    // Finalize setup for writers
    output_manager_->FinalizeWritersSetup();

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
    // Only need to write mesh output once (at the start)
    output_manager_->WriteMeshFiles(0);

    // Output field and particles at time zero
    output_manager_->WriteFieldFiles(0);
    output_manager_->WriteParticlesFiles(0);

    // Loop steps
    for (std::size_t step = 0; step < max_step_; ++step) {
        // TODO -- explicit solve and update solution

        // Output is at the end of the iteration, so we're actually
        // at (step + 1) iterations
        if ((step + 1) % output_step_ == 0) {
            const std::size_t width = Utilities::CountDigits(max_step_);
            std::string step_str = std::to_string(step + 1);
            Utilities::PadString(step_str, width);

            std::stringstream info_message;
            info_message << "Step " << step_str << " of " << max_step_;
            Utilities::PrintInfoOnRoot(info_message.str());

            // output field and particles at current step
            output_manager_->WriteFieldFiles(step + 1);
            output_manager_->WriteParticlesFiles(step + 1);
        }
    }

    // TODO -- write time series

}  // SolverExplicit::Run_

}  // namespace pwr
