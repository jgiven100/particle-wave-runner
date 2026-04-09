#include "solver_explicit.h"

#include <iostream>  // TODO remove me

namespace pwr {

// ----------------------------------------------------------------------------
// Setup explicit solver
// ----------------------------------------------------------------------------
void SolverExplicit::Setup_() {
    // TODO
    CheckSetup_();

}  // SolverExplicit::Setup_

// ----------------------------------------------------------------------------
// Check explicit solver
// ----------------------------------------------------------------------------
void SolverExplicit::CheckSetup_() {
    // Setup is complete
    setup_complete_ = true;

}  // SolverExplicit::CheckSetup_

// ----------------------------------------------------------------------------
// Step
// ----------------------------------------------------------------------------
void SolverExplicit::Step_() {
    // TODO
    std::cout << "Inside SolverExplicit::Step_() " << std::endl;

}  // SolverExplicit::Step_

}  // namespace pwr
