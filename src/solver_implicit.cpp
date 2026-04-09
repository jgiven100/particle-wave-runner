#include "solver_implicit.h"

#include <iostream>  // TODO remove me

namespace pwr {

// ----------------------------------------------------------------------------
// Setup implicit solver
// ----------------------------------------------------------------------------
void SolverImplicit::Setup_() {
    // TODO
    CheckSetup_();

}  // SolverImplicit::Setup_

// ----------------------------------------------------------------------------
// Check implicit solver
// ----------------------------------------------------------------------------
void SolverImplicit::CheckSetup_() {
    // Setup is complete
    setup_complete_ = true;

}  // SolverImplicit::CheckSetup_

// ----------------------------------------------------------------------------
// Step
// ----------------------------------------------------------------------------
void SolverImplicit::Step_() {
    // TODO
    std::cout << "Inside SolverImplicit::Step_()" << std::endl;

}  // SolverExplicit::Step_

}  // namespace pwr
