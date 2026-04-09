#ifndef PWR_SOLVER_EXPLICIT_H
#define PWR_SOLVER_EXPLICIT_H

#include <cassert>

#include "solver_base.h"

namespace pwr {

// Explicit solver class
class SolverExplicit : public SolverBase {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    SolverExplicit() { Setup_(); }

    // ------------------------------------------------------------------------
    // Destructor
    // ------------------------------------------------------------------------
    ~SolverExplicit() override = default;

    // ------------------------------------------------------------------------
    // Delete copy constructor
    // ------------------------------------------------------------------------
    SolverExplicit(const SolverExplicit&) = delete;

    // ------------------------------------------------------------------------
    // Delete copy assignment
    // ------------------------------------------------------------------------
    SolverExplicit& operator=(const SolverExplicit&) = delete;

    // ------------------------------------------------------------------------
    // Step
    // ------------------------------------------------------------------------
    void Step() override {
        assert(setup_complete_);
        Step_();
    }

   private:
    // ------------------------------------------------------------------------
    // Setup explicit solver
    // Calls:
    //   TODO
    //   CheckSetup_()
    // ------------------------------------------------------------------------
    void Setup_();

    // ------------------------------------------------------------------------
    // Check explicit solver setup
    // Making it this far means that `setup_complete_` should be set to `true`
    // ------------------------------------------------------------------------
    void CheckSetup_();

    // ------------------------------------------------------------------------
    // Step
    // ------------------------------------------------------------------------
    void Step_();

    // CheckSetup_() has been successfully called
    bool setup_complete_ = false;
};

}  // namespace pwr

#endif  // PWR_SOLVER_EXPLICIT_H
