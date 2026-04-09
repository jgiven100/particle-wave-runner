#ifndef PWR_SOLVER_IMPLICIT_H
#define PWR_SOLVER_IMPLICIT_H

#include <cassert>

#include "solver_base.h"

namespace pwr {

// Implicit solver class
class SolverImplicit : public SolverBase {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    SolverImplicit() { Setup_(); }

    // ------------------------------------------------------------------------
    // Destructor
    // ------------------------------------------------------------------------
    ~SolverImplicit() override = default;

    // ------------------------------------------------------------------------
    // Delete copy constructor
    // ------------------------------------------------------------------------
    SolverImplicit(const SolverImplicit&) = delete;

    // ------------------------------------------------------------------------
    // Delete copy assignment
    // ------------------------------------------------------------------------
    SolverImplicit& operator=(const SolverImplicit&) = delete;

    // ------------------------------------------------------------------------
    // Step
    // ------------------------------------------------------------------------
    void Step() override {
        assert(setup_complete_);
        Step_();
    }

   private:
    // ------------------------------------------------------------------------
    // Setup implicit solver
    // Calls:
    //   TODO
    //   CheckSetup_()
    // ------------------------------------------------------------------------
    void Setup_();

    // ------------------------------------------------------------------------
    // Check implicit solver setup
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

#endif  // PWR_SOLVER_IMPLICIT_H
