#ifndef PWR_SOLVER_BASE_H
#define PWR_SOLVER_BASE_H

namespace pwr {

// Solver base class
class SolverBase {
   public:
    virtual ~SolverBase() = default;

    // Step
    virtual void Step() = 0;
};

}  // namespace pwr

#endif  // PWR_SOLVER_BASE_H
