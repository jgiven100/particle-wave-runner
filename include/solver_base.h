#ifndef PWR_SOLVER_BASE_H
#define PWR_SOLVER_BASE_H

namespace pwr {

// Solver base class
class SolverBase {
   public:
    virtual ~SolverBase() = default;

    // Run
    virtual void Run() = 0;
};

}  // namespace pwr

#endif  // PWR_SOLVER_BASE_H
