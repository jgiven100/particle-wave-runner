#ifndef PWR_SOLVER_EXPLICIT_H
#define PWR_SOLVER_EXPLICIT_H

#include <cassert>
#include <memory>
#include <string>
#include <vector>

#include "field.h"
#include "mesh.h"
#include "output_manager.h"
#include "particle_base.h"
#include "solver_base.h"

namespace pwr {

// Forward declare MeshBase
// This is sufficient since the header only stores a pointer-like object
// and the complete type
class MeshBase;

// Explicit solver class
class SolverExplicit : public SolverBase {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    SolverExplicit(
        const std::shared_ptr<const pwr::MeshBase>& mesh,
        const std::vector<std::shared_ptr<pwr::ParticleBase>>& particles,
        const std::vector<std::shared_ptr<pwr::Field>>& fields,
        const std::vector<std::string>& fields_names)
        : mesh_(mesh),
          particles_(particles),
          fields_(fields),
          fields_names_(fields_names) {
        // Sanity check: pointer to mesh is not null
        assert(mesh_);

        Setup_();
    }

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
    // Run
    // ------------------------------------------------------------------------
    void Run() override {
        assert(setup_complete_);
        Run_();
    }

   private:
    // ------------------------------------------------------------------------
    // Setup explicit solver
    // Calls:
    //   SetOutputManager_()
    //   CheckSetup_()
    // ------------------------------------------------------------------------
    void Setup_();

    // ------------------------------------------------------------------------
    // Set output manager
    // Create shared pointer to OutputManager
    // ------------------------------------------------------------------------
    void SetOutputManager_();

    // ------------------------------------------------------------------------
    // Check explicit solver setup
    // Making it this far means that `setup_complete_` should be set to `true`
    // ------------------------------------------------------------------------
    void CheckSetup_();

    // ------------------------------------------------------------------------
    // Run
    // ------------------------------------------------------------------------
    void Run_();

    // CheckSetup_() has been successfully called
    bool setup_complete_ = false;

    // Shared pointer to mesh object
    const std::shared_ptr<const pwr::MeshBase> mesh_;

    // Particles
    const std::vector<std::shared_ptr<pwr::ParticleBase>> particles_;

    // Fields
    const std::vector<std::shared_ptr<pwr::Field>> fields_;

    // Fields names
    const std::vector<std::string> fields_names_;

    // Shared pointer to output manager object
    std::shared_ptr<pwr::OutputManager> output_manager_;

    // Max step
    std::size_t max_step_ = 100;  // TODO : read from file
};

}  // namespace pwr

#endif  // PWR_SOLVER_EXPLICIT_H
