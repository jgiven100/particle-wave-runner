#ifndef PWR_PARTICLE_POISSON_H
#define PWR_PARTICLE_POISSON_H

#include <array>
#include <cassert>
#include <cstddef>
#include <memory>

#include "mesh_base.h"
#include "particle_base.h"

namespace pwr {

// Particle class for Poisson problem
class ParticlePoisson : public ParticleBase {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    ParticlePoisson(std::shared_ptr<const MeshBase> mesh, std::size_t pid_local,
                    std::size_t pid_global, bool owned,
                    const std::array<double, 3>& coords_global)
        : mesh_(mesh),
          pid_local_(pid_local),
          pid_global_(pid_global),
          owned_(owned),
          coords_global_(coords_global) {
        Setup_();
    }

    // ------------------------------------------------------------------------
    // Destructor
    // ------------------------------------------------------------------------
    ~ParticlePoisson() override = default;

    // ------------------------------------------------------------------------
    // Delete copy constructor
    // ------------------------------------------------------------------------
    ParticlePoisson(const ParticlePoisson&) = delete;

    // ------------------------------------------------------------------------
    // Delete copy assignment
    // ------------------------------------------------------------------------
    ParticlePoisson& operator=(const ParticlePoisson&) = delete;

    // ------------------------------------------------------------------------
    // Getters
    // ------------------------------------------------------------------------

    // Get local id (particle)
    std::size_t GetIdLocal() const override {
        assert(setup_complete_);
        return pid_local_;
    }

    // Get global id (particle)
    std::size_t GetIdGlobal() const override {
        assert(setup_complete_);
        return pid_global_;
    }

    // Get global coordinates
    const std::array<double, 3>& GetCoordsGlobal() const override {
        assert(setup_complete_);
        return coords_global_;
    }

    // Get local coordinates
    const std::array<double, 3>& GetCoordsLocal() const override {
        assert(setup_complete_);
        return coords_local_;
    }

    // Get solution
    double GetSolution() const override {
        assert(setup_complete_);
        return u_;
    }

    // ------------------------------------------------------------------------
    // Setters
    // ------------------------------------------------------------------------

    // Set solution
    void SetSolution(double u) override {
        assert(setup_complete_);
        u_ = u;
    }

   private:
    // ------------------------------------------------------------------------
    // Setup particle
    // Calls:
    //   SetContainingElement_() -- TODO
    //   SetLocalCoordinates_() -- TODO
    //   SetConnectedNodes_() -- TODO
    //   ComputeShapeFunctions_() -- TODO
    //   ComputeShapeFunctionGradientsLocal_() -- TODO
    //   ComputeJacobianTerms_() -- TODO
    //   ComputeShapeFunctionGradientsGlobal_() -- TODO
    //   ComputeVolume_() -- TODO
    //   CheckSetup_()
    // ------------------------------------------------------------------------
    void Setup_();

    // ------------------------------------------------------------------------
    // Set containing element
    // ------------------------------------------------------------------------
    void SetContainingElement_();

    // ------------------------------------------------------------------------
    // Compute Jacobian terms
    // Calls:
    //   ComputeJacobian_() -- TODO
    //   ComputeJacobianInverse_() -- TODO
    //   ComputeJacobianDeterminate_() -- TODO
    // ------------------------------------------------------------------------

    // ------------------------------------------------------------------------
    // Check setup
    // Making it this far means that `setup_complete_` should be set to `true`
    // ------------------------------------------------------------------------
    void CheckSetup_();

    // CheckSetup_() has been successfully called
    bool setup_complete_ = false;

    // Mesh
    const std::shared_ptr<const MeshBase> mesh_;

    // Local particle id (partition)
    const std::size_t pid_local_;

    // Gobal particle id
    const std::size_t pid_global_;

    // Owned particle
    const bool owned_;

    // Local coordinates
    std::array<double, 3> coords_local_;

    // Global coordinates
    std::array<double, 3> coords_global_;

    // Local element id (partition) for containing element
    std::size_t eid_local_;

    // Global element id for containing element
    std::size_t eid_global_;

    // Solution
    double u_;
};

}  // namespace pwr

#endif  // PWR_PARTICLE_POISSON_H
