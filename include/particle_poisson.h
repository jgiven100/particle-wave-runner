#ifndef PWR_PARTICLE_POISSON_H
#define PWR_PARTICLE_POISSON_H

#include <array>
#include <cassert>
#include <cstddef>

#include "particle_base.h"

namespace pwr {

// Particle class for Poisson problem
class ParticlePoisson : public ParticleBase {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    ParticlePoisson(std::size_t id, std::array<double, 3> coords_global)
        : id_(id), coords_global_(coords_global) {
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

    // Get id
    std::size_t GetId() const override {
        assert(setup_complete_);
        return id_;
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

   private:
    // ------------------------------------------------------------------------
    // Setup particle
    // Calls:
    //   TODO
    //   CheckSetup_()
    // ------------------------------------------------------------------------
    void Setup_();

    // ------------------------------------------------------------------------
    // Check setup
    // Making it this far means that `setup_complete_` should be set to `true`
    // ------------------------------------------------------------------------
    void CheckSetup_();

    // CheckSetup_() has been successfully called
    bool setup_complete_ = false;

    // Id
    const std::size_t id_;

    // Local coordinates
    std::array<double, 3> coords_local_;

    // Global coordinates
    std::array<double, 3> coords_global_;

    // Solution
    double u_ = 101.;  // TODO update
};

}  // namespace pwr

#endif  // PWR_PARTICLE_POISSON_H
