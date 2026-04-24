#ifndef PWR_PARTICLE_POISSON_H
#define PWR_PARTICLE_POISSON_H

#include <cassert>
#include <cstddef>
#include <iostream>  // TODO

#include "particle_base.h"

namespace pwr {

// Particle class for Poisson problem
class ParticlePoisson : public ParticleBase {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    ParticlePoisson() {
        std::cout << "Inside ParticlePoisson constructor" << std::endl;
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
};

}  // namespace pwr

#endif  // PWR_PARTICLE_POISSON_H
