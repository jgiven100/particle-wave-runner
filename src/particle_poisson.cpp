#include "particle_poisson.h"

#include <cassert>

namespace pwr {

// ----------------------------------------------------------------------------
// Setup particle
// ----------------------------------------------------------------------------
void ParticlePoisson::Setup_() {
    // TODO
    CheckSetup_();

}  // ParticlePoisson::Setup_

// ----------------------------------------------------------------------------
// Check particle
// ----------------------------------------------------------------------------
void ParticlePoisson::CheckSetup_() {
    // Setup is complete
    setup_complete_ = true;

}  // ParticlePoisson::CheckSetup_

}  // namespace pwr
