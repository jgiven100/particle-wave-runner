#include "particle_poisson.h"

#include <cassert>

namespace pwr {

// ----------------------------------------------------------------------------
// Setup particle
// ----------------------------------------------------------------------------
void ParticlePoisson::Setup_() {
    SetContainingElement_();
    CheckSetup_();

}  // ParticlePoisson::Setup_

// ----------------------------------------------------------------------------
// Set containing element
// ----------------------------------------------------------------------------
void ParticlePoisson::SetContainingElement_() {
    eid_local_ = 0;   // mesh->GetContainingElementLocal(coords_global_);
    eid_global_ = 0;  // mesh->GetContainingElementGlobal(coords_global_);
}

// ----------------------------------------------------------------------------
// Check particle
// ----------------------------------------------------------------------------
void ParticlePoisson::CheckSetup_() {
    // Setup is complete
    setup_complete_ = true;

}  // ParticlePoisson::CheckSetup_

}  // namespace pwr
