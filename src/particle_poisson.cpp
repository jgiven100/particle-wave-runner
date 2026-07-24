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
    // Grab global element id from mesh
    eid_global_ = mesh_->FindContainingElemIdGlobal(coords_global_);

    // TODO
    eid_local_ = 0;  // mesh->GetContainingElementLocal(coords_global_);
}

// ----------------------------------------------------------------------------
// Check particle
// ----------------------------------------------------------------------------
void ParticlePoisson::CheckSetup_() {
    // Setup is complete
    setup_complete_ = true;

}  // ParticlePoisson::CheckSetup_

}  // namespace pwr
