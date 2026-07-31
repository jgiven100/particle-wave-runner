#include "particle_poisson.h"

#include <cassert>

namespace pwr {

// ----------------------------------------------------------------------------
// Setup particle
// ----------------------------------------------------------------------------
void ParticlePoisson::Setup_() {
    SetContainingElement_();
    SetLocalCoordinates_();
    CheckSetup_();

}  // ParticlePoisson::Setup_

// ----------------------------------------------------------------------------
// Set containing element
// ----------------------------------------------------------------------------
void ParticlePoisson::SetContainingElement_() {
    // Grab global element id from mesh
    eid_global_ = mesh_->FindContainingElemIdGlobal(coords_global_);

    // Placeholder for local element id (currently not used!)
    eid_local_ = 0;

}  // ParticlePoisson::SetContainingElement_

// ----------------------------------------------------------------------------
// Set local coordinates
// ----------------------------------------------------------------------------
void ParticlePoisson::SetLocalCoordinates_() {
    // Compute local coordinates in mesh (mesh knowns own element shape/type)
    mesh_->ComputeLocalCoordinates(eid_global_, coords_global_, coords_local_);

}  // ParticlePoisson::SetLocalCoordinates_

// ----------------------------------------------------------------------------
// Check particle
// ----------------------------------------------------------------------------
void ParticlePoisson::CheckSetup_() {
    // Setup is complete
    setup_complete_ = true;

}  // ParticlePoisson::CheckSetup_

}  // namespace pwr
