#include "particle_poisson.h"

#include <cassert>
#include <unordered_map>

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

    // Grab map from global-to-local element id
    const auto& elem_id_local = mesh_->GetElemIdLocal();

    // Find corresponding local element id
    const auto it = elem_id_local.find(eid_global_);
    assert(it != elem_id_local.end());
    eid_local_ = it->second;

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
