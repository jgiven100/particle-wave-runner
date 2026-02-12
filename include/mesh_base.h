#ifndef PWR_MESH_BASE_H
#define PWR_MESH_BASE_H

#include <cstddef>

namespace pwr {

// Mesh base class
class MeshBase {
   public:
    virtual ~MeshBase() = default;

    // ------------------------------------------------------------------------
    // Getters
    // ------------------------------------------------------------------------

    // Get number of elements partition
    virtual std::size_t GetNumElemPartition() const = 0;

    // Get number of elements ghost
    virtual std::size_t GetNumElemGhost() const = 0;
};
}  // namespace pwr

#endif  // PWR_MESH_BASE_H
