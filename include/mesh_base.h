#ifndef PWR_MESH_BASE_H
#define PWR_MESH_BASE_H

#include <cstddef>
#include <unordered_map>
#include <vector>

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

    // Get number of elements partition + ghost
    virtual std::size_t GetNumElemTotal() const = 0;

    // Get number of active nodes
    virtual std::size_t GetNumNodesActive() const = 0;

    // Get global element indices (in each direction)
    virtual const std::vector<std::size_t>& GetElemIndexGlobal() const = 0;

    // Get global element id
    virtual const std::vector<std::size_t>& GetElemIdGlobal() const = 0;

    // Get global element-wise connectivity
    virtual const std::vector<std::size_t>& GetElemConnGlobal() const = 0;

    // Get local element-wise connectivity
    virtual const std::vector<std::size_t>& GetElemConnLocal() const = 0;

    // Get nodal ownership
    virtual const std::vector<char>& GetNodalOwnership() const = 0;

    // Get nodal coordinates
    virtual const std::vector<double>& GetNodalCoordinates() const = 0;

    // ------------------------------------------------------------------------
    // Find containing global element id
    // ------------------------------------------------------------------------
    virtual std::size_t FindContainingElemIdGlobal(
        const std::vector<double>& coords) const = 0;
};

}  // namespace pwr

#endif  // PWR_MESH_BASE_H
