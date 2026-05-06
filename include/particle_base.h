#ifndef PWR_PARTICLE_BASE_H
#define PWR_PARTICLE_BASE_H

#include <array>
#include <cstddef>

namespace pwr {

// Particle bass class
class ParticleBase {
   public:
    virtual ~ParticleBase() = default;

    // ------------------------------------------------------------------------
    // Getters
    // ------------------------------------------------------------------------

    // Get id
    virtual std::size_t GetId() const = 0;

    // Get global coordinates
    virtual const std::array<double, 3>& GetCoordsGlobal() const = 0;

    // Get local coordinates
    virtual const std::array<double, 3>& GetCoordsLocal() const = 0;

    // Get solution
    virtual double GetSolution() const = 0;
};

}  // namespace pwr

#endif  // PWR_PARTICLE_BASE_H
