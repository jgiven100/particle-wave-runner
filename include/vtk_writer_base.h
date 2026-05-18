#ifndef VTK_WRITER_BASE_H
#define VTK_WRITER_BASE_H

#include <memory>
#include <string>
#include <vector>

namespace pwr {

// Foward declare ParticleBase, Field
// This is sufficient since the header only stores a pointer-like object
// and not the complete type
class ParticleBase;
class Field;

// VTK writer base class
class VTKWriterBase {
   public:
    virtual ~VTKWriterBase() = default;

    // Finalize setup
    virtual void FinalizeSetup() = 0;

    // Add field
    virtual void AddField(const std::shared_ptr<pwr::Field>& field,
                          const std::string& field_name) = 0;

    // Add particles
    virtual void AddParticles(
        const std::vector<std::shared_ptr<pwr::ParticleBase>>& particles,
        const std::string& particles_name) = 0;

    // Write output file
    virtual void Write(const std::string& filename) const = 0;

    // Write parallel output file
    virtual void WriteParallel(
        const std::string& filename,
        const std::vector<std::string>& piece_filenames) const = 0;

    // Write time output file
    virtual void WriteTime(
        const std::string& filename,
        const std::vector<std::string>& piece_filenames) const = 0;
};

}  // namespace pwr

#endif  // VTK_WRITER_BASE_H
