#ifndef VTK_WRITER_BASE_H
#define VTK_WRITER_BASE_H

#include <string>
#include <vector>

namespace pwr {

// Foward declare Field
// This is sufficient since the header only stores a pointer-like object
// and not the complete type
class Field;

// VTK writer base class
class VTKWriterBase {
   public:
    virtual ~VTKWriterBase() = default;

    // Add field
    virtual void AddField(const std::shared_ptr<pwr::Field>& field,
                          const std::string& field_name) = 0;

    // Write output file
    virtual void Write(const std::string& filename) const = 0;

    // Write parallel output file
    virtual void WriteParallel(
        const std::string& filename,
        const std::vector<std::string>& piece_filenames) const = 0;
};

}  // namespace pwr

#endif  // VTK_WRITER_BASE_H
