#ifndef VTK_WRITER_BASE_H
#define VTK_WRITER_BASE_H

#include <string>
#include <vector>

namespace pwr {

// VTK writer base class
class VTKWriterBase {
   public:
    virtual ~VTKWriterBase() = default;

    // Write output file
    virtual void Write(const std::string& filename) const = 0;

    // Write parallel output file
    virtual void WriteParallel(
        const std::string& filename,
        const std::vector<std::string>& piece_filenames) const = 0;
};

}  // namespace pwr

#endif  // VTK_WRITER_BASE_H
