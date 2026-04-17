#ifndef PWR_FIELD_H
#define PWR_FIELD_H

#include <cassert>
#include <cstddef>
#include <vector>

namespace pwr {

// Field class
class Field {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    Field(std::size_t size) : field_(size) { Setup_(); }

    // ------------------------------------------------------------------------
    // Destructor
    // ------------------------------------------------------------------------
    ~Field() = default;

    // ------------------------------------------------------------------------
    // Delete copy constructor
    // ------------------------------------------------------------------------
    Field(const Field&) = delete;

    // ------------------------------------------------------------------------
    // Delete copy assignment
    // ------------------------------------------------------------------------
    Field& operator=(const Field&) = delete;

    // ------------------------------------------------------------------------
    // Overload [] operator
    // ------------------------------------------------------------------------
    double& operator[](std::size_t i) { return field_[i]; }
    const double& operator[](std::size_t i) const { return field_[i]; }

    // ------------------------------------------------------------------------
    // Getters
    // ------------------------------------------------------------------------

    // Get field
    const std::vector<double>& GetField() const {
        assert(setup_complete_);
        return field_;
    }

    // Get field size
    std::size_t GetFieldSize() const {
        assert(setup_complete_);
        return field_.size();
    }

   private:
    // ------------------------------------------------------------------------
    // Setup field
    // Calls:
    //   TODO
    //   CheckSetup_()
    // ------------------------------------------------------------------------
    void Setup_();

    // ------------------------------------------------------------------------
    // Check field setup
    // Making it this far means that `setup_complete_` should be set to `true`
    // ------------------------------------------------------------------------
    void CheckSetup_();

    // CheckSetup_() has been successfully called
    bool setup_complete_ = false;

    // Field
    std::vector<double> field_;
};

}  // namespace pwr

#endif  // PWR_FIELD_H
