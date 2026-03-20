#ifndef PWR_UTILITIES_H
#define PWR_UTILITIES_H

#include <type_traits>

namespace pwr {

// Utility class
class Utilities {
   public:
    // ------------------------------------------------------------------------
    // Count digits in number
    // ------------------------------------------------------------------------
    template <typename T>
    static int CountDigits(T n) {
        // Sanity check: type is integer-like
        static_assert(std::is_integral_v<T>,
                      "CountDigits requires an integral type");

        // Early exit for 0
        if (n == 0) {
            return 1;
        }

        // Cast into unsigned and handle negative values
        using UnsignedT = std::make_unsigned_t<T>;
        UnsignedT n_unsigned = (n < 0) ? static_cast<UnsignedT>(-(n + 1)) + 1
                                       : static_cast<UnsignedT>(n);

        // Count digits
        int digits = 1;
        while (n_unsigned >= 10u) {
            n_unsigned /= 10u;
            ++digits;
        }

        return digits;

    }  // Utilities::CountDigits
};

}  // namespace pwr

#endif  // PWR_UTILITIES_H
