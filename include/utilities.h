#ifndef PWR_UTILITIES_H
#define PWR_UTILITIES_H

#include <cassert>
#include <iostream>
#include <string>
#include <type_traits>

#include "mpi_utilities.h"

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

    // ------------------------------------------------------------------------
    // Pad string
    // ------------------------------------------------------------------------
    static void PadString(std::string &str, const int width,
                          const char c = '0') {
        // Padding width should be greater or equal to string size
        if (width <= static_cast<int>(str.size())) {
            return;
        }

        // Pad char `c` on left
        str.insert(0, width - str.size(), c);

    }  // Utilities::PadString

    // ------------------------------------------------------------------------
    // Print error on root
    // ------------------------------------------------------------------------
    static void PrintErrorOnRoot(const std::string &str, const int root = 0) {
        // Only print on root (rank 0)
        if (pwr::MPIUtilities::Rank() == root) {
            std::cerr << "\x1b[1;31m[ERROR]\x1b[0m " << str << std::endl;
        }

    }  // PrintErrorOnRoot

    // ------------------------------------------------------------------------
    // Print info on root
    // ------------------------------------------------------------------------
    static void PrintInfoOnRoot(const std::string &str, const int root = 0) {
        // Only print on root (rank 0)
        if (pwr::MPIUtilities::Rank() == root) {
            std::cout << "\x1b[1;32m[INFO]\x1b[0m " << str << std::endl;
        }

    }  // PrintInfoOnRoot

    // ------------------------------------------------------------------------
    // Print warning on root
    // ------------------------------------------------------------------------
    static void PrintWarningOnRoot(const std::string &str, const int root = 0) {
        // Only print on root (rank 0)
        if (pwr::MPIUtilities::Rank() == root) {
            std::cout << "\x1b[1;33m[WARNING]\x1b[0m " << str << std::endl;
        }

    }  // PrintWarningOnRoot
};

}  // namespace pwr

#endif  // PWR_UTILITIES_H
