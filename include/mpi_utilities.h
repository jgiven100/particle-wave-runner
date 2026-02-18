#ifndef PWR_MPI_UTILITIES_H
#define PWR_MPI_UTILITIES_H

#include <mpi.h>

namespace pwr {

// MPI Utility class
// Place for useful MPI functions to simplify remainder of codebase
class MPIUtilities {
   public:
    // ------------------------------------------------------------------------
    // Getters
    // ------------------------------------------------------------------------

    // Get MPI rank
    static int Rank(MPI_Comm comm = MPI_COMM_WORLD) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        return rank;
    }

    // Get MPI size
    static int Size(MPI_Comm comm = MPI_COMM_WORLD) {
        int size;
        MPI_Comm_size(comm, &size);
        return size;
    }
};

}  // namespace pwr

#endif  // PWR_MPI_UTILITIES_H
