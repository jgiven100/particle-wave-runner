#ifndef PWR_MPI_ENVIRONMENT_H
#define PWR_MPI_ENVIRONMENT_H

#include <mpi.h>

namespace pwr {

// MPI Environment class
// RAII wrapper for MPI_Init and MPI_Finalize
class MPIEnvironment {
   public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    MPIEnvironment(int& argc, char**& argv) {
        int initialized;
        MPI_Initialized(&initialized);

        if (!initialized) {
            MPI_Init(&argc, &argv);
        }
    }

    // ------------------------------------------------------------------------
    // Destructor
    // ------------------------------------------------------------------------
    ~MPIEnvironment() noexcept {
        int finalized;
        MPI_Finalized(&finalized);

        if (!finalized) {
            MPI_Finalize();
        }
    }

    // ------------------------------------------------------------------------
    // Delete copy constructor
    // ------------------------------------------------------------------------
    MPIEnvironment(const MPIEnvironment&) = delete;

    // ------------------------------------------------------------------------
    // Delete copy assignment
    // ------------------------------------------------------------------------
    MPIEnvironment& operator=(const MPIEnvironment&) = delete;
};

}  // namespace pwr

#endif  // PWR_MPI_ENVIRONMENT_H
