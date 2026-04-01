#ifndef PWR_MPI_UTILITIES_H
#define PWR_MPI_UTILITIES_H

#include <mpi.h>

namespace pwr {

template <typename T>
struct MPIType {
    static_assert(sizeof(T) == 0, "Unsupported MPI datatype");
};

// MPI Utility class
// Place for useful MPI functions to simplify remainder of codebase
class MPIUtilities {
   public:
    // ------------------------------------------------------------------------
    // Get MPI rank
    // ------------------------------------------------------------------------
    static int Rank(MPI_Comm comm = MPI_COMM_WORLD) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        return rank;
    }

    // ------------------------------------------------------------------------
    // Get MPI size
    // ------------------------------------------------------------------------
    static int Size(MPI_Comm comm = MPI_COMM_WORLD) {
        int size;
        MPI_Comm_size(comm, &size);
        return size;
    }

    // ------------------------------------------------------------------------
    // MPI barrier
    // ------------------------------------------------------------------------
    static void Barrier(MPI_Comm comm = MPI_COMM_WORLD) { MPI_Barrier(comm); }

    // ------------------------------------------------------------------------
    // Get MPI type
    // ------------------------------------------------------------------------
    template <typename T>
    static MPI_Datatype Type() {
        return MPIType<T>::value;
    }

    // ------------------------------------------------------------------------
    // MPI all reduce
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllReduce(const T& sendbuf, T& recvbuf, MPI_Op op,
                          MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce maximum
    // ------------------------------------------------------------------------
    template <typename T>
    static T AllReduceMax(const T& sendbuf, MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce minimum
    // ------------------------------------------------------------------------
    template <typename T>
    static T AllReduceMin(const T& sendbuf, MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce sum
    // ------------------------------------------------------------------------
    template <typename T>
    static T AllReduceSum(const T& sendbuf, MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce in place
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllReduceInPlace(T& sendbuf, MPI_Op op,
                                 MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce maximum in place
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllReduceMaxInPlace(T& sendbuf, MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce minimum in place
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllReduceMinInPlace(T& sendbuf, MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce sum in place
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllReduceSumInPlace(T& sendbuf, MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI reduce
    // ------------------------------------------------------------------------
    template <typename T>
    static void Reduce(const T& sendbuf, T& recvbuf, MPI_Op op, int root,
                       MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI reduce maximum
    // ------------------------------------------------------------------------
    template <typename T>
    static void ReduceMax(const T& sendbuf, T& recvbuf, int root,
                          MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI reduce maximum
    // ------------------------------------------------------------------------
    template <typename T>
    static void ReduceMin(const T& sendbuf, T& recvbuf, int root,
                          MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI reduce sum
    // ------------------------------------------------------------------------
    template <typename T>
    static void ReduceSum(const T& sendbuf, T& recvbuf, int root,
                          MPI_Comm comm = MPI_COMM_WORLD);
};

}  // namespace pwr

#include "mpi_utilities.tcc"

#endif  // PWR_MPI_UTILITIES_H
