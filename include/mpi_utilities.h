#ifndef PWR_MPI_UTILITIES_H
#define PWR_MPI_UTILITIES_H

#include <mpi.h>

#include <vector>

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
        return MPIType<T>::type;
    }

    // ------------------------------------------------------------------------
    // MPI all reduce (scalar)
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllReduce(const T& sendbuf, T& recvbuf, MPI_Op op,
                          MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce (vector)
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllReduce(const std::vector<T>& sendbuf,
                          std::vector<T>& recvbuf, MPI_Op op,
                          MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce (scalar) maximum
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllReduceMax(const T& sendbuf, T& recvbuf,
                             MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce (vector) maximum
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllReduceMax(const std::vector<T>& sendbuf,
                             std::vector<T>& recvbuf,
                             MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce (scalar) minimum
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllReduceMin(const T& sendbuf, T& recvbuf,
                             MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce (vector) minimum
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllReduceMin(const std::vector<T>& sendbuf,
                             std::vector<T>& recvbuf,
                             MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce (scalar) sum
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllReduceSum(const T& sendbuf, T& recvbuf,
                             MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all reduce (vector) sum
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllReduceSum(const std::vector<T>& sendbuf,
                             std::vector<T>& recvbuf,
                             MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI reduce (scalar)
    // ------------------------------------------------------------------------
    template <typename T>
    static void Reduce(const T& sendbuf, T& recvbuf, MPI_Op op, int root,
                       MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI reduce (scalar) maximum
    // ------------------------------------------------------------------------
    template <typename T>
    static void ReduceMax(const T& sendbuf, T& recvbuf, int root,
                          MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI reduce (scalar) maximum
    // ------------------------------------------------------------------------
    template <typename T>
    static void ReduceMin(const T& sendbuf, T& recvbuf, int root,
                          MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI reduce (scalar) sum
    // ------------------------------------------------------------------------
    template <typename T>
    static void ReduceSum(const T& sendbuf, T& recvbuf, int root,
                          MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all gather (scalar)
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllGather(const T& sendbuf, std::vector<T>& recvbuf,
                          MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all gather (equal-sized vector)
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllGather(const std::vector<T>& sendbuf,
                          std::vector<T>& recvbuf,
                          MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI all gather (variable-sized vector)
    // ------------------------------------------------------------------------
    template <typename T>
    static void AllGatherV(const std::vector<T>& sendbuf,
                           std::vector<T>& recvbuf,
                           MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI gather (scalar)
    // ------------------------------------------------------------------------
    template <typename T>
    static void Gather(const T& sendbuf, std::vector<T>& recvbuf, int root,
                       MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI gather (equal-sized vector)
    // ------------------------------------------------------------------------
    template <typename T>
    static void Gather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                       int root, MPI_Comm comm = MPI_COMM_WORLD);

    // ------------------------------------------------------------------------
    // MPI gather (variable-sized vector)
    // ------------------------------------------------------------------------
    template <typename T>
    static void GatherV(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                        int root, MPI_Comm comm = MPI_COMM_WORLD);
};

}  // namespace pwr

#include "mpi_utilities.tcc"

#endif  // PWR_MPI_UTILITIES_H
