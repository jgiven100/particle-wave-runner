#include <mpi.h>

#include <cassert>
#include <numeric>
#include <type_traits>
#include <vector>

#include "mpi_utilities.h"

namespace pwr {

// ----------------------------------------------------------------------------
// Get MPI datatype
// ----------------------------------------------------------------------------
template <>
struct MPIType<char> {
    static inline MPI_Datatype type = MPI_CHAR;
};

template <>
struct MPIType<int> {
    static inline MPI_Datatype type = MPI_INT;
};

template <>
struct MPIType<unsigned> {
    static inline MPI_Datatype type = MPI_UNSIGNED;
};

template <>
struct MPIType<unsigned long> {
    static inline MPI_Datatype type = MPI_UNSIGNED_LONG;
};

template <>
struct MPIType<unsigned long long> {
    static inline MPI_Datatype type = MPI_UNSIGNED_LONG_LONG;
};

template <>
struct MPIType<float> {
    static inline MPI_Datatype type = MPI_FLOAT;
};

template <>
struct MPIType<double> {
    static inline MPI_Datatype type = MPI_DOUBLE;
};

// ----------------------------------------------------------------------------
// MPI all reduce
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::AllReduce(const T& sendbuf, T& recvbuf, MPI_Op op,
                             MPI_Comm comm) {
    // Sanity check: compatible MPI_Allreduce type
    static_assert(std::is_trivially_copyable_v<T>,
                  "MPI_Allreduce requires trivially copyable type");

    int err = MPI_Allreduce(    // All reduce
        &sendbuf,               // Send buffer
        &recvbuf,               // Receive buffer
        1,                      // Count
        pwr::MPIType<T>::type,  // MPI datatype
        op,                     // MPI operation
        comm);                  // MPI communicator

    // Check for successful MPI_Allreduce
    if (err != MPI_SUCCESS) {
        MPI_Abort(comm, err);
    }

}  // MPIUtilities::AllReduce

// ----------------------------------------------------------------------------
// MPI all reduce maximum
// ----------------------------------------------------------------------------
template <typename T>
T MPIUtilities::AllReduceMax(const T& sendbuf, MPI_Comm comm) {
    // Receive buffer
    T recvbuf;

    // Call `AllReduce` with MPI_MAX
    AllReduce(sendbuf, recvbuf, MPI_MAX, comm);

    return recvbuf;

}  // MPIUtilities::AllReduceMax

// ----------------------------------------------------------------------------
// MPI all reduce minimum
// ----------------------------------------------------------------------------
template <typename T>
T MPIUtilities::AllReduceMin(const T& sendbuf, MPI_Comm comm) {
    // Receive buffer
    T recvbuf;

    // Call `AllReduce` with MPI_MIN
    AllReduce(sendbuf, recvbuf, MPI_MIN, comm);

    return recvbuf;

}  // MPIUtilities::AllReduceMin

// ----------------------------------------------------------------------------
// MPI all reduce sum
// ----------------------------------------------------------------------------
template <typename T>
T MPIUtilities::AllReduceSum(const T& sendbuf, MPI_Comm comm) {
    // Receive buffer
    T recvbuf;

    // Call `AllReduce` with MPI_SUM
    AllReduce(sendbuf, recvbuf, MPI_SUM, comm);

    return recvbuf;

}  // MPIUtilities::AllReduceSum

// ----------------------------------------------------------------------------
// MPI all reduce in place
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::AllReduceInPlace(T& send_recv_buf, MPI_Op op,
                                    MPI_Comm comm) {
    // Sanity check: compatible MPI_Allreduce type
    static_assert(std::is_trivially_copyable_v<T>,
                  "MPI_Allreduce requires trivially copyable type");

    int err = MPI_Allreduce(    // All reduce
        MPI_IN_PLACE,           // MPI_IN_PLACE value
        &send_recv_buf,         // Send and receive buffer
        1,                      // Count
        pwr::MPIType<T>::type,  // MPI datatype
        op,                     // MPI operation
        comm);                  // MPI communicator

    // Check for successful MPI_Allreduce
    if (err != MPI_SUCCESS) {
        MPI_Abort(comm, err);
    }

}  // MPIUtilities::AllReduceInPlace

// ----------------------------------------------------------------------------
// MPI all reduce maximum in place
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::AllReduceMaxInPlace(T& send_recv_buf, MPI_Comm comm) {
    // Call `AllReduceInPlace` with MPI_MAX
    AllReduceInPlace(send_recv_buf, MPI_MAX, comm);

}  // MPIUtilities::AllReduceMaxInPlace

// ----------------------------------------------------------------------------
// MPI all reduce minimum in place
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::AllReduceMinInPlace(T& send_recv_buf, MPI_Comm comm) {
    // Call `AllReduceInPlace` with MPI_MIN
    AllReduceInPlace(send_recv_buf, MPI_MIN, comm);

}  // MPIUtilities::AllReduceMinInPlace

// ----------------------------------------------------------------------------
// MPI all reduce sum in place
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::AllReduceSumInPlace(T& send_recv_buf, MPI_Comm comm) {
    // Call `AllReduceInPlace` with MPI_SUM
    AllReduceInPlace(send_recv_buf, MPI_SUM, comm);

}  // MPIUtilities::AllReduceSumInPlace

// ----------------------------------------------------------------------------
// MPI reduce
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::Reduce(const T& sendbuf, T& recvbuf, MPI_Op op, int root,
                          MPI_Comm comm) {
    // Sanity check: compatible MPI_Reduce type
    static_assert(std::is_trivially_copyable_v<T>,
                  "MPI_Reduce requires trivially copyable type");

    int err = MPI_Reduce(       // Reduce
        &sendbuf,               // Send buffer
        &recvbuf,               // Receive buffer
        1,                      // Count
        pwr::MPIType<T>::type,  // MPI datatype
        op,                     // MPI operation
        root,                   // Root
        comm);                  // MPI communicator

    // Check for successful MPI_Reduce
    if (err != MPI_SUCCESS) {
        MPI_Abort(comm, err);
    }

}  // MPIUtilities::Reduce

// ----------------------------------------------------------------------------
// MPI reduce maximum
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::ReduceMax(const T& sendbuf, T& recvbuf, int root,
                             MPI_Comm comm) {
    // Call `Reduce` with MPI_MAX
    Reduce(sendbuf, recvbuf, MPI_MAX, root, comm);

}  // MPIUtilities::ReduceMax

// ----------------------------------------------------------------------------
// MPI reduce minimum
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::ReduceMin(const T& sendbuf, T& recvbuf, int root,
                             MPI_Comm comm) {
    // Call `Reduce` with MPI_MIN
    Reduce(sendbuf, recvbuf, MPI_MIN, root, comm);

}  // MPIUtilities::ReduceMin

// ----------------------------------------------------------------------------
// MPI reduce sum
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::ReduceSum(const T& sendbuf, T& recvbuf, int root,
                             MPI_Comm comm) {
    // Call `Reduce` with MPI_SUM
    Reduce(sendbuf, recvbuf, MPI_SUM, root, comm);

}  // MPIUtilities::ReduceSum

// ----------------------------------------------------------------------------
// MPI all gather
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::AllGather(const T& sendbuf, std::vector<T>& recvbuf,
                             MPI_Comm comm) {
    // Sanity check: compatible MPI_Allgather type
    static_assert(std::is_trivially_copyable_v<T>,
                  "MPI_Allgather requires trivially copyable type");

    // Get size
    int size = Size(comm);

    // Resize recvbuf
    recvbuf.resize(size);

    int err = MPI_Allgather(    // All gather
        &sendbuf,               // Send buffer
        1,                      // Send count
        pwr::MPIType<T>::type,  // Send MPI type
        recvbuf.data(),         // Receive buffer
        1,                      // Receive count per rank
        pwr::MPIType<T>::type,  // Receive MPI type
        comm);                  // MPI communicator

    // Check for successful MPI_Allgather
    if (err != MPI_SUCCESS) {
        MPI_Abort(comm, err);
    }

}  // MPIUtilities::AllGather

// ----------------------------------------------------------------------------
// MPI all gather (equal-sized vector)
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::AllGather(const std::vector<T>& sendbuf,
                             std::vector<T>& recvbuf, MPI_Comm comm) {
    // Sanity check: compatible MPI_Allgather type
    static_assert(std::is_trivially_copyable_v<T>,
                  "MPI_Allgather requires trivially copyable type");

    // Get size
    int size = Size(comm);
    int local_n = static_cast<int>(sendbuf.size());

    // Collect min and max vector sizes
    int min_n = AllReduceMin(local_n, comm);
    int max_n = AllReduceMax(local_n, comm);

    // Check for equal-sized vectors across all ranks
    if (min_n != max_n) {
        MPI_Abort(comm, MPI_ERR_INTERN);
    }

    // Resize recvbuf (each rank contributes equal-sized vector)
    recvbuf.resize(size * local_n);

    int err = MPI_Allgather(    // All gather
        sendbuf.data(),         // Send buffer
        local_n,                // Send count
        pwr::MPIType<T>::type,  // Send MPI type
        recvbuf.data(),         // Receive buffer
        local_n,                // Receive count per rank
        pwr::MPIType<T>::type,  // Receive MPI type
        comm);                  // MPI communicator

    // Check for successful MPI_Allgather
    if (err != MPI_SUCCESS) {
        MPI_Abort(comm, err);
    }

}  // MPIUtilities::AllGather

// ----------------------------------------------------------------------------
// MPI all gather variable
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::AllGatherV(const std::vector<T>& sendbuf,
                              std::vector<T>& recvbuf, MPI_Comm comm) {
    // Sanity check: compatible MPI_Allgatherv type
    static_assert(std::is_trivially_copyable_v<T>,
                  "MPI_Allgatherv requires trivially copyable type");

    // Get local size
    int local_n = static_cast<int>(sendbuf.size());

    // Gather recvcounts
    std::vector<int> recvcounts;
    AllGather(local_n, recvcounts, comm);

    // Compute displacements
    std::vector<int> displs(recvcounts.size());
    int offset = 0;
    for (std::size_t i = 0; i < recvcounts.size(); ++i) {
        displs[i] = offset;
        offset += recvcounts[i];
    }

    // Resize recvbuf
    recvbuf.resize(offset);

    int err = MPI_Allgatherv(   // All gather variable
        sendbuf.data(),         // Send buffer
        local_n,                // Send count
        pwr::MPIType<T>::type,  // Send MPI type
        recvbuf.data(),         // Receive buffer
        recvcounts.data(),      // Receive counts
        displs.data(),          // Displacements
        pwr::MPIType<T>::type,  // Receive MPI type
        comm);                  // MPI communicator

    // Check for successful MPI_Allgatherv
    if (err != MPI_SUCCESS) {
        MPI_Abort(comm, err);
    }

}  // MPIUtilities::AllGatherV

// ----------------------------------------------------------------------------
// MPI gather
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::Gather(const T& sendbuf, std::vector<T>& recvbuf, int root,
                          MPI_Comm comm) {
    // Sanity check: compatible MPI_Gather type
    static_assert(std::is_trivially_copyable_v<T>,
                  "MPI_Gather requires trivially copyable type");

    // Get rank and size
    int rank = Rank(comm);
    int size = Size(comm);

    // Resize recvbuf
    if (rank == root) {
        recvbuf.resize(size);
    }

    int err = MPI_Gather(                           // Gather
        &sendbuf,                                   // Send buffer
        1,                                          // Send count
        pwr::MPIType<T>::type,                      // Send MPI type
        (rank == root ? recvbuf.data() : nullptr),  // Receive buffer
        1,                                          // Receive count per rank
        pwr::MPIType<T>::type,                      // Receive MPI type
        root,                                       // Root
        comm);                                      // MPI communicator

    // Check for successful MPI_Gather
    if (err != MPI_SUCCESS) {
        MPI_Abort(comm, err);
    }

}  // MPIUtilities::Gather

// ----------------------------------------------------------------------------
// MPI gather (equal-sized vector)
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::Gather(const std::vector<T>& sendbuf,
                          std::vector<T>& recvbuf, int root, MPI_Comm comm) {
    // Sanity check: compatible MPI_Gather type
    static_assert(std::is_trivially_copyable_v<T>,
                  "MPI_Gather requires trivially copyable type");

    // Get rank and size
    int rank = Rank(comm);
    int size = Size(comm);
    int local_n = static_cast<int>(sendbuf.size());

    // Collect min and max vector sizes
    int min_n = AllReduceMin(local_n, comm);
    int max_n = AllReduceMax(local_n, comm);

    // Check for equal-sized vectors across all ranks
    if (min_n != max_n) {
        MPI_Abort(comm, MPI_ERR_INTERN);
    }

    // Resize recvbuf
    if (rank == root) {
        recvbuf.resize(size * local_n);
    }

    int err = MPI_Gather(                           // Gather
        sendbuf.data(),                             // Send buffer
        local_n,                                    // Send count
        pwr::MPIType<T>::type,                      // Send MPI type
        (rank == root ? recvbuf.data() : nullptr),  // Receive buffer
        local_n,                                    // Receive count per rank
        pwr::MPIType<T>::type,                      // Receive MPI type
        root,                                       // Root
        comm);                                      // MPI communicator

    // Check for successful MPI_Gather
    if (err != MPI_SUCCESS) {
        MPI_Abort(comm, err);
    }

}  // MPIUtilities::Gather

// ----------------------------------------------------------------------------
// MPI gather variable
// ----------------------------------------------------------------------------
template <typename T>
void MPIUtilities::GatherV(const std::vector<T>& sendbuf,
                           std::vector<T>& recvbuf, int root, MPI_Comm comm) {
    // Sanity check: compatible MPI_Gatherv type
    static_assert(std::is_trivially_copyable_v<T>,
                  "MPI_Gatherv requires trivially copyable type");

    // Get rank and local size
    int rank = Rank(comm);
    int local_n = static_cast<int>(sendbuf.size());

    // Gather recvcounts
    std::vector<int> recvcounts;
    Gather(local_n, recvcounts, root, comm);

    // Displacements
    std::vector<int> displs;

    if (rank == root) {
        // Compute displacements
        displs.resize(recvcounts.size());
        int offset = 0;
        for (std::size_t i = 0; i < recvcounts.size(); ++i) {
            displs[i] = offset;
            offset += recvcounts[i];
        }

        // Resize recvbuf
        recvbuf.resize(offset);
    }

    int err = MPI_Gatherv(                             // Gather variable
        sendbuf.data(),                                // Send buffer
        local_n,                                       // Send count
        pwr::MPIType<T>::type,                         // Send MPI type
        (rank == root ? recvbuf.data() : nullptr),     // Receive buffer
        (rank == root ? recvcounts.data() : nullptr),  // Receive counts
        (rank == root ? displs.data() : nullptr),      // Displacements
        pwr::MPIType<T>::type,                         // Receive MPI type
        root,                                          // Root
        comm);                                         // MPI communicator

    // Check for successful MPI_Gatherv
    if (err != MPI_SUCCESS) {
        MPI_Abort(comm, err);
    }

}  // MPIUtilities::GatherV

}  // namespace pwr
