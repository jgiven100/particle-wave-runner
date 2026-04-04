#include <mpi.h>

#include <cassert>
#include <type_traits>

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
    // Sanity check: compatable MPI_Allreduce type
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
    // Sanity check: compatable MPI_Allreduce type
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
    // Sanity check: compatable MPI_Reduce type
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

}  // namespace pwr
