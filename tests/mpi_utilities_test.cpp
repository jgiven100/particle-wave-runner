#include "mpi_utilities.h"

#include <catch2/catch_all.hpp>

TEST_CASE("MPI Utilities", "[mpi]") {
    // Grab rank and size
    const int rank = pwr::MPIUtilities::Rank();
    const int size = pwr::MPIUtilities::Size();

    // ------------------------------------------------------------------------
    // MPI datatype
    // ------------------------------------------------------------------------
    {
        INFO("MPI datatype");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // -- Check type: MPI_CHAR --------------------------------------------
        if (pwr::MPIType<char>::type != MPI_CHAR) {
            num_failed_tests_local++;
        }

        // -- Check type: MPI_INT ---------------------------------------------
        if (pwr::MPIType<int>::type != MPI_INT) {
            num_failed_tests_local++;
        }

        // -- Check type: MPI_UNSIGNED ----------------------------------------
        if (pwr::MPIType<unsigned>::type != MPI_UNSIGNED) {
            num_failed_tests_local++;
        }

        // -- Check type: MPI_UNSIGNED_LONG -----------------------------------
        if (pwr::MPIType<unsigned long>::type != MPI_UNSIGNED_LONG) {
            num_failed_tests_local++;
        }

        // -- Check type: MPI_UNSIGNED_LONG_LONG ------------------------------
        if (pwr::MPIType<unsigned long long>::type != MPI_UNSIGNED_LONG_LONG) {
            num_failed_tests_local++;
        }

        // -- Check type: MPI_FLOAT -------------------------------------------
        if (pwr::MPIType<float>::type != MPI_FLOAT) {
            num_failed_tests_local++;
        }

        // -- Check type: MPI_DOUBLE ------------------------------------------
        if (pwr::MPIType<double>::type != MPI_DOUBLE) {
            num_failed_tests_local++;
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root (rank 0)
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI datatype

    // ------------------------------------------------------------------------
    // MPI all reduce maximum
    // ------------------------------------------------------------------------
    {
        INFO("MPI all reduce maximum");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Set global values
        const std::array<int, 4> values = {40, 30, 20, 10};

        // Grab send value for this rank
        const int send = values[rank];

        // MPI all reduce maximum
        const int recv = pwr::MPIUtilities::AllReduceMax(send);

        // -- Check if recv matches solution ----------------------------------
        if (recv != 40) {
            num_failed_tests_local++;
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root (rank 0)
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all reduce maximum

    // ------------------------------------------------------------------------
    // MPI all reduce minimum
    // ------------------------------------------------------------------------
    {
        INFO("MPI all reduce minimum");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Set global values
        const std::array<int, 4> values = {10, 20, 30, 40};

        // Grab send value for this rank
        const int send = values[rank];

        // MPI all reduce minimum
        const int recv = pwr::MPIUtilities::AllReduceMin(send);

        // -- Check if recv matches solution ----------------------------------
        if (recv != 10) {
            num_failed_tests_local++;
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root (rank 0)
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all reduce minimum

    // ------------------------------------------------------------------------
    // MPI all reduce sum
    // ------------------------------------------------------------------------
    {
        INFO("MPI all reduce sum");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Set global solution
        const std::array<int, 4> solutions = {10, 20, 30, 40};

        // Set send value for all ranks
        const int send = 10;

        // MPI all reduce sum
        const int recv = pwr::MPIUtilities::AllReduceSum(send);

        // -- Check if recv matches solution ----------------------------------
        if (recv != solutions[size - 1]) {
            num_failed_tests_local++;
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root (rank 0)
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all reduce sum

    // ------------------------------------------------------------------------
    // MPI all reduce maximum in place
    // ------------------------------------------------------------------------
    {
        INFO("MPI all reduce maximum in place");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Set global values
        const std::array<int, 4> values = {40, 30, 20, 10};

        // Grab send/recv for this rank
        int value = values[rank];

        // MPI all reduce maximum in place
        pwr::MPIUtilities::AllReduceMaxInPlace(value);

        // -- Check if send/recv matches solution -----------------------------
        if (value != 40) {
            num_failed_tests_local++;
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root (rank 0)
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all reduce maximum in place

    // ------------------------------------------------------------------------
    // MPI all reduce minimum in place
    // ------------------------------------------------------------------------
    {
        INFO("MPI all reduce minimum in place");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Set global values
        const std::array<int, 4> values = {10, 20, 30, 40};

        // Grab send/recv value for this rank
        int value = values[rank];

        // MPI all reduce minimum in place
        pwr::MPIUtilities::AllReduceMinInPlace(value);

        // -- Check if send/recv matches solution -----------------------------
        if (value != 10) {
            num_failed_tests_local++;
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root (rank 0)
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all reduce minimum in place

    // ------------------------------------------------------------------------
    // MPI all reduce sum in place
    // ------------------------------------------------------------------------
    {
        INFO("MPI all reduce sum in place");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Set global solution
        const std::array<int, 4> values = {10, 20, 30, 40};

        // Set send/recv for all ranks
        int value = 10;

        // MPI all reduce sum in place
        pwr::MPIUtilities::AllReduceSumInPlace(value);

        // -- Check if send/recv matches solution -----------------------------
        if (value != values[size - 1]) {
            num_failed_tests_local++;
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root (rank 0)
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all reduce sum in place

    // ------------------------------------------------------------------------
    // MPI reduce maximum
    // ------------------------------------------------------------------------
    {
        INFO("MPI reduce maximum");

        // Set global values
        const std::array<int, 4> values = {40, 30, 20, 10};

        // Grab send value for this rank
        const int send = values[rank];

        // Set root (rank 0)
        const int root = 0;

        // MPI reduce maximum
        int recv = 0;
        pwr::MPIUtilities::ReduceMax(send, recv, root);

        // Only check on root (rank 0)
        if (rank == root) {
            REQUIRE(recv == 40);
        }

    }  // MPI reduce maximum

    // ------------------------------------------------------------------------
    // MPI reduce minimum
    // ------------------------------------------------------------------------
    {
        INFO("MPI reduce minimum");

        // Set global values
        const std::array<int, 4> values = {10, 20, 30, 40};

        // Grab send value for this rank
        const int send = values[rank];

        // Set root (rank 0)
        const int root = 0;

        // MPI reduce minimum
        int recv = 0;
        pwr::MPIUtilities::ReduceMin(send, recv, root);

        // Only check on root (rank 0)
        if (rank == root) {
            REQUIRE(recv == 10);
        }

    }  // MPI reduce minimum

    // ------------------------------------------------------------------------
    // MPI reduce sum
    // ------------------------------------------------------------------------
    {
        INFO("MPI reduce sum");

        // Set global solution
        const std::array<int, 4> values = {10, 20, 30, 40};

        // Set send for all ranks
        const int send = 10;

        // Set root (rank 0)
        const int root = 0;

        // MPI reduce sum
        int recv = 0;
        pwr::MPIUtilities::ReduceSum(send, recv, root);

        // Only check on root (rank 0)
        if (rank == root) {
            REQUIRE(recv == values[size - 1]);
        }

    }  // MPI reduce sum

}  // TEST_CASE("MPI Utilities")
