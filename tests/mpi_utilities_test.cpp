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

    // ------------------------------------------------------------------------
    // MPI all gather
    // ------------------------------------------------------------------------
    {
        INFO("MPI all gather");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Each rank sends itself
        const int send = rank;

        // Buffer to receive
        std::vector<int> recv;

        // MPI all gather
        pwr::MPIUtilities::AllGather(send, recv);

        // -- Check if size matches gather vector -----------------------------
        if (recv.size() != static_cast<std::size_t>(size)) {
            num_failed_tests_local++;
        }

        // Loop each rank
        for (int i = 0; i < size; ++i) {
            // -- Check if recv value matches solution ------------------------
            if (i != recv[i]) {
                num_failed_tests_local++;
            }
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

    }  // MPI all gather

    // ------------------------------------------------------------------------
    // MPI all gather (equal-sized vector)
    // ------------------------------------------------------------------------
    {
        INFO("MPI all gather (equal-sized vector)");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Each rank sends {itself, itself + offset}
        const int base = rank;
        std::vector<int> send = {base, base + 100};

        // Buffer to receive
        std::vector<int> recv;

        // MPI all gather
        pwr::MPIUtilities::AllGather(send, recv);

        // -- Check if size matches gather vector -----------------------------
        if (recv.size() != static_cast<std::size_t>(2 * size)) {
            num_failed_tests_local++;
        }

        // Loop each rank
        for (int i = 0; i < size; ++i) {
            // -- Check if recv value matches solution ------------------------
            if (i != recv[2 * i + 0]) {
                num_failed_tests_local++;
            }
            if ((i + 100) != recv[2 * i + 1]) {
                num_failed_tests_local++;
            }
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

    }  // MPI all gather (equal-sized vector)

    // ------------------------------------------------------------------------
    // MPI all gather variable
    // ------------------------------------------------------------------------
    {
        INFO("MPI all gather variable");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Each rank sends variable length vector
        std::vector<int> send(rank + 1);
        for (int i = 0; i < (rank + 1); ++i) {
            send[i] = 100 * rank + i;
        }

        // Expected solution
        std::vector<int> expected;
        for (int r = 0; r < size; ++r) {
            for (int i = 0; i < (r + 1); ++i) {
                expected.push_back(100 * r + i);
            }
        }

        // Buffer to receive
        std::vector<int> recv;

        // MPI all gather variable
        pwr::MPIUtilities::AllGatherV(send, recv);

        // -- Check if expected size matches gather vector --------------------
        if (recv.size() != expected.size()) {
            num_failed_tests_local++;
        }

        // Loop recv vector
        for (int i = 0; i < recv.size(); ++i) {
            // -- Check if recv value matches expected ------------------------
            if (recv[i] != expected[i]) {
                num_failed_tests_local++;
            }
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

    }  // MPI all gather variable

    // ------------------------------------------------------------------------
    // MPI gather
    // ------------------------------------------------------------------------
    {
        INFO("MPI gather");

        // Each rank sends itself
        const int send = rank;

        // Buffer to receive
        std::vector<int> recv;

        // MPI gather
        const int root = 0;
        pwr::MPIUtilities::Gather(send, recv, root);

        // Only check on root
        if (rank == root) {
            // Keep track of failed tests on this rank
            int num_failed_tests_local = 0;

            // -- Check if size matches gather vector -------------------------
            if (recv.size() != static_cast<std::size_t>(size)) {
                num_failed_tests_local++;
            }

            // Loop each rank
            for (int i = 0; i < size; ++i) {
                // -- Check if recv value matches solution --------------------
                if (i != recv[i]) {
                    num_failed_tests_local++;
                }
            }

            // Only check on root (rank 0)
            REQUIRE(num_failed_tests_local == 0);
        }

    }  // MPI gather

    // ------------------------------------------------------------------------
    // MPI gather (equal-sized vector)
    // ------------------------------------------------------------------------
    {
        INFO("MPI gather (equal-sized vector)");

        // Each rank sends {itself, itself + offset}
        const int base = rank;
        std::vector<int> send = {base, base + 100};

        // Buffer to receive
        std::vector<int> recv;

        // MPI gather
        const int root = 0;
        pwr::MPIUtilities::Gather(send, recv, root);

        // Only check on root
        if (rank == root) {
            // Keep track of failed tests on this rank
            int num_failed_tests_local = 0;

            // -- Check if size matches gather vector -------------------------
            if (recv.size() != static_cast<std::size_t>(2 * size)) {
                num_failed_tests_local++;
            }

            // Loop each rank
            for (int i = 0; i < size; ++i) {
                // -- Check if recv value matches solution --------------------
                if (i != recv[2 * i + 0]) {
                    num_failed_tests_local++;
                }
                if ((i + 100) != recv[2 * i + 1]) {
                    num_failed_tests_local++;
                }
            }

            // Only check on root (rank 0)
            REQUIRE(num_failed_tests_local == 0);
        }

    }  // MPI gather (equal-sized vector)

    // ------------------------------------------------------------------------
    // MPI gather variable
    // ------------------------------------------------------------------------
    {
        INFO("MPI gather variable");

        // Each rank sends variable length vector
        std::vector<int> send(rank + 1);
        for (int i = 0; i < (rank + 1); ++i) {
            send[i] = 100 * rank + i;
        }

        // Buffer to receive
        std::vector<int> recv;

        // MPI gather variable
        const int root = 0;
        pwr::MPIUtilities::GatherV(send, recv, root);

        // Only check on root
        if (rank == root) {
            // Keep track of failed tests on this rank
            int num_failed_tests_local = 0;

            // Expected solution
            std::vector<int> expected;
            for (int r = 0; r < size; ++r) {
                for (int i = 0; i < (r + 1); ++i) {
                    expected.push_back(100 * r + i);
                }
            }

            // -- Check if expected size matches gather vector ----------------
            if (recv.size() != expected.size()) {
                num_failed_tests_local++;
            }

            // Loop recv vector
            for (int i = 0; i < recv.size(); ++i) {
                // -- Check if recv value matches expected --------------------
                if (recv[i] != expected[i]) {
                    num_failed_tests_local++;
                }
            }

            // Only check on root (rank 0)
            REQUIRE(num_failed_tests_local == 0);
        }

    }  // MPI gather variable

}  // TEST_CASE("MPI Utilities")
