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
    // MPI all reduce (scalar) maximum
    // ------------------------------------------------------------------------
    {
        INFO("MPI all reduce (scalar) maximum");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Each rank sends itself
        const int send = rank;

        // Buffer to receive
        int recv;

        // MPI all reduce maximum
        pwr::MPIUtilities::AllReduceMax(send, recv);

        // -- Check if recv matches solution ----------------------------------
        if (recv != (size - 1)) {
            num_failed_tests_local++;
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all reduce scalar maximum

    // ------------------------------------------------------------------------
    // MPI all reduce (vector) maximum
    // ------------------------------------------------------------------------
    {
        INFO("MPI all reduce (vector) maximum");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Each rank sends {itself, itself + offset, ...}
        const int base = rank;
        const std::vector<int> send = {base, base + 1, base + 2};

        // Buffer to receive
        std::vector<int> recv;

        // MPI all reduce maximum
        pwr::MPIUtilities::AllReduceMax(send, recv);

        // Loop each element
        for (std::size_t i = 0; i < recv.size(); ++i) {
            // -- Check if recv matches solution ------------------------------
            if (recv[i] != (size - 1 + i)) {
                num_failed_tests_local++;
            }
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all reduce (vector) maximum

    // ------------------------------------------------------------------------
    // MPI all reduce (scalar) minimum
    // ------------------------------------------------------------------------
    {
        INFO("MPI all reduce (scalar) minimum");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Each rank sends itself
        const int send = rank;

        // Buffer to receive
        int recv;

        // MPI all reduce minimum
        pwr::MPIUtilities::AllReduceMin(send, recv);

        // -- Check if recv matches solution ----------------------------------
        if (recv != 0) {
            num_failed_tests_local++;
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all reduce (scalar) minimum

    // ------------------------------------------------------------------------
    // MPI all reduce (vector) minimum
    // ------------------------------------------------------------------------
    {
        INFO("MPI all reduce (vector) minimum");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Each rank sends {itself, itself + offset, ...}
        const int base = rank;
        const std::vector<int> send = {base, base + 1, base + 2};

        // Buffer to receive
        std::vector<int> recv;

        // MPI all reduce minimum
        pwr::MPIUtilities::AllReduceMin(send, recv);

        // Loop each element
        for (std::size_t i = 0; i < recv.size(); ++i) {
            // -- Check if recv matches solution ------------------------------
            if (recv[i] != i) {
                num_failed_tests_local++;
            }
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all reduce (vector) minimum

    // ------------------------------------------------------------------------
    // MPI all reduce (scalar) sum
    // ------------------------------------------------------------------------
    {
        INFO("MPI all reduce (scalar) sum");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Each rank sends itself
        const int send = rank;

        // Buffer to receive
        int recv;

        // MPI all reduce sum
        pwr::MPIUtilities::AllReduceSum(send, recv);

        // Expected solution
        int expected = 0;
        for (int r = 0; r < size; ++r) {
            expected += r;
        }

        // -- Check if recv matches solution ----------------------------------
        if (recv != expected) {
            num_failed_tests_local++;
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all reduce (scalar) sum

    // ------------------------------------------------------------------------
    // MPI all reduce (vector) sum
    // ------------------------------------------------------------------------
    {
        INFO("MPI all reduce (vector) sum");

        // Keep track of failed tests on this rank
        int num_failed_tests_local = 0;

        // Each rank sends {itself, itself + offset, ...}
        const int base = rank;
        const std::vector<int> send = {base, base + 1, base + 2};

        // Buffer to receive
        std::vector<int> recv;

        // MPI all reduce sum
        pwr::MPIUtilities::AllReduceSum(send, recv);

        // Expected solution
        std::vector<int> expected = {0, 0, 0};
        for (int r = 0; r < size; ++r) {
            expected[0] += r + 0;
            expected[1] += r + 1;
            expected[2] += r + 2;
        }

        // Loop each element
        for (std::size_t i = 0; i < recv.size(); ++i) {
            // -- Check if recv matches solution ------------------------------
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

        // Only check on root
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all reduce (vector) sum

    // ------------------------------------------------------------------------
    // MPI reduce (scalar) maximum
    // ------------------------------------------------------------------------
    {
        INFO("MPI reduce (scalar) maximum ");

        // Each rank sends itself
        const int send = rank;

        // Buffer to receive
        int recv;

        // Set root (rank 0)
        const int root = 0;

        // MPI reduce maximum
        pwr::MPIUtilities::ReduceMax(send, recv, root);

        // Only check on root
        if (rank == root) {
            REQUIRE(recv == (size - 1));
        }

    }  // MPI reduce (scalar) maximum

    // ------------------------------------------------------------------------
    // MPI reduce (vector) maximum
    // ------------------------------------------------------------------------
    {
        INFO("MPI reduce (vector) maximum");

        // Each rank sends {itself, itself + offset, ...}
        const int base = rank;
        const std::vector<int> send = {base, base + 1, base + 2};

        // Buffer to receive
        std::vector<int> recv;

        // Set root (rank 0)
        const int root = 0;

        // MPI reduce maximum
        pwr::MPIUtilities::ReduceMax(send, recv, root);

        // Only check on root
        if (rank == root) {
            // Keep track of failed tests on this rank
            int num_failed_tests_local = 0;

            // Loop recv vector
            for (int i = 0; i < recv.size(); ++i) {
                // -- Check if recv value matches expected --------------------
                if (recv[i] != (size - 1 + i)) {
                    num_failed_tests_local++;
                }
            }

            // Only check on root
            REQUIRE(num_failed_tests_local == 0);
        }

    }  // MPI reduce (vector) maximum

    // ------------------------------------------------------------------------
    // MPI reduce (scalar) minimum
    // ------------------------------------------------------------------------
    {
        INFO("MPI reduce (scalar) minimum");

        // Each rank sends itself
        const int send = rank;

        // Buffer to receive
        int recv;

        // Set root (rank 0)
        const int root = 0;

        // MPI reduce minimum
        pwr::MPIUtilities::ReduceMin(send, recv, root);

        // Only check on root
        if (rank == root) {
            REQUIRE(recv == 0);
        }

    }  // MPI reduce (scalar) minimum

    // ------------------------------------------------------------------------
    // MPI reduce (vector) minimum
    // ------------------------------------------------------------------------
    {
        INFO("MPI reduce (vector) minimum");

        // Each rank sends {itself, itself + offset, ...}
        const int base = rank;
        const std::vector<int> send = {base, base + 1, base + 2};

        // Buffer to receive
        std::vector<int> recv;

        // Set root (rank 0)
        const int root = 0;

        // MPI reduce minimum
        pwr::MPIUtilities::ReduceMin(send, recv, root);

        // Only check on root
        if (rank == root) {
            // Keep track of failed tests on this rank
            int num_failed_tests_local = 0;

            // Loop recv vector
            for (int i = 0; i < recv.size(); ++i) {
                // -- Check if recv value matches expected --------------------
                if (recv[i] != i) {
                    num_failed_tests_local++;
                }
            }

            // Only check on root
            REQUIRE(num_failed_tests_local == 0);
        }

    }  // MPI reduce (vector) minimum

    // ------------------------------------------------------------------------
    // MPI reduce (scalar) sum
    // ------------------------------------------------------------------------
    {
        INFO("MPI reduce (scalar) sum");

        // Each rank sends itself
        const int send = rank;

        // Buffer to receive
        int recv;

        // Set root (rank 0)
        const int root = 0;

        // MPI reduce sum
        pwr::MPIUtilities::ReduceSum(send, recv, root);

        // Only check on root
        if (rank == root) {
            // Expected solution
            int expected = 0;
            for (int r = 0; r < size; ++r) {
                expected += r;
            }

            // Only check on root
            REQUIRE(recv == expected);
        }

    }  // MPI reduce (scalar) sum

    // ------------------------------------------------------------------------
    // MPI reduce (vector) sum
    // ------------------------------------------------------------------------
    {
        INFO("MPI reduce (vector) sum");

        // Each rank sends {itself, itself + offset, ...}
        const int base = rank;
        const std::vector<int> send = {base, base + 1, base + 2};

        // Buffer to receive
        std::vector<int> recv;

        // Set root (rank 0)
        const int root = 0;

        // MPI reduce sum
        pwr::MPIUtilities::ReduceSum(send, recv, root);

        // Only check on root
        if (rank == root) {
            // Keep track of failed tests on this rank
            int num_failed_tests_local = 0;

            // Expected solution
            std::vector<int> expected = {0, 0, 0};
            for (int r = 0; r < size; ++r) {
                expected[0] += r + 0;
                expected[1] += r + 1;
                expected[2] += r + 2;
            }

            // Loop recv vector
            for (int i = 0; i < recv.size(); ++i) {
                // -- Check if recv value matches expected --------------------
                if (recv[i] != expected[i]) {
                    num_failed_tests_local++;
                }
            }

            // Only check on root
            REQUIRE(num_failed_tests_local == 0);
        }

    }  // MPI reduce (vector) sum

    // ------------------------------------------------------------------------
    // MPI all gather (scalar)
    // ------------------------------------------------------------------------
    {
        INFO("MPI all gather (scalar)");

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
            if (recv[i] != i) {
                num_failed_tests_local++;
            }
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all gather (scalar)

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
            if (recv[2 * i + 0] != i) {
                num_failed_tests_local++;
            }
            if (recv[2 * i + 1] != (i + 100)) {
                num_failed_tests_local++;
            }
        }

        // Set root (rank 0)
        const int root = 0;

        // Collect results
        int num_failed_tests_global;
        pwr::MPIUtilities::ReduceSum(num_failed_tests_local,
                                     num_failed_tests_global, root);

        // Only check on root
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all gather (equal-sized vector)

    // ------------------------------------------------------------------------
    // MPI all gather (variable-sized vector)
    // ------------------------------------------------------------------------
    {
        INFO("MPI all gather (variable-sized vector)");

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

        // Only check on root
        if (rank == root) {
            REQUIRE(num_failed_tests_global == 0);
        }

    }  // MPI all gather (variable-sized vector)

    // ------------------------------------------------------------------------
    // MPI gather (scalar)
    // ------------------------------------------------------------------------
    {
        INFO("MPI gather (scalar)");

        // Each rank sends itself
        const int send = rank;

        // Buffer to receive
        std::vector<int> recv;

        // Set root (rank 0)
        const int root = 0;

        // MPI gather
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
                if (recv[i] != i) {
                    num_failed_tests_local++;
                }
            }

            // Only check on root
            REQUIRE(num_failed_tests_local == 0);
        }

    }  // MPI gather (scalar)

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

        // Set root (rank 0)
        const int root = 0;

        // MPI gather
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
                if (recv[2 * i + 0] != i) {
                    num_failed_tests_local++;
                }
                if (recv[2 * i + 1] != (i + 100)) {
                    num_failed_tests_local++;
                }
            }

            // Only check on root
            REQUIRE(num_failed_tests_local == 0);
        }

    }  // MPI gather (equal-sized vector)

    // ------------------------------------------------------------------------
    // MPI gather (variable-sized vector)
    // ------------------------------------------------------------------------
    {
        INFO("MPI gather (variable-sized vector)");

        // Each rank sends variable length vector
        std::vector<int> send(rank + 1);
        for (int i = 0; i < (rank + 1); ++i) {
            send[i] = 100 * rank + i;
        }

        // Buffer to receive
        std::vector<int> recv;

        // Set root (rank 0)
        const int root = 0;

        // MPI gather variable
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

            // Only check on root
            REQUIRE(num_failed_tests_local == 0);
        }

    }  // MPI gather (variable-sized vector)

}  // TEST_CASE("MPI Utilities")
