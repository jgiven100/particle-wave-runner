#include "particle_poisson.h"

#include <catch2/catch_all.hpp>

#include "mpi_utilities.h"
#include "particle_base.h"

TEST_CASE("Particle Poisson", "[particle_poisson]") {
    // Grab rank and size
    const int rank = pwr::MPIUtilities::Rank();
    const int size = pwr::MPIUtilities::Size();

    REQUIRE(true);

}  // TEST_CASE("Particle Poisson")
