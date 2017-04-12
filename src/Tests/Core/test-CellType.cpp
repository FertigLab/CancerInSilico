#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <testthat.h>
#pragma GCC diagnostic pop

#include <Rcpp.h>

#include "../../Core/Random.h"
#include "../../Core/CellType.h"

CATCH_TEST_CASE("Test CellType.h")
{
/*    Random::setSeed(1234);

    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::S4 testCellType = pkgEnv.find("testCellType");

    CellType type = CellType(0, testCellType);

    CATCH_REQUIRE(type.id() == 0);
    CATCH_REQUIRE(type.size() == 2.0);
    CATCH_REQUIRE(type.name() == "TestType");
    CATCH_REQUIRE(type.cycleLength() == 36);
    CATCH_REQUIRE(type.minCycleLength() == 12);*/
}


