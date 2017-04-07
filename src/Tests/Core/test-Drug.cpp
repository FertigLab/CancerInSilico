#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <testthat.h>
#pragma GCC diagnostic pop

#include <Rcpp.h>

#include "../../Core/Random.h"
#include "../../Core/Drug.h"

CATCH_TEST_CASE("Test Drug.h")
{
    Random::setSeed(1234);

    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::S4 testCellType = pkgEnv.find("testCellType");
    Rcpp::S4 testDrug = pkgEnv.find("testDrug");

    CellType type = CellType(0, testCellType);
    Drug drug = Drug(0, testDrug);

    CATCH_REQUIRE(drug.id() == 0);
    CATCH_REQUIRE(drug.timeAdded() == 10);
    CATCH_REQUIRE(drug.cycleLengthEffect(type, 48, INTERPHASE) == 12);
    CATCH_REQUIRE(drug.cycleLengthEffect(type, 24, INTERPHASE) == 12);
    CATCH_REQUIRE(drug.cycleLengthEffect(type, 24, MITOSIS) == 24);

}


