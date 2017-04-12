#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <testthat.h>
#pragma GCC diagnostic pop

#include <Rcpp.h>

#include "../../Core/Random.h"
#include "../../Core/Cell.h"

CATCH_TEST_CASE("Test Cell.h")
{
/*    Random::setSeed(1234);

    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::S4 testCellType = pkgEnv.find("testCellType");
    Rcpp::S4 testDrug = pkgEnv.find("testDrug");

    CellType type = CellType(0, testCellType);
    Drug drug = Drug(0, testDrug);

    Cell cell1 = Cell(&type);
    Cell cell2 = Cell(&type);

    CATCH_REQUIRE(cell1.phase() == INTERPHASE);
    CATCH_REQUIRE(cell1.cellType()->name() == "testType");
    CATCH_REQUIRE(cell1.cycleLength() == 24);
    CATCH_REQUIRE(cell2.cycleLength() == 36);
    CATCH_REQUIRE(!cell1.drugApplied(drug.id()));
    CATCH_REQUIRE(!cell2.drugApplied(drug.id()));

    cell2.setPhase(MITOSIS);
    cell2.setCycleLength(64);
    cell2.markDrugAsApplied(drug.id());
   
    CATCH_REQUIRE(cell2.phase() == MITOSIS);
    CATCH_REQUIRE(cell2.cycleLength() == 64);
    CATCH_REQUIRE(cell2.drugApplied(drug.id()));*/
}




