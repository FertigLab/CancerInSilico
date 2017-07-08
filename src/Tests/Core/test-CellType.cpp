#include <Rcpp.h>

#include "../TestHeader.h"
#include "../../Core/Random.h"
#include "../../Core/CellType.h"

CATCH_TEST_CASE("Test CellType.h")
{
    // load R environment and R model containing sample cell types
    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::S4 model = pkgEnv.find("modDrugs");
    Rcpp::List types = model.slot("cellTypes");

    CATCH_SECTION("default cell properties")
    {
        CellType type0 (0, types(0));
        CATCH_REQUIRE(type0.name() == "DEFAULT");
        CATCH_REQUIRE(type0.id() == 0);
        CATCH_REQUIRE(type0.size() == 1.0);
        CATCH_REQUIRE(type0.minCycle() == 48);
        CATCH_REQUIRE(type0.cycleLength() == 48);
    }

    CATCH_SECTION("double the size of a normal cell")
    {
        CellType type1 (1, types(1));
        CATCH_REQUIRE(type1.name() == "DOUBLE_SIZE");
        CATCH_REQUIRE(type1.id() == 1);
        CATCH_REQUIRE(type1.size() == 2.0);
        CATCH_REQUIRE(type1.minCycle() == 48);
        CATCH_REQUIRE(type1.cycleLength() == 48);
    }
}


