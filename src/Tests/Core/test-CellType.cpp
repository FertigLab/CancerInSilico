#include <Rcpp.h>

#include "../catch.h"
#include "../../Core/Random.h"
#include "../../Core/CellType.h"

TEST_CASE("Test CellType.h")
{
    // load R environment and R model containing sample cell types
    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::S4 model = pkgEnv.find("modDrugs");
    Rcpp::List types = model.slot("cellTypes");

    SECTION("default cell properties")
    {
        CellType type0 (0, types(0));
        REQUIRE(type0.name() == "DEFAULT");
        REQUIRE(type0.id() == 0);
        REQUIRE(type0.size() == 1.0);
        REQUIRE(type0.minCycle() == 48);
        REQUIRE(type0.cycleLength() == 48);
    }

    SECTION("double the size of a normal cell")
    {
        CellType type1 (1, types(1));
        REQUIRE(type1.name() == "DOUBLE_SIZE");
        REQUIRE(type1.id() == 1);
        REQUIRE(type1.size() == 2.0);
        REQUIRE(type1.minCycle() == 48);
        REQUIRE(type1.cycleLength() == 48);
    }
}


