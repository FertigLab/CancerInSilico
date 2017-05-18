#include <Rcpp.h>

#include "../TestHeader.h"
#include "../../Core/Random.h"
#include "../../Core/CellType.h"

CATCH_TEST_CASE("Test CellType.h")
{
    // load R environment and R model containing sample cell types
    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::S4 model = pkgEnv.find("modCellTypes");
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

    CATCH_SECTION("different, but constant, cycle length")
    {
        CellType type2 (2, types(2));
        CATCH_REQUIRE(type2.name() == "SHORT_CYCLE");
        CATCH_REQUIRE(type2.id() == 2);
        CATCH_REQUIRE(type2.size() == 1.0);
        CATCH_REQUIRE(type2.minCycle() == 6);
        CATCH_REQUIRE(type2.cycleLength() == 6);
    }

    CATCH_SECTION("random cycle length")
    {
        CellType type3 (3, types(3));
        CATCH_REQUIRE(type3.name() == "RANDOM_CYCLE");
        CATCH_REQUIRE(type3.id() == 3);
        CATCH_REQUIRE(type3.size() == 1.0);
        CATCH_REQUIRE(type3.minCycle() == 10);

        Random::setSeed(0);
        double min = 20.0, max = 10.0;
        for (unsigned i = 0; i < 10000; ++i)
        {
            min = std::min(type3.cycleLength(), min);
            max = std::max(type3.cycleLength(), max);
        }
        CATCH_REQUIRE(min >= 10.0); // test bounds of cycle length
        CATCH_REQUIRE(max <= 20.0); // distribution
    }
}


