#include <Rcpp.h>

#include "../TestHeader.h"
#include "../../Core/Random.h"
#include "../../Core/CellType.h"

CATCH_TEST_CASE("Test CellType.h")
{
    Random::setSeed(0);

    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::S4 model = pkgEnv.find("modCellTypes");

    Rcpp::List types = model.slot("cellTypes");

    CellType type0 (0, types(0)); //DEFAULT
    CellType type1 (1, types(1)); //DOUBLE_SIZE
    CellType type2 (2, types(2)); //SHORT_CYCLE
    CellType type3 (3, types(3)); //RANDOM_CYCLE

    CATCH_REQUIRE(type0.name() == "DEFAULT");
    CATCH_REQUIRE(type0.id() == 0);
    CATCH_REQUIRE(type0.size() == 1.0);
    CATCH_REQUIRE(type0.minCycle() == 48);
    CATCH_REQUIRE(type0.cycleLength() == 48);

    CATCH_REQUIRE(type1.name() == "DOUBLE_SIZE");
    CATCH_REQUIRE(type1.id() == 1);
    CATCH_REQUIRE(type1.size() == 2.0);
    CATCH_REQUIRE(type1.minCycle() == 48);
    CATCH_REQUIRE(type1.cycleLength() == 48);

    CATCH_REQUIRE(type2.name() == "SHORT_CYCLE");
    CATCH_REQUIRE(type2.id() == 2);
    CATCH_REQUIRE(type2.size() == 1.0);
    CATCH_REQUIRE(type2.minCycle() == 6);
    CATCH_REQUIRE(type2.cycleLength() == 6);

    CATCH_REQUIRE(type3.name() == "RANDOM_CYCLE");
    CATCH_REQUIRE(type3.id() == 3);
    CATCH_REQUIRE(type3.size() == 1.0);
    CATCH_REQUIRE(type3.minCycle() == 10);

    double min = 20.0, max = 10.0;
    for (unsigned i = 0; i < 10000; ++i)
    {
        min = std::min(type3.cycleLength(), min);
        max = std::max(type3.cycleLength(), max);
    }
}


