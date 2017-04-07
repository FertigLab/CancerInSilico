#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <testthat.h>
#pragma GCC diagnostic pop

#include <Rcpp.h>

#include "../../Core/Random.h"
#include "../../Core/Drug.h"
#include "../../Core/Parameters.h"

CATCH_TEST_CASE("Test Parameters.h")
{
    Random::setSeed(1234);

    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::List rParams = pkgEnv.find("testParams");

    Parameters params (rParams);

    CATCH_REQUIRE(params.initialNum() == 100);
    CATCH_REQUIRE(params.runTime() == 168);
    CATCH_REQUIRE(params.density() == 0.05);
    CATCH_REQUIRE(!params.boundary());
    CATCH_REQUIRE(params.randSeed() == 40);
    CATCH_REQUIRE(params.syncCycles());
    CATCH_REQUIRE(params.outputIncrement() == 0.1);
    CATCH_REQUIRE(params.recordIncrement() == 0.01);
    CATCH_REQUIRE(params.timeIncrement() == 0.001);

    std::vector<Drug>::iterator it = params.drugsBegin();
    CATCH_REQUIRE(it->id() == 0);
    ++it;
    CATCH_REQUIRE(it->id() == 1);
    ++it;
    CATCH_REQUIRE(it->id() == 2);
    ++it;
    CATCH_REQUIRE(it == params.drugsEnd());

    CATCH_REQUIRE(params.randomCellType().name() == "testType1");
}
