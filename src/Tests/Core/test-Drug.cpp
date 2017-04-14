#include <Rcpp.h>

#include "../TestHeader.h"
#include "../../Core/Random.h"
#include "../../Core/Drug.h"

CATCH_TEST_CASE("Test Drug.h")
{
    Random::setSeed(0);

    Rcpp::Environment pkgEnv;
    pkgEnv = Rcpp::Environment::namespace_env("CancerInSilico");
    Rcpp::S4 model = pkgEnv.find("drugsInSystem");

    Rcpp::List types = model.slot("cellTypes");
    Rcpp::List drugs = model.slot("drugs");

    CellType type0 = CellType(0, types(0)); //DEFAULT
    CellType type1 = CellType(1, types(1)); //DOUBLE_SIZE

    Drug drug0 = Drug(0, drugs(0)); // NO_EFFECT
    Drug drug1 = Drug(1, drugs(1)); // HALF_CYCLE_LENGTH
    Drug drug2 = Drug(2, drugs(2)); // HALF_DEFAULT_TYPE
    Drug drug3 = Drug(3, drugs(3)); // ADD_LATE

    CATCH_REQUIRE(drug0.id() == 0);
    CATCH_REQUIRE(drug0.timeAdded() == 0);
    CATCH_REQUIRE(drug1.id() == 1);
    CATCH_REQUIRE(drug1.timeAdded() == 0);
    CATCH_REQUIRE(drug2.id() == 2);
    CATCH_REQUIRE(drug2.timeAdded() == 0);
    CATCH_REQUIRE(drug3.id() == 3);
    CATCH_REQUIRE(drug3.timeAdded() == 6);

    CATCH_REQUIRE(drug0.cycleLengthEffect(type0, 48) == 48);
    CATCH_REQUIRE(drug1.cycleLengthEffect(type0, 48) == 24);
    CATCH_REQUIRE(drug2.cycleLengthEffect(type0, 48) == 24);
    CATCH_REQUIRE(drug3.cycleLengthEffect(type0, 48) == 48);

    CATCH_REQUIRE(drug0.cycleLengthEffect(type1, 48) == 48);
    CATCH_REQUIRE(drug1.cycleLengthEffect(type1, 48) == 24);
    CATCH_REQUIRE(drug2.cycleLengthEffect(type1, 48) == 48);
    CATCH_REQUIRE(drug3.cycleLengthEffect(type1, 48) == 48);
}
