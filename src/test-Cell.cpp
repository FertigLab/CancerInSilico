#include <testthat.h>
#include <Rcpp.h>

#include "Cell.h"
#include "Parameters.h"
#include "Point.h"

#define TEST_APPROX(x) Approx(x).epsilon(0.01)

CATCH_TEST_CASE("Test Cell") {

    Rcpp::Environment env;
    env = Rcpp::Environment::namespace_env("CancerInSilico");

    Rcpp::List Rparams = env.find("testParams");

    CATCH_REQUIRE(Rparams.size() == 17);

    Parameters* params;

    CATCH_REQUIRE_NOTHROW(params = new Parameters(pow(2, 0.5), Rparams));

    std::vector<Cell> cells;

    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(40);

    cells.push_back(Cell(Point(0,0), params));
    CATCH_REQUIRE(cells[0].GetGrowth() == TEST_APPROX(0.000334));
    cells[0].SetGrowth(0.03);    

    CATCH_SECTION("new cell constructor works properly") {
    
        CATCH_REQUIRE(cells[0].GetCoord() == Point(0,0));
        CATCH_REQUIRE(cells[0].GetRadius() == TEST_APPROX(1.0));
        CATCH_REQUIRE(cells[0].GetAxisLength() == TEST_APPROX(2.0));
        CATCH_REQUIRE(cells[0].GetAxisAngle() == TEST_APPROX(4.30));
        CATCH_REQUIRE(cells[0].GetGrowth() == 0.03);
        CATCH_REQUIRE(!cells[0].ReadyToDivide());
        
    }

    cells.push_back(Cell(Point(6,8), params));
    
    CATCH_SECTION("test simple distance calculation") {

        CATCH_REQUIRE(cells[0].CellDistance(cells[1]) ==
                    10 - cells[0].GetRadius() - cells[1].GetRadius());

    }

    double r0, r1, a0, a1;

    while (cells[0].GetRadius() < params->maxRadius()
            || cells[1].GetRadius() < params->maxRadius()) {

        r0 = cells[0].GetRadius();
        r1 = cells[1].GetRadius();
        cells[0].Growth();
        cells[1].Growth();
        CATCH_REQUIRE(cells[0].GetRadius() < r0 + cells[0].GetGrowth());
        CATCH_REQUIRE(cells[1].GetRadius() < r1 + cells[1].GetGrowth());

    }

    CATCH_SECTION("test cell growth") {

        CATCH_REQUIRE(cells[0].GetRadius() == params->maxRadius());
        CATCH_REQUIRE(cells[1].GetRadius() == params->maxRadius());

    }

    CATCH_SECTION("test cell movement") {

        cells[0].Translation();
        cells[1].Translation();

    }

    while (cells[0].GetAxisLength() < 4.0 || cells[1].GetAxisLength() < 4.0) {

        a0 = cells[0].GetAxisLength();
        a1 = cells[1].GetAxisLength();
        cells[0].Deformation();
        cells[1].Deformation();
        CATCH_REQUIRE(cells[0].GetAxisLength() > a0 - params->maxDeform());
        CATCH_REQUIRE(cells[1].GetAxisLength() > a1 - params->maxDeform());

    }

    CATCH_SECTION("test cell deformation") {
    
        CATCH_REQUIRE(cells[0].GetRadius() == TEST_APPROX(1.0));
        CATCH_REQUIRE(cells[1].GetRadius() == TEST_APPROX(1.0));
        CATCH_REQUIRE(cells[0].GetAxisAngle() == TEST_APPROX(5.60));
        CATCH_REQUIRE(cells[1].GetAxisAngle() == TEST_APPROX(5.94));

    }

    setSeed(15);
    cells[0].Rotation();
    cells[1].Rotation();

    CATCH_SECTION("test cell rotation") {

        CATCH_REQUIRE(cells[0].GetAxisAngle() == TEST_APPROX(5.66));
        CATCH_REQUIRE(cells[1].GetAxisAngle() == TEST_APPROX(5.76));

    }

    cells[0].Translation();
    cells[1].Translation();

    CATCH_SECTION("test dividing cell movement") {

        CATCH_REQUIRE(cells[0].GetCoord() == Point(-0.05733,-0.07986));
        CATCH_REQUIRE(cells[1].GetCoord() == Point(6.06044,7.99576));

    }

    CATCH_SECTION("test complicated distance calculation") {

        CATCH_REQUIRE(cells[0].CellDistance(cells[1])
                        == TEST_APPROX(8.03));

    }
    
    while (!cells[0].ReadyToDivide()) {

        cells[0].Deformation();

    }

    cells.push_back(cells[0].Divide(false));

    CATCH_SECTION("test cell division") {

        CATCH_REQUIRE(cells[0].GetRadius() == 1);
        CATCH_REQUIRE(cells[0].GetAxisLength() == 2.0);
        CATCH_REQUIRE(cells[0].GetAxisAngle() == TEST_APPROX(5.122));
        CATCH_REQUIRE(cells[0].GetGrowth() == 0.03);
        CATCH_REQUIRE(!cells[0].ReadyToDivide());

        CATCH_REQUIRE(cells[2].GetRadius() == 1);
        CATCH_REQUIRE(cells[2].GetAxisLength() == 2.0);
        CATCH_REQUIRE(cells[2].GetAxisAngle() == TEST_APPROX(1.596));
        CATCH_REQUIRE(cells[2].GetGrowth() == TEST_APPROX(0.000333));
        CATCH_REQUIRE(!cells[2].ReadyToDivide());

        CATCH_REQUIRE(cells[0].CellDistance(cells[2]) == TEST_APPROX(0));

    }

}
