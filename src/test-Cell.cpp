#include <testthat.h>
#include <Rcpp.h>

#include "Cell.h"
#include "Parameters.h"
#include "Point.h"

#define TEST_APPROX(x) Approx(x).epsilon(0.01)

CATCH_TEST_CASE("Test Cell") {

	Parameters params = Parameters(pow(2,0.5));

	params.SetMeanGrowth(0.03);
	params.SetVarGrowth(0.0);
	params.SetMaxMigration(0.3);
	params.SetMaxRotate(2.0);
	params.SetMaxDeform(0.1);

	std::vector<Cell> cells;

    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(40);
	cells.push_back(Cell(Point(0,0), &params));
	
	CATCH_SECTION("new cell constructor works properly") {
	
		CATCH_REQUIRE(cells[0].GetCoord() == Point(0,0));
		CATCH_REQUIRE(cells[0].GetRadius() == TEST_APPROX(1.28));
		CATCH_REQUIRE(cells[0].GetAxisLength() == 0);
		CATCH_REQUIRE(cells[0].GetAxisAngle() == 0);
		CATCH_REQUIRE(cells[0].GetGrowth() == 0.03);
		CATCH_REQUIRE(!cells[0].ReadyToDivide());
		
	}

	cells.push_back(Cell(Point(6,8), &params));
	
	CATCH_SECTION("test simple distance calculation") {

		CATCH_REQUIRE(cells[0].CellDistance(cells[1]) ==
					10 - cells[0].GetRadius() - cells[1].GetRadius());

	}

	double r0, r1;

	while (cells[0].GetRadius() < params.GetMaxRadius()
			|| cells[1].GetRadius() < params.GetMaxRadius()) {

		r0 = cells[0].GetRadius();
		r1 = cells[1].GetRadius();
		cells[0].Growth();
		cells[1].Growth();
		CATCH_REQUIRE(cells[0].GetRadius() < r0 + cells[0].GetGrowth());
		CATCH_REQUIRE(cells[1].GetRadius() < r1 + cells[1].GetGrowth());

	}

	CATCH_SECTION("test cell growth") {

		CATCH_REQUIRE(cells[0].GetRadius() == params.GetMaxRadius());
		CATCH_REQUIRE(cells[1].GetRadius() == params.GetMaxRadius());

	}

	CATCH_SECTION("test cell movement") {

		cells[0].Migration();
		cells[1].Migration();
		CATCH_REQUIRE(cells[0].GetCoord() == Point(-0.05895,-0.15666));
		CATCH_REQUIRE(cells[1].GetCoord() == Point(6.012694,7.75392));

	}

	while (cells[0].GetRadius() > 1.1 || cells[1].GetRadius() > 1.1) {

		r0 = cells[0].GetRadius();
		r1 = cells[1].GetRadius();
		cells[0].Deformation();
		cells[1].Deformation();
		CATCH_REQUIRE(cells[0].GetRadius() > r0 - params.GetMaxDeform());
		CATCH_REQUIRE(cells[1].GetRadius() > r1 - params.GetMaxDeform());

	}

	CATCH_SECTION("test cell deformation") {
	
		CATCH_REQUIRE(cells[0].GetRadius() <= 1.1);
		CATCH_REQUIRE(cells[1].GetRadius() <= 1.1);
		CATCH_REQUIRE(cells[0].GetAxisLength() >= 3.3);
		CATCH_REQUIRE(cells[1].GetAxisLength() >= 3.3);
		CATCH_REQUIRE(cells[0].GetAxisAngle() == 0);
		CATCH_REQUIRE(cells[1].GetAxisAngle() == 0);

	}

	setSeed(15);
	cells[0].Rotation();
	cells[1].Rotation();

	CATCH_SECTION("test cell rotation") {

		CATCH_REQUIRE(cells[0].GetAxisAngle() == TEST_APPROX(0.408));
		CATCH_REQUIRE(cells[1].GetAxisAngle() == TEST_APPROX(-1.220));

	}

	params.SetMaxMigration(2.0);
	cells[0].Migration();
	cells[1].Migration();

	CATCH_SECTION("test dividing cell movement") {

		CATCH_REQUIRE(cells[0].GetCoord() == Point(-1.14662,-1.59722));
		CATCH_REQUIRE(cells[1].GetCoord() == Point(7.20876,7.91525));

	}

	CATCH_SECTION("test complicated distance calculation") {

		CATCH_REQUIRE(cells[0].CellDistance(cells[1])
						== TEST_APPROX(9.36));

	}
	
	while (!cells[0].ReadyToDivide()) {

		cells[0].Deformation();

	}

	cells.push_back(cells[0].Divide());

	CATCH_SECTION("test cell division") {

		CATCH_REQUIRE(cells[0].GetCoord() == Point(-0.22888,-1.20002));
		CATCH_REQUIRE(cells[0].GetRadius() == 1);
		CATCH_REQUIRE(cells[0].GetAxisLength() == 0);
		CATCH_REQUIRE(cells[0].GetAxisAngle() == 0);
		CATCH_REQUIRE(cells[0].GetGrowth() == 0.03);
		CATCH_REQUIRE(!cells[0].ReadyToDivide());

		CATCH_REQUIRE(cells[2].GetCoord() == Point(-2.06435,-1.99441));
		CATCH_REQUIRE(cells[2].GetRadius() == 1);
		CATCH_REQUIRE(cells[2].GetAxisLength() == 0);
		CATCH_REQUIRE(cells[2].GetAxisAngle() == 0);
		CATCH_REQUIRE(cells[2].GetGrowth() == 0.03);
		CATCH_REQUIRE(!cells[2].ReadyToDivide());

		CATCH_REQUIRE(cells[0].CellDistance(cells[2]) == 0);

	}

}
