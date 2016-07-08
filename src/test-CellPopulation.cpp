#include <cmath>

#include "Cell.h"
#include "CellPopulation.h"
#include "Parameters.h"
#include "SpatialHash.h"
#include "Point.h"

class TestCellPopulation {

  private:

  public:

	SpatialHash<Cell>* hash;

    TestCellPopulation(CellPopulation* in_pop) {

        hash = &(in_pop->m_population);

    }

    Cell* GetRandomCell() {

        return hash->GetRandomValue();

    }

	void TestAddKey(Point pt, Cell* cell) {
	
		hash->AddKey(pt, cell);

	}

	void TestRemoveKey(Point pt) {
	
		hash->RemoveKey(pt);

	}


};

CATCH_TEST_CASE("Test Cell Population") {

    Parameters *params = new Parameters(pow(2,0.5));

    params->SetEnergyConstant(1);
    params->SetMeanGrowth(0.1);
    params->SetVarGrowth(0.0);
    params->SetMaxMigration(1.0);
    params->SetMaxDeform(0.15);
    params->SetMaxRotate(2.0);
    params->SetResistanceEPSILON(0.05);
    params->SetCompressionDELTA(5.0);

    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(0);

	CellPopulation pop = CellPopulation(params, 100, 0.05);
	TestCellPopulation test_pop(&pop);
	

	CATCH_SECTION("Test Constructor") {

		CATCH_REQUIRE(pop.size() == 100);

	}

	CATCH_SECTION("Test AddDrug") {

		CATCH_REQUIRE_NOTHROW(pop.AddDrug());

	}

	Cell* rand_cell = test_pop.GetRandomCell();

	CATCH_SECTION("Test Get Random Cell") {

		CATCH_REQUIRE(rand_cell->GetRadius() > 0);
		CATCH_REQUIRE_NOTHROW(rand_cell->GetCoord());
		CATCH_REQUIRE(rand_cell->GetAxisLength() >= 0);
		CATCH_REQUIRE(rand_cell->GetAxisAngle() >= 0);
		CATCH_REQUIRE(rand_cell->GetGrowth() > 0);
		CATCH_REQUIRE(!rand_cell->ReadyToDivide());

	}

	Point existing_loc;
	CATCH_REQUIRE_NOTHROW(existing_loc = rand_cell->GetCoord());

	CATCH_SECTION("Test ValidCellPlacement") {

		CATCH_REQUIRE(!pop.ValidCellPlacement(existing_loc));
		existing_loc.x += 1.99;
		CATCH_REQUIRE(!pop.ValidCellPlacement(existing_loc));

	}

	CATCH_SECTION("Test CheckMitosis") {

		while (!rand_cell->ReadyToDivide()) {

			Point old_key = rand_cell->GetCoord();
			rand_cell->DoTrial();
			test_pop.hash->Update(old_key, rand_cell->GetCoord());

		}

		CATCH_REQUIRE(rand_cell->ReadyToDivide());

		Point old_key = rand_cell->GetCoord();
		Cell* daughter_cell = new Cell(rand_cell->Divide());
		test_pop.hash->Insert(daughter_cell->GetCoord(), daughter_cell);
		test_pop.hash->Update(old_key, rand_cell->GetCoord());

		CATCH_REQUIRE_THROWS(test_pop.TestAddKey(rand_cell->GetCoord(), rand_cell));
		CATCH_REQUIRE_THROWS(test_pop.TestAddKey(daughter_cell->GetCoord(), daughter_cell));
		CATCH_REQUIRE_THROWS(test_pop.TestRemoveKey(old_key));
		CATCH_REQUIRE(daughter_cell->GetRadius() > 0);
		CATCH_REQUIRE(daughter_cell->GetGrowth() == rand_cell->GetGrowth());
	
		//delete daughter_cell;	

	}

	CATCH_SECTION("Test CalculateTotalInteraction") {

		CATCH_REQUIRE_NOTHROW(pop.CalculateTotalInteraction(rand_cell));

	}

	CATCH_SECTION("Test AttempTrial") {

		CATCH_REQUIRE_NOTHROW(pop.AttemptTrial(rand_cell));

	}

	CATCH_SECTION("Test Update") {

		CATCH_REQUIRE_NOTHROW(pop.Update());

	}

	CATCH_SECTION("Test OneTimeStep") {

		CATCH_REQUIRE_NOTHROW(pop.OneTimeStep());

	}

	delete params;

}
