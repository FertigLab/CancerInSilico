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

	params->StoreGrowthRates(std::vector<double> (100,0.1));
	CATCH_REQUIRE_NOTHROW(pop.AddDrug());	

	CATCH_SECTION("Test Add Drug") {

		SpatialHash<Cell>::full_iterator iter = test_pop.hash->begin();
				
		for (; iter != test_pop.hash->end(); ++iter) {
		
			CATCH_REQUIRE((*iter).GetGrowth() == 0.1);

		}

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
	Cell temp = Cell(existing_loc, params);

	CATCH_SECTION("Test ValidCellPlacement") {

		CATCH_REQUIRE(!pop.ValidCellPlacement(&temp));
		existing_loc.x += 1.99;
		temp.SetCoord(existing_loc);
		CATCH_REQUIRE(!pop.ValidCellPlacement(&temp));

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

		for (unsigned int i = 0; i < 100; ++i) {

			CATCH_REQUIRE_NOTHROW(pop.OneTimeStep());

		}

		SpatialHash<Cell>::full_iterator iter1 = test_pop.hash->begin();
		SpatialHash<Cell>::full_iterator iter2 = test_pop.hash->begin();
		
		for (; iter1 != test_pop.hash->end(); ++iter1) {

			for (; iter2 != test_pop.hash->end(); ++iter2) {

				if (*iter1 != *iter2) {

					CATCH_REQUIRE((*iter1).CellDistance(*iter2) > 0);

				}


			}

		}

	}

	delete params;

}
