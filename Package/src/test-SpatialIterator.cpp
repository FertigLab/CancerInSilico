#include "Cell.hpp"
#include "SpatialHash.hpp"
#include "Parameters.hpp"

#include <Rcpp.h>
#include <testthat.h>

class TestSpatialIterator {

private:

  SpatialIterator* iter;

public:

  TestSpatialIterator(SpatialIterator* in_iter) {

    iter = in_iter;

  }

  std::vector<Point>* GetPointList() {

    return &(iter->search_points);

  }

};

CATCH_TEST_CASE("Test Spatial Iterator") {

    Parameters params;
    params.SetMinRadius(1);

    SpatialHash hash = SpatialHash(1);
    SpatialIterator full_iter = hash.getFullIterator();
    TestSpatialIterator test_full_iter = TestSpatialIterator(&full_iter);

    CATCH_SECTION("Test empty iterator") {
    
      CATCH_REQUIRE(!full_iter.Next());     
      CATCH_REQUIRE(test_full_iter.GetPointList()->empty());

    }
    
    Cell* cell_1 = new Cell(std::make_pair(0,0), &params, 0);  
    hash.Insert(cell_1);

    SpatialIterator circ_iter = hash.getCircularIterator(cell_1,1);
    TestSpatialIterator test_circ_iter = TestSpatialIterator(&circ_iter);
  
    full_iter = hash.getFullIterator();

    CATCH_SECTION("Test single cell iterator") {

      CATCH_REQUIRE(full_iter.Next());
      CATCH_REQUIRE(!full_iter.Next());

      CATCH_REQUIRE(test_full_iter.GetPointList()->size() == 1);

      CATCH_REQUIRE(test_circ_iter.GetPointList()->empty());
      CATCH_REQUIRE(!circ_iter.Next());

    }

    for (int i = 1; i < 10; ++i) {

      CATCH_REQUIRE_NOTHROW(hash.Insert(new Cell(std::make_pair(i,i), &params, 0)));

    }

    full_iter = hash.getFullIterator();
    circ_iter = hash.getCircularIterator(cell_1,13);
    Cell* temp_cell;

    CATCH_SECTION("Test 10 cell iterator") {
      
      CATCH_SECTION("Test full iterator") {

        int count = 0;
        while (full_iter.Next()) {

          CATCH_REQUIRE_NOTHROW(temp_cell = full_iter.getCell());
          count++;
        
        }

        CATCH_REQUIRE(count == hash.size());

      }

      CATCH_SECTION("Test circ iterator") {

        int count = 0;
        while (circ_iter.Next()) {

          CATCH_REQUIRE_NOTHROW(temp_cell = circ_iter.getCell());
          count++;
        
        }

        CATCH_REQUIRE(count == hash.size() - 1);

      }

    }

    for (int i = 10; i < 50; ++i) {

      CATCH_REQUIRE_NOTHROW(hash.Insert(new Cell(std::make_pair(i,i), &params, 0)));

    }

    full_iter = hash.getFullIterator();
    circ_iter = hash.getCircularIterator(cell_1,40);

    CATCH_SECTION("Test 50 cell iterator") {

      CATCH_SECTION("Test full iterator") {

        int count = 0;
        while (full_iter.Next()) {

          CATCH_REQUIRE_NOTHROW(temp_cell = full_iter.getCell());
          count++;
        
        }

        CATCH_REQUIRE(count == hash.size());

      }

      CATCH_SECTION("Test circ iterator") {

        int count = 0;
        while (circ_iter.Next()) {

          CATCH_REQUIRE_NOTHROW(temp_cell = circ_iter.getCell());
          count++;
        
        }

        CATCH_REQUIRE(count == 29);

      }

    }

}




