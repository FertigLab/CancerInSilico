// [[Rcpp::depends(BH)]]

#include "Cell.h"
#include "SpatialHash.h"
#include "Parameters.h"

#include <vector>
#include <testthat.h>
#include <Rcpp.h>
#include <boost/unordered_map.hpp>

#define TEST_APPROX(x) Approx(x).epsilon(0.01)

class TestSpatialHash {

  private:

    SpatialHash *hash;

  public:

    TestSpatialHash(SpatialHash *in_hash) {
        hash = in_hash;
    }

    std::vector<Cell *> *GetCellList() {
        return &(hash->m_cell_list);
    }

    bh_map *GetHashMap() {
        return &(hash->m_hash_map);
    }

    double GetBucketSize() {
        return hash->m_bucket_size;
    }

    Point TestHash(Cell *cell) {
        return hash->Hash(cell);
    }

    void TestRemoveKey(Cell *cell) {
        hash->RemoveKey(cell);
    }

    void TestAddKey(Cell *cell) {
        hash->AddKey(cell);
    }

};

CATCH_TEST_CASE("Test Spatial Hash") {
    CATCH_SECTION("Point Relations") {
        Point p1 = {1, 1};
        Point p2 = {0, 2};
        Point p3 = { -5, 0};
        Point p4 = { -5, 0};
        Point p5 = {0, 0};
        expect_true(p1 < p2);
        expect_true(p2 < p3);
        expect_true(p5 < p1);
        expect_true(p3 == p4);
    }
    CATCH_SECTION("Test Hash Map") {
        Parameters params;
        params.SetMinRadius(1);
        params.SetMaxMigration(10);
        SpatialHash hash = SpatialHash(1);
        TestSpatialHash test_hash = TestSpatialHash(&hash);
        bh_map *cell_map = test_hash.GetHashMap();
        std::vector<Cell *> *cell_list = test_hash.GetCellList();
        CATCH_REQUIRE(test_hash.GetBucketSize() == TEST_APPROX(0.697));
        Cell cell_1 = Cell(std::make_pair(0, 0), &params, 0);
        Cell cell_2 = Cell(std::make_pair(3, 0), &params, 0);
        Cell cell_3 = Cell(std::make_pair(-7, -2), &params, 0);
        Point pt_1 = test_hash.TestHash(&cell_1);
        Point pt_2 = test_hash.TestHash(&cell_2);
        Point pt_3 = test_hash.TestHash(&cell_3);
        hash.Insert(&cell_1);
        hash.Insert(&cell_2);
        hash.Insert(&cell_3);
        CATCH_SECTION("hash cells to points") {
            CATCH_REQUIRE(pt_1.x == TEST_APPROX(0.348));
            CATCH_REQUIRE(pt_1.y == TEST_APPROX(0.348));
            CATCH_REQUIRE(pt_2.x == TEST_APPROX(3.137));
            CATCH_REQUIRE(pt_2.y == TEST_APPROX(0.348));
            CATCH_REQUIRE(pt_3.x == TEST_APPROX(-7.319));
            CATCH_REQUIRE(pt_3.y == TEST_APPROX(-1.743));
        }
        CATCH_SECTION("add cells to hash map") {
            CATCH_REQUIRE(hash.size() == 3);
            CATCH_REQUIRE(cell_map->size() == 3);
            CATCH_REQUIRE(cell_list->size() == 3);
            CATCH_REQUIRE(cell_map->count(pt_1) == 1);
            CATCH_REQUIRE(cell_map->at(pt_2) == &cell_2);
            bh_map::const_iterator it = cell_map->find(pt_3);
            CATCH_REQUIRE(it != cell_map->end());
        }
        test_hash.TestRemoveKey(&cell_1);
        test_hash.TestRemoveKey(&cell_2);
        test_hash.TestRemoveKey(&cell_3);
        CATCH_SECTION("remove keys") {
            CATCH_REQUIRE_THROWS(cell_map->at(pt_1));
            CATCH_REQUIRE_THROWS(test_hash.TestRemoveKey(&cell_3));
            CATCH_REQUIRE(cell_map->size() == 0);
            CATCH_REQUIRE(cell_list->size() == 3);
            CATCH_REQUIRE_THROWS(hash.size());
        }
        test_hash.TestAddKey(&cell_1);
        test_hash.TestAddKey(&cell_2);
        test_hash.TestAddKey(&cell_3);
        CATCH_SECTION("add back keys") {
            CATCH_REQUIRE_THROWS(test_hash.TestAddKey(&cell_1));
            CATCH_REQUIRE_NOTHROW(cell_map->at(pt_2));
            CATCH_REQUIRE(cell_map->size() == 3);
            CATCH_REQUIRE(cell_list->size() == 3);
            CATCH_REQUIRE(hash.size() == 3);
        }
        Cell orig_1 = cell_1;
        Cell orig_2 = cell_2;
        Rcpp::Environment baseEnv("package:base");
        Rcpp::Function setSeed = baseEnv["set.seed"];
        setSeed(40); //makes sure cells move far enough to have new hash value
        cell_1.Migration();
        cell_2.Migration();
        Point new_pt_1 = test_hash.TestHash(&cell_1);
        Point new_pt_2 = test_hash.TestHash(&cell_2);
        CATCH_REQUIRE_NOTHROW(hash.Update(orig_1, cell_1));
        CATCH_REQUIRE_NOTHROW(hash.Update(orig_2, cell_2));
        CATCH_SECTION("update cells") {
            CATCH_REQUIRE_THROWS(hash.Update(orig_1, cell_1));
            CATCH_REQUIRE(cell_map->size() == 3);
            CATCH_REQUIRE(cell_list->size() == 3);
            CATCH_REQUIRE(cell_map->count(pt_1) == 0);
            CATCH_REQUIRE(cell_map->count(pt_2) == 0);
            CATCH_REQUIRE(cell_map->count(new_pt_1) == 1);
            CATCH_REQUIRE(cell_map->count(new_pt_2) == 1);
        }
        CATCH_SECTION("delete cells") {
            CATCH_REQUIRE_NOTHROW(hash.Delete(&cell_2));
            CATCH_REQUIRE_NOTHROW(hash.Delete(&cell_3));
            CATCH_REQUIRE_THROWS(hash.Delete(&cell_3));
            CATCH_REQUIRE(cell_map->size() == 1);
            CATCH_REQUIRE(cell_list->size() == 1);
            CATCH_REQUIRE(cell_map->count(new_pt_2) == 0);
            CATCH_REQUIRE(cell_map->count(pt_3) == 0);
        }
        CATCH_SECTION("get random cell") {
            Cell *p_cell;
            CATCH_REQUIRE_NOTHROW(p_cell = hash.GetRandomCell());
            CATCH_REQUIRE_NOTHROW(hash.Delete(p_cell));
        }
    }
}


