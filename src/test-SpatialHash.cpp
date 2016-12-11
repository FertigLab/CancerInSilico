// [[Rcpp::depends(BH)]]

#include <R.h>
#include <vector>
#include <testthat.h>
#include <boost/unordered_map.hpp>

#include "SpatialHash.h"
#include "Point.h"

#define TEST_APPROX(x) Approx(x).epsilon(0.01)

typedef struct test_object {

    Point coord;
    int value;

    test_object(Point c, int v) {
        coord = c;
        value = v;
    }

} TestObject;

typedef boost::unordered_map<Point, TestObject*, ihash, iequal_to> bh_map;

class TestSpatialHash {

  private:

    SpatialHash<TestObject>* hash;

  public:

    TestSpatialHash(SpatialHash<TestObject>* in_hash) {

        hash = in_hash;

    }

    std::vector<TestObject*>* GetValueList() {

        return &(hash->m_value_list);

    }

    bh_map* GetHashMap() {

        return &(hash->m_hash_map);

    }

    double GetBucketSize() {

        return hash->m_bucket_size;

    }

    Point TestHash(TestObject t) {

        return hash->Hash(t.coord);

    }

    void TestRemoveKey(TestObject t) {

        hash->RemoveKey(t.coord);

    }

    void TestAddKey(TestObject t) {

        hash->AddKey(t.coord, &t);

    }

};

CATCH_TEST_CASE("Test Spatial Hash and its components") {

    CATCH_SECTION("Point Functions") {

        Point p1 = Point(1,1);
        Point p2 = Point(0,2);
        Point p3 = Point(-5,0);
        Point p4 = Point(-5,0);
        Point p5 = Point(0,0);
        expect_true(p1 < p2);
        expect_true(p2 < p3);
        expect_true(p5 != p1);
        expect_true(p3 == p4);

    }

    CATCH_SECTION("Test Hash Map") {

        SpatialHash<TestObject> hash = SpatialHash<TestObject>(1.0);
        TestSpatialHash test_hash = TestSpatialHash(&hash);

        bh_map* internal_map = test_hash.GetHashMap();
        std::vector<TestObject*>* internal_list = test_hash.GetValueList();

        CATCH_REQUIRE(test_hash.GetBucketSize() == TEST_APPROX(0.697));

        TestObject obj_1 = TestObject(Point(0,0), 1);
        TestObject obj_2 = TestObject(Point(3, 0), 2);
        TestObject obj_3 = TestObject(Point(-7, -2), 3);

        Point pt_1 = test_hash.TestHash(obj_1);
        Point pt_2 = test_hash.TestHash(obj_2);
        Point pt_3 = test_hash.TestHash(obj_3);

        hash.Insert(obj_1.coord, &obj_1);
        hash.Insert(obj_2.coord, &obj_2);
        hash.Insert(obj_3.coord, &obj_3);

        CATCH_SECTION("hash function") {

            CATCH_REQUIRE(pt_1.x == 0);
            CATCH_REQUIRE(pt_1.y == 0);
            CATCH_REQUIRE(pt_2.x == 4);
            CATCH_REQUIRE(pt_2.y == 0);
            CATCH_REQUIRE(pt_3.x == -10);
            CATCH_REQUIRE(pt_3.y == -3);

        }

       CATCH_SECTION("add objects to hash map") {

            CATCH_REQUIRE(hash.size() == 3);
            CATCH_REQUIRE(internal_map->size() == 3);
            CATCH_REQUIRE(internal_list->size() == 3);
            CATCH_REQUIRE(internal_map->count(pt_1) == 1);
            CATCH_REQUIRE(internal_map->at(pt_2) == &obj_2);
            bh_map::const_iterator it = internal_map->find(pt_3);
            CATCH_REQUIRE(it != internal_map->end());

        }

        test_hash.TestRemoveKey(obj_1);
        test_hash.TestRemoveKey(obj_2);
        test_hash.TestRemoveKey(obj_3);

        CATCH_SECTION("remove keys") {

            CATCH_REQUIRE_THROWS(internal_map->at(pt_1));
            CATCH_REQUIRE_THROWS(test_hash.TestRemoveKey(obj_3));
            CATCH_REQUIRE(internal_map->size() == 0);
            CATCH_REQUIRE(internal_list->size() == 3);

        }

        test_hash.TestAddKey(obj_1);
        test_hash.TestAddKey(obj_2);
        test_hash.TestAddKey(obj_3);

        CATCH_SECTION("add back keys") {

            CATCH_REQUIRE_THROWS(test_hash.TestAddKey(obj_1));
            CATCH_REQUIRE_NOTHROW(internal_map->at(pt_2));
            CATCH_REQUIRE(internal_map->size() == 3);
            CATCH_REQUIRE(internal_list->size() == 3);
            CATCH_REQUIRE(hash.size() == 3);

        }

        Point orig_1 = obj_1.coord;
        Point orig_2 = obj_2.coord;

        obj_1.coord.x = -43.1;
        obj_1.coord.y = 0.02;

        obj_2.coord.x = 12.2;
        obj_2.coord.y = -4.87;

        Point new_1 = test_hash.TestHash(obj_1);
        Point new_2 = test_hash.TestHash(obj_2);

        CATCH_REQUIRE_NOTHROW(hash.Update(orig_1, obj_1.coord));
        CATCH_REQUIRE_NOTHROW(hash.Update(orig_2, obj_2.coord));

        CATCH_SECTION("update objects") {

            CATCH_REQUIRE_THROWS(hash.Update(orig_1, obj_1.coord));
            CATCH_REQUIRE(internal_map->size() == 3);
            CATCH_REQUIRE(internal_list->size() == 3);
            CATCH_REQUIRE(internal_map->count(pt_1) == 0);
            CATCH_REQUIRE(internal_map->count(pt_2) == 0);
            CATCH_REQUIRE(internal_map->count(new_1) == 1);
            CATCH_REQUIRE(internal_map->count(new_2) == 1);

        }

        CATCH_SECTION("delete objects") {

            CATCH_REQUIRE_NOTHROW(hash.Delete(obj_2.coord, &obj_2));
            CATCH_REQUIRE_NOTHROW(hash.Delete(obj_3.coord, &obj_3));
            CATCH_REQUIRE_THROWS(hash.Delete(obj_3.coord, &obj_3));
            CATCH_REQUIRE(internal_map->size() == 1);
            CATCH_REQUIRE(internal_list->size() == 1);
            CATCH_REQUIRE(internal_map->count(new_2) == 0);
            CATCH_REQUIRE(internal_map->count(pt_3) == 0);

        }

        CATCH_SECTION("get random object") {

            TestObject *p_obj;
            CATCH_REQUIRE(hash.size() == 3);
            CATCH_REQUIRE_NOTHROW(p_obj = hash.GetRandomValue());
            CATCH_REQUIRE_NOTHROW(hash.Delete(p_obj->coord,p_obj));

        }

    }

    CATCH_SECTION("Test SpatialHash Iterator") {

        SpatialHash<TestObject> hash = SpatialHash<TestObject>(1.0);
        SpatialHash<TestObject>::full_iterator full_iter = hash.begin();

        CATCH_SECTION("Test empty iterator") {

            full_iter == hash.end();
            CATCH_REQUIRE_NOTHROW(TestObject temp = *full_iter);

        }

        TestObject* first = new TestObject(Point(0,0),0);
        CATCH_REQUIRE_NOTHROW(hash.Insert(first->coord, first));

        SpatialHash<TestObject>::circular_iterator circ_iter = hash.begin(first->coord, 0);
        full_iter = hash.begin();
        
        CATCH_SECTION("Test single object iterator: initialization") {

            CATCH_REQUIRE(full_iter != hash.end());
            CATCH_REQUIRE(circ_iter == hash.end(first->coord, 0));

        }

        CATCH_SECTION("Test single object iterator: pre-increment") {
    
            CATCH_REQUIRE(++full_iter == hash.end());

        }
        
        CATCH_SECTION("Test single object iterator: post-increment") {
    
            CATCH_REQUIRE(full_iter++ != hash.end());
            CATCH_REQUIRE(full_iter == hash.end());

        }

        double rad = 2.01, ang;
        TestObject* t;

        for (unsigned int i = 1; i <= 10; i++) {

            ang = R::runif(0, 2 * M_PI);            
            t = new TestObject(Point(rad * cos(ang), rad * sin(ang)), i);
            rad += 2.01;
            CATCH_REQUIRE_NOTHROW(hash.Insert(t->coord, t));

        }

        full_iter = hash.begin();

        CATCH_SECTION("Test 10 object iterator - full iterator") {

            int count = 0;
            for (; full_iter != hash.end(); ++full_iter) {

                count++;
                CATCH_REQUIRE_NOTHROW(*full_iter);
                CATCH_REQUIRE_NOTHROW(&full_iter);

            }
            CATCH_REQUIRE(count == 11);

        }

        CATCH_SECTION("Test 10 object iterator - circular iterator") {

            int count = 0;
            circ_iter = hash.begin(first->coord, 2.1);
            for (; circ_iter != hash.end(first->coord, 2.1); ++circ_iter) {

                count++;
                CATCH_REQUIRE_NOTHROW(*circ_iter);
                CATCH_REQUIRE_NOTHROW(&circ_iter);

            }
            CATCH_REQUIRE(count >= 1);

            count = 0;
            circ_iter = hash.begin(first->coord, 4.1);
            for (; circ_iter != hash.end(first->coord, 4.1); ++circ_iter) {

                count++;
                CATCH_REQUIRE_NOTHROW(*circ_iter);
                CATCH_REQUIRE_NOTHROW(&circ_iter);

            }
            CATCH_REQUIRE(count >= 2);

            count = 0;
            circ_iter = hash.begin(first->coord, 20.2);
            for (; circ_iter != hash.end(first->coord, 20.2); ++circ_iter) {

                count++;
                CATCH_REQUIRE_NOTHROW(*circ_iter);
                CATCH_REQUIRE_NOTHROW(&circ_iter);

            }
            CATCH_REQUIRE(count == 10);
        
        }

        SpatialHash<TestObject>::full_iterator del_iter = hash.begin();
        for (; del_iter != hash.end(); ++del_iter) {

            delete &del_iter;            
    
        }

    }

}

