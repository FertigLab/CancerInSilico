#include <Rcpp.h>
#include <testthat.h>

#include "SquareLattice.h"
#include "Point.h"

#define TEST_APPROX(x) Approx(x).epsilon(0.01)

struct TestObject
{
    Point coord;
    int value;

    TestObject(Point c, int v) {
        coord = c;
        value = v;
    }

};

CATCH_TEST_CASE("Test SquareLattice.h (and Lattice.h) with doubles")
{
    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(20);

    SquareLattice<double> testLat1 (1.0);

    // test insert
    CATCH_REQUIRE_NOTHROW(testLat1.insert(Point(0.5, 0), 0.0));
    CATCH_REQUIRE_NOTHROW(testLat1.insert(Point(1.5, 0), 1.0));
    CATCH_REQUIRE_NOTHROW(testLat1.insert(Point(2.5, 0), 2.0));

    CATCH_REQUIRE(testLat1.size() == 3);
    CATCH_REQUIRE_THROWS(testLat1.insert(Point(0,0), 1.0));

    /// test iterator
    double sum = 0.0;
    SquareLattice<double>::iterator it = testLat1.begin();
    for (; it != testLat1.end(); ++it) {sum += *it;}
    CATCH_REQUIRE(sum == 3.0);

    // test erase
    CATCH_REQUIRE_NOTHROW(testLat1.erase(Point(1.5, 0)));
    CATCH_REQUIRE_THROWS(testLat1.erase(Point(1.5, 0)));
    CATCH_REQUIRE(testLat1.size() == 2);

    sum = 0.0;
    it = testLat1.begin();
    for (; it != testLat1.end(); ++it) {sum += *it;}
    CATCH_REQUIRE(sum == 2.0);

    // test update
    CATCH_REQUIRE_NOTHROW(testLat1.update(Point(2.5, 0), Point(0, 5)));
    CATCH_REQUIRE(testLat1.size() == 2);
    CATCH_REQUIRE_THROWS(testLat1.update(Point(2.5, 0), Point(10, 5)));
/*    CATCH_REQUIRE_THROWS(testLat1.erase(Point(2.5, 0)));
/*    CATCH_REQUIRE_THROWS(testLat1.insert(Point(0, 5), 7.0));

    sum = 0.0;
    it = testLat1.begin();
    for (; it != testLat1.end(); ++it) {sum += *it;}
    CATCH_REQUIRE(sum == 2.0);

    // test randomValue
    CATCH_REQUIRE(testLat1.randomValue() == 2.0);
    testLat1.randomValue() = 4.0;

    sum = 0.0;
    it = testLat1.begin();
    for (; it != testLat1.end(); ++it) {sum += *it;}
    CATCH_REQUIRE(sum == 6.0);*/

}

CATCH_TEST_CASE("Test SquareLattice.h (and Lattice.h) with objects")
{
    SquareLattice<TestObject> testLat2 (0.5);

    TestObject obj0 (Point(0, 0), 0);
    TestObject obj1 (Point(0.75, 0), 1);
    TestObject obj2 (Point(1.25, 0), 2);
    TestObject obj3 (Point(1.75, 0), 3);

    CATCH_REQUIRE_NOTHROW(testLat2.insert(obj0.coord, obj0));
    CATCH_REQUIRE_NOTHROW(testLat2.insert(obj1.coord, obj1));
    CATCH_REQUIRE_NOTHROW(testLat2.insert(obj2.coord, obj2));
    CATCH_REQUIRE_NOTHROW(testLat2.insert(obj3.coord, obj3));
}

