#include <Rcpp.h>
#include <testthat.h>

#include "SquareLattice.h"
#include "Point.h"

struct TestObject
{
    Point<double> coords;
    int val;

    TestObject(double x, double y, int v)
        : coords(Point<double>(x,y)), val(v)
    {}
};

CATCH_TEST_CASE("Test SquareLattice.h (and Lattice.h) with doubles")
{
    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(20);

    SquareLattice<TestObject> testLat1 (1.0);

    // test insert
    CATCH_REQUIRE_NOTHROW(testLat1.insert(TestObject(0.5, 0, 0)));
    CATCH_REQUIRE_NOTHROW(testLat1.insert(TestObject(1.5, 0, 1)));
    CATCH_REQUIRE_NOTHROW(testLat1.insert(TestObject(2.5, 0, 2)));

    CATCH_REQUIRE(testLat1.size() == 3);
    CATCH_REQUIRE_NOTHROW(testLat1.insert(TestObject(0.5, 0, 0)));

    /// test iterator
    int sum = 0;
    SquareLattice<TestObject>::iterator it = testLat1.begin();
    for (; it != testLat1.end(); ++it) {sum += (*it).val;}
    CATCH_REQUIRE(sum == 3);

    // test erase
    CATCH_REQUIRE_NOTHROW(testLat1.erase(Point<double>(1.5, 0)));
    CATCH_REQUIRE_THROWS(testLat1.erase(Point<double>(1.5, 0)));
    CATCH_REQUIRE(testLat1.size() == 2);

    sum = 0;
    it = testLat1.begin();
    for (; it != testLat1.end(); ++it) {sum += (*it).val;}
    CATCH_REQUIRE(sum == 2);

    // test update
    CATCH_REQUIRE_NOTHROW(testLat1.update(Point<double>(2.5, 0),
        Point<double>(0, 5)));
    CATCH_REQUIRE(testLat1.size() == 2);
    CATCH_REQUIRE_THROWS(testLat1.erase(Point<double>(2.5, 0)));

    sum = 0;
    it = testLat1.begin();
    for (; it != testLat1.end(); ++it) {sum += (*it).val;}
    CATCH_REQUIRE(sum == 2);

    // test randomValue
    CATCH_REQUIRE(testLat1.randomValue() == 2);
    testLat1.randomValue().val = 4;

    sum = 0;
    it = testLat1.begin();
    for (; it != testLat1.end(); ++it) {sum += (*it).val;}
    CATCH_REQUIRE(sum == 4);
}


