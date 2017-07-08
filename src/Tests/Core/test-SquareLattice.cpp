#include <Rcpp.h>

#include "../TestHeader.h"
#include "../../Core/SquareLattice.h"

struct TestObject
{
    Point<double> coords;
    int val;

    TestObject(double x, double y, int v)
        : coords(Point<double>(x,y)), val(v)
    {}
};

#define TEST_OBJ(x,y,v) Point<double>(x,y), TestObject(x,y,v)

CATCH_TEST_CASE("Test SquareLattice.h")
{
    Random::setSeed(0);

    CATCH_SECTION("Square Lattice Constructor - default")
    {
        SquareLattice<TestObject> testLat;

        CATCH_REQUIRE(testLat.size() == 0);
        CATCH_REQUIRE_THROWS(testLat.randomValue());
        CATCH_REQUIRE_THROWS(testLat.at(Point<double>(0,0)));
    }

    CATCH_SECTION("Square Lattice Constructor - width")
    {
        SquareLattice<TestObject> testLat(1.0);

        CATCH_REQUIRE(testLat.size() == 0);
        CATCH_REQUIRE_THROWS(testLat.randomValue());
        CATCH_REQUIRE_THROWS(testLat.at(Point<double>(0,0)));
    }

    CATCH_SECTION("Square Lattice with one element")
    {
        SquareLattice<TestObject> testLat(1.0);

        // insert object
        testLat.insert(TEST_OBJ(0,0,34));
        CATCH_REQUIRE(testLat.size() == 1);
        CATCH_REQUIRE(testLat.at(Point<double>(0,0)).val == 34);
        CATCH_REQUIRE(testLat.randomValue().val == 34);
    
        // handle collision on insert
        CATCH_REQUIRE_THROWS(testLat.insert(TEST_OBJ(0,0,0)));
        CATCH_REQUIRE_THROWS(testLat.insert(TEST_OBJ(0.5,0,0)));
        CATCH_REQUIRE_THROWS(testLat.insert(TEST_OBJ(-0.5,0,0)));
        CATCH_REQUIRE_THROWS(testLat.insert(TEST_OBJ(0,0.5,0)));
        CATCH_REQUIRE_THROWS(testLat.insert(TEST_OBJ(0,-0.5,0)));

        // move object
        testLat.update(Point<double>(0,0), Point<double>(1,1));
        CATCH_REQUIRE(testLat.size() == 1);
        CATCH_REQUIRE(testLat.at(Point<double>(1,1)).val == 34);
        CATCH_REQUIRE_THROWS(testLat.at(Point<double>(0,0)));

        // erase object
        testLat.erase(Point<double>(1,1));
        CATCH_REQUIRE(testLat.size() == 0);
        CATCH_REQUIRE_THROWS(testLat.at(Point<double>(1,1)));
    }

    CATCH_SECTION("Square Lattice with two elements")
    {
        SquareLattice<TestObject> testLat(1.0);

        testLat.insert(TEST_OBJ(0,0,34));
        testLat.insert(TEST_OBJ(2,2,42));

    }

    CATCH_SECTION("Square Lattice with many elements")
    {


    }


#if 0

    SquareLattice<TestObject> testLat (1.0);

    // test objects
    TestObject obj0 (0.5, 0, 0);
    TestObject obj1 (1.5, 0, 1);
    TestObject obj2 (2.5, 0, 2);
    TestObject objNA (2.0, 0, 2);

    // test insert
    CATCH_REQUIRE_NOTHROW(testLat.insert(obj0.coords, obj0));
    CATCH_REQUIRE_NOTHROW(testLat.insert(obj1.coords, obj1));
    CATCH_REQUIRE_NOTHROW(testLat.insert(obj2.coords, obj2));

    CATCH_REQUIRE(testLat.size() == 3);
    CATCH_REQUIRE_THROWS(testLat.insert(objNA.coords, objNA));
    CATCH_REQUIRE(testLat.size() == 3);

    // test full iterator
    SquareLattice<TestObject>::iterator it1 = testLat.begin();
    int v = 0;
    for (; it1 != testLat.end(); ++it1)
    {
        CATCH_REQUIRE((*it1).val == v++);
    }

    // test local iterator
    SquareLattice<TestObject>::local_iterator it2 = testLat.lbegin(
        Point<double>(0, 0), 1.4);
    v = 0;
    for (; it2 != testLat.lend(Point<double>(0,0), 1.4); ++it2)
    {
        CATCH_REQUIRE((*it2).val == v++);
    }
       
    // test erase
    CATCH_REQUIRE_NOTHROW(testLat.erase(obj1.coords));
    CATCH_REQUIRE(testLat.size() == 2);

    SquareLattice<TestObject>::iterator it3 = testLat.begin();
    CATCH_REQUIRE((*it3).val == 0);
    ++it3;
    CATCH_REQUIRE((*it3).val == 2);

    // test random value - read
    double sum = 0.0;
    for (unsigned i = 0; i < 1000; ++i)
    {
        sum += testLat.randomValue().val;
    }
    CATCH_REQUIRE(sum == 992);

    // test random value - write
    testLat.randomValue().val = 10;
    SquareLattice<TestObject>::iterator it4 = testLat.begin();
    CATCH_REQUIRE((*it4).val == 10);

    // test update
    testLat.update(obj0.coords, Point<double>(10.0,10.0));
    SquareLattice<TestObject>::iterator it5 = testLat.begin();
    CATCH_REQUIRE((*it5).val == 10);
    ++it5;
    CATCH_REQUIRE((*it5).val == 2);
    TestObject objNA2 (10.5, 10.5, 0);
    CATCH_REQUIRE_THROWS(testLat.insert(objNA2.coords, objNA2));
#endif
}


