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

CATCH_TEST_CASE("Test SquareLattice.h")
{
    Random::setSeed(0);

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
//    CATCH_REQUIRE(v == 2);
       
    // test erase
    CATCH_REQUIRE_NOTHROW(testLat.erase(obj1.coords));
    CATCH_REQUIRE(testLat.size() == 2);

    SquareLattice<TestObject>::iterator it3 = testLat.begin();
    CATCH_REQUIRE((*it3).val == 0);
    ++it3;
    CATCH_REQUIRE((*it3).val == 2);

    // test random value
    CATCH_REQUIRE(testLat.randomValue().val == 2);
    testLat.randomValue().val = 10;
    
    SquareLattice<TestObject>::iterator it4 = testLat.begin();
    ++it4;
    CATCH_REQUIRE((*it4).val == 10);

    // test update

}


