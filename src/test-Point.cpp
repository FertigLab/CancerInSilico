#include <testthat.h>

#include "Point.h"

#define TEST_APPROX(x) Approx(x).epsilon(0.01)

CATCH_TEST_CASE("Testing Point.h")
{
    CATCH_SECTION("test comparison of points")
    {
        Point p1 = Point(1,1);
        Point p2 = Point(0,2);
        Point p3 = Point(-5,0);
        Point p4 = Point(-5,0);
        Point p5 = Point(0,0);
        Point p6 = Point(0.00005, 0);

        CATCH_REQUIRE(p1 < p2);
        CATCH_REQUIRE(p2 < p3);
        CATCH_REQUIRE(p5 != p1);
        CATCH_REQUIRE(p3 == p4);
        CATCH_REQUIRE(p5 == p6);
    }

    CATCH_SECTION("distance between points")
    {
        Point p1 = Point(0, 0);
        Point p2 = Point(3, 4);
        
        CATCH_REQUIRE(p1.distance(p2) == 5);
    }
}
