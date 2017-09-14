#include "../catch.h"

#include "../../Core/Point.h"

TEST_CASE("Testing Point.h")
{
    SECTION("test comparison of points")
    {
        Point<double> p1(1,1);
        Point<double> p2(0, 2);
        Point<double> p3(-5, 0);
        Point<double> p4(-5, 0);
        Point<double> p5(0, 0);
        Point<double> p6(0.00005, 0);

        REQUIRE(p1 < p2);
        REQUIRE(p2 < p3);
        REQUIRE(p5 != p1);
        REQUIRE(p3 == p4);
        REQUIRE(p5 != p6);
    }

    SECTION("distance between points")
    {
        Point<double> p1(0,0);
        Point<double> p2(3,4);
        
        REQUIRE(p1.distance(p2) == 5);
    }
}
