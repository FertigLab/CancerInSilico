#include "../TestHeader.h"
#include "../../Core/Point.h"

CATCH_TEST_CASE("Testing Point.h")
{
    CATCH_SECTION("test comparison of points")
    {
        Point<double> p1(1,1);
        Point<double> p2(0, 2);
        Point<double> p3(-5, 0);
        Point<double> p4(-5, 0);
        Point<double> p5(0, 0);
        Point<double> p6(0.00005, 0);

        CATCH_REQUIRE(p1 < p2);
        CATCH_REQUIRE(p2 < p3);
        CATCH_REQUIRE(p5 != p1);
        CATCH_REQUIRE(p3 == p4);
        CATCH_REQUIRE(p5 != p6);
    }

    CATCH_SECTION("distance between points")
    {
        Point<double> p1(0,0);
        Point<double> p2(3,4);
        
        CATCH_REQUIRE(p1.distance(p2) == 5);
    }
}
