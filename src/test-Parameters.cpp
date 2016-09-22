#include <iostream>
#include <cmath>

#include "Parameters.h"

CATCH_TEST_CASE("Test Parameters and Radius Solver") {
/*
    CATCH_SECTION("Test max error of Radius Solver") {

        Parameters params = Parameters(pow(2,0.5));

        double radius, axis_len, theta, error, max_error = 0.0;
        for (int i = 282843; i <= 400000; ++i) {

            axis_len = (double) i / 100000;
            theta = params.GetThetaSlow(axis_len);
            radius = params.GetRadius(axis_len);
            error = fabs(radius - pow(2 * M_PI / (2 * M_PI - theta + sin(theta)), 0.5));
            if (error > max_error) { max_error = error;}
            error = fabs(radius - axis_len / (2 + 2 * cos(theta / 2)));
            if (error > max_error) { max_error = error;}

        }

        CATCH_REQUIRE(max_error < 0.001);

    }
*/
}
