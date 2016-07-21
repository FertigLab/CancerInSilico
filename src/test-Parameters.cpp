#include <iostream>
#include <cmath>

#include "Parameters.h"

CATCH_TEST_CASE("Test Parameters and Radius Solver") {

    CATCH_SECTION("Test max error of Radius Solver") {

        Parameters params = Parameters(pow(2,0.5));

        double radius, theta, error, max_error = 0.0;
        for (int i = 100000; i <= 141421; ++i) {

            radius = (double) i / 100000;
            theta = params.GetTheta(radius);
            error = fabs(radius - params.GetMaxRadius() * pow(M_PI / (2 * M_PI - theta + sin(theta)), 0.5));
            if (error > max_error) { max_error = error;}

        }

        CATCH_REQUIRE(max_error < 0.001);

    }

}
