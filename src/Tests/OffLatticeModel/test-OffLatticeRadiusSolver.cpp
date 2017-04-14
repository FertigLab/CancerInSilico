#include <Rcpp.h>

#include "../TestHeader.h"
#include "../../OffLatticeModel/OffLatticeRadiusSolver.h"

CATCH_TEST_CASE("Test OffLatticeRadiusSolver") {

    OffLatticeRadiusSolver solver;
    double radius, theta, error, max_error = 0.0;

    for (double axis = 2 * sqrt(2); axis <= 4.0; axis += 0.00001)
    {
        radius = solver.radius(axis);
        theta = 2 * acos(axis / (2 * radius) - 1);

        error = fabs(axis - 2 * radius - 2 * radius * cos(theta / 2));
        if (error > max_error) { max_error = error;}
    }
    CATCH_REQUIRE(max_error < 0.001);
}
