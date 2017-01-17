#include <iostream>
#include <cmath>
#include <testthat.h>
#include <Rcpp.h>
#include "Parameters.h"

#define TEST_APPROX(x) Approx(x).epsilon(0.01)

CATCH_TEST_CASE("parameters class") {

    Rcpp::Environment env;
    env = Rcpp::Environment::namespace_env("CancerInSilico");

    Rcpp::List Rparams = env.find("testParams");

    CATCH_REQUIRE(Rparams.size() == 16);

    Parameters* par;

    CATCH_REQUIRE_NOTHROW(par = new Parameters(pow(2, 0.5), Rparams));

    CATCH_SECTION("parameters constructor") {
    
        CATCH_REQUIRE(par->initialNum() == 10);
        CATCH_REQUIRE(par->runTime() == 24);
        CATCH_REQUIRE(par->density() == 0.01);
        CATCH_REQUIRE(!par->inheritGrowth());
        CATCH_REQUIRE(par->drugTime() == 0.0);
        CATCH_REQUIRE(par->boundary());
        CATCH_REQUIRE(par->randSeed() == 0);
        CATCH_REQUIRE(par->syncCycles());
        CATCH_REQUIRE(par->outputIncrement() == 6);
        CATCH_REQUIRE(par->recordIncrement() == 0.25);
        CATCH_REQUIRE(par->nG() == 24);
        CATCH_REQUIRE(par->epsilon() == 10);
        CATCH_REQUIRE(par->delta() == 0.2);
        CATCH_REQUIRE(par->maxRadius() == pow(2, 0.5));
        CATCH_REQUIRE(par->maxDeform() == 0.1);
        CATCH_REQUIRE(par->maxTranslation() == 0.1);
        CATCH_REQUIRE(par->maxRotate() == TEST_APPROX(0.3095));
        CATCH_REQUIRE(par->timeIncrement() == TEST_APPROX(0.008));

        CATCH_REQUIRE(par->GetRandomGrowthRate('A') == TEST_APPROX(0.0003));
        CATCH_REQUIRE(par->GetDrugEffect(0.01) == TEST_APPROX(5.6019));

    }

    CATCH_SECTION("max error of radius solver") {

        double radius, axis_len, theta, error, max_error = 0.0;

        for (int i = 282843; i <= 400000; ++i) {

            axis_len = (double) i / 100000;
            theta = par->GetThetaSlow(axis_len);
            radius = par->GetRadius(axis_len);

            error = fabs(radius - pow(2 * M_PI / (2 * M_PI - theta
                         + sin(theta)), 0.5));
            if (error > max_error) { max_error = error;}

            error = fabs(radius - axis_len / (2 + 2 * cos(theta / 2)));
            if (error > max_error) { max_error = error;}

        }

        CATCH_REQUIRE(max_error < 0.001);

    }

}
