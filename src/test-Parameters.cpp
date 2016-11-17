#include <iostream>
#include <cmath>
#include <testthat.h>

#include "Parameters.h"

CATCH_TEST_CASE("Test Parameters with Cell Type") {

    /* Create list of CellType objects */
    Rcpp::List cell_types;

    /* Create 2 CellType objects */
    Rcpp::S4 c1 = Rcpp::S4("CellType");
    Rcpp::S4 c2 = Rcpp::S4("CellType");
    Rcpp::S4 c3 = Rcpp::S4("CellType");
    
    /* Name celltypes, set slots */
    c1.slot("mType") = "NORMAL";
    c2.slot("mType") = "CANCER";
    c3.slot("mType") = "RESISTANT";

    /* Add to list */
    cell_types.push_back(c1);
    cell_types.push_back(c2);
    cell_types.push_back(c3);

    // Create NumericVector of Cell Type Distributions
    Rcpp::NumericVector cell_type_dist
            = Rcpp::NumericVector::create(0.25,0.25,0.5);

    // Set a growth rates
    Rcpp::NumericVector growth_dist = Rcpp::NumericVector::create(1,2);
    
    /* Create DrugEffect List */
    Rcpp::List drug_effect;

    /* add distributions for the two growth rates */
    drug_effect.push_back(Rcpp::NumericVector::create(1, 0.2));
    drug_effect.push_back(Rcpp::NumericVector::create(2, 0.3));

    /* Create a Parameters object */
    Parameters params = Parameters(pow(2,0.5));
    params.StoreCellTypes(cell_types);
    params.StoreCellTypeDistribution(cell_type_dist);
    params.StoreGrowthDistribution(growth_dist);
    params.StoreDrugEffect(drug_effect);

    // tests GetThetaSlow() and GetRadius()
    CATCH_SECTION("Test max error of Radius Solver") {

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

    CATCH_SECTION("Test GetRandomGrowthRate") {

        Rcpp::Environment baseEnv("package:base");
        Rcpp::Function setSeed = baseEnv["set.seed"];
        setSeed(25);

        CATCH_REQUIRE(params.GetRandomGrowthRate() == 1);
        CATCH_REQUIRE(params.GetRandomGrowthRate() == 2);

    }

    CATCH_SECTION("Test GetRandomCellType") {

        Rcpp::Environment baseEnv("package:base");
        Rcpp::Function setSeed = baseEnv["set.seed"];
        setSeed(20);

        // std::vector<int> sellected_cell_types(3);
        int sum_cell_types = 0;

        for (int i = 0 ; i < 1000; i++) {

            //sellected_cell_types[params.GetRandomCellType()]++;
            sum_cell_types += params.GetRandomCellType();

        }

        //std::cout << sum_cell_types;
        CATCH_REQUIRE(sum_cell_types++ == 1259);
    }

    CATCH_SECTION("Test StoreCellTypes") {

        std::vector<Rcpp::S4> list_of_cell_types = params.GetCellTypes();

        CATCH_REQUIRE(list_of_cell_types.size() == 3);
        CATCH_REQUIRE(list_of_cell_types[0].hasSlot("mType"));
        CATCH_REQUIRE(list_of_cell_types[1].hasSlot("mType"));
        CATCH_REQUIRE(list_of_cell_types[2].hasSlot("mType"));

    }


}