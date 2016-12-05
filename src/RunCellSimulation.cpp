#include "Parameters.h"
#include "Simulation.h"

#include <Rcpp.h>
#include <iostream>
#include <boost/unordered_map.hpp>

// [[Rcpp::export]]
Rcpp::List runCellSimulation(

    Rcpp::List Rparams

) {

    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(Rparams["randSeed"]);

    Parameters* Cparams = new Parameters(pow(2, 0.5), Rparams);

    Simulation mainSim = Simulation(Cparams);

    mainSim.Run();
    Rcpp::List cellData = mainSim.GetCellsAsList();
    
    delete Cparams;
    return cellData;

}




