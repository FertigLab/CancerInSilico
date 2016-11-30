#include "Parameters.h"
#include "Simulation.h"

#include <Rcpp.h>
#include <iostream>
#include <boost/unordered_map.hpp>

// [[Rcpp::export]]
Rcpp::List CellModel(

    Rcpp::List params

) {

    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(params["randSeed"]);

    Parameters* parameters = new Parameters(pow(2, 0.5));

    parameters->SetMaxTranslation(params["maxTranslation"]);
    parameters->SetMaxDeform(params["maxDeform"]);
    parameters->SetMaxRotate(params["maxRotate"]);
    parameters->SetResistanceEPSILON(params["epsilon"]);
    parameters->SetCompressionDELTA(params["delta"]);
	parameters->StoreGrowthDistribution(params["maxGrowth"]);
	parameters->SetInheritGrowth(params["inheritGrowth"]);
	parameters->SetNG(params["nG"]);
    parameters->StoreDrugEffect(params["drugEffect"]);
    parameters->SetDrugTime(params["drugTime"]);
    parameters->SetBoundary(params["boundary"]);
    parameters->SetSyncCellCycle(params["syncCellCycle"]);

    Simulation main_sim = Simulation(parameters, params["initialNum"],
        params["density"]);

    main_sim.Run(params);
    Rcpp::List ret_val = main_sim.GetCellsAsList();
    
    delete parameters;
    return ret_val;

}




