#include "Parameters.h"
#include "Simulation.h"

#include <Rcpp.h>
#include <iostream>
#include <boost/unordered_map.hpp>

// [[Rcpp::export]]
Rcpp::List CellModel(

    int initialNum,
    int numMCSteps,
    double density,
    double maxTranslation,
    double maxDeform,
    double maxRotate,
    double epsilon,
    double delta,
    int outIncrement,
    int randSeed,
    Rcpp::List drugEffect,
	Rcpp::NumericVector growthRates,
	bool inheritGrowth,
	double nG,
	double timeIncrement,
	int recordIncrement,
    double drugTime,
    double boundary,
    bool syncCellCycle,
    Rcpp::List cellTypes

) {

    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(randSeed);

    Parameters* params = new Parameters(pow(2, 0.5));

    params->SetMaxTranslation(maxTranslation);
    params->SetMaxDeform(maxDeform);
    params->SetMaxRotate(maxRotate);
    params->SetResistanceEPSILON(epsilon);
    params->SetCompressionDELTA(delta);
	params->StoreGrowthDistribution(growthRates);
	params->SetInheritGrowth(inheritGrowth);
	params->SetNG(nG);
    params->StoreDrugEffect(drugEffect);
    params->SetDrugTime(drugTime);
    params->SetBoundary(boundary);
    params->SetSyncCellCycle(syncCellCycle);
    params->StoreCellTypes(cellTypes);

    Simulation main_sim = Simulation(params, initialNum, density);

    main_sim.Run(numMCSteps, outIncrement, timeIncrement, recordIncrement);
    Rcpp::List ret_val = main_sim.GetCellsAsList();
    
    delete params;
    return ret_val;

}




