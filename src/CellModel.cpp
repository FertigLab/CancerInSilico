#include "Parameters.h"
#include "Simulation.h"

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix CellModel(

    int initialNum,
    int runTime,
    double density,
    double maxMigration,
    double maxDeform,
    double maxRotate,
    double epsilon,
    double delta,
    int outIncrement,
    int randSeed,
	Rcpp::NumericVector growthRates,
	bool inheritGrowth

) {

	std::vector<double> gr_rates;
	for (unsigned int i = 0; i < growthRates.length(); ++i) {
	
		gr_rates.push_back(growthRates[i]);

	}		

    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(randSeed);

    Parameters *params = new Parameters(pow(2,0.5));

    double apoptosisRate = 0.0;

    params->SetEnergyConstant(1);
    params->SetInitialNumCells(initialNum);
    params->SetInitialDensity(density);
    params->SetApoptosisRate(apoptosisRate);
    params->SetMaxMigration(maxMigration);
    params->SetMaxDeform(maxDeform);
    params->SetMaxRotate(maxRotate);
    params->SetResistanceEPSILON(epsilon);
    params->SetCompressionDELTA(delta);
	params->SetInheritGrowth(inheritGrowth);

    Simulation main_sim = Simulation(params);
    main_sim.Run(runTime, outIncrement, gr_rates);
    Rcpp::NumericMatrix ret_val = main_sim.GetCellsAsMatrix();

    delete params;
    return ret_val;

}




