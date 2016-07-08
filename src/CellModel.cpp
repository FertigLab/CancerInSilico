#include "Parameters.h"
#include "Simulation.h"

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix CellModel(

    int initialNum,
    int runTime,
    double density,
    double meanGrowth,
    double varGrowth,
    double maxMigration,
    double maxDeform,
    double maxRotate,
    double epsilon,
    double delta,
    int outIncrement,
    int randSeed

) {

    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(randSeed);

    Parameters *params = new Parameters(pow(2,0.5));

    double apoptosisRate = 0.0;

    params->SetEnergyConstant(1);
    params->SetInitialNumCells(initialNum);
    params->SetInitialDensity(density);
    params->SetMeanGrowth(meanGrowth);
    params->SetVarGrowth(varGrowth);
    params->SetApoptosisRate(apoptosisRate);
    params->SetMaxMigration(maxMigration);
    params->SetMaxDeform(maxDeform);
    params->SetMaxRotate(maxRotate);
    params->SetResistanceEPSILON(epsilon);
    params->SetCompressionDELTA(delta);

    Simulation main_sim = Simulation(params);
    main_sim.Run(runTime, outIncrement);
    Rcpp::NumericMatrix ret_val = main_sim.GetCellsAsMatrix();

    delete params;
    return ret_val;

}




