#include "Parameters.hpp"
#include "Simulation.hpp"

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
  double delta

) {

  Parameters* params = new Parameters();

  double apoptosisRate = 0.0;

	params->SetMinRadius(1);
	params->SetMaxRadius(pow(2,0.5) * params->GetMinRadius());
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

  main_sim.Run(runTime);

  Rcpp::NumericMatrix ret_val = main_sim.GetCellsAsMatrix();

  delete params;

  return ret_val;

}

