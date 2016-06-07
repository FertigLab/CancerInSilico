#include "Parameters.hpp"
#include "Simulation.hpp"

#include <Rcpp.h>


// [[Rcpp::export]]
Rcpp::NumericMatrix CellModel(

  int initialNum,
  int runTime

) {

  double density = 0.05;
  double meanGrowth = 0.15;
  double varGrowth = 0.0;
  double apoptosisRate = 0.0;
  double maxMigration = 0.5;
  double maxDeform = 0.075;
  double maxRotate = 0.3;
  double epsilon = 0.05;
  double delta = 5.0;

  Parameters* params = new Parameters();

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

