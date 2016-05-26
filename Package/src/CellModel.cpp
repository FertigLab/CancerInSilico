#include "Parameters.hpp"
#include "Simulation.hpp"

#include <Rcpp.h>


// [[Rcpp::export]]
Rcpp::NumericMatrix CellModel(SEXP init_num, SEXP run_time, SEXP density) {

  Parameters* params = new Parameters();

  params->SetInitialNumCells(Rcpp::as<int>(init_num));
	params->SetMinRadius(1);
	params->SetMaxRadius(pow(2,0.5) * params->GetMinRadius());
	params->SetMeanGrowth(0.15);
	params->SetVarGrowth(0.05);
  params->SetDrugProportion(0.5);
	params->SetMaxMigration(0.5);
	params->SetMaxRotate(0.3);
	params->SetMaxDeform(0.075);
  params->SetApoptosisRate(0.0);
	params->SetEnergyConstant(1);
	params->SetResistanceEPSILON(0.05);
	params->SetCompressionDELTA(5);
  params->SetInitialDensity(Rcpp::as<float>(density));

  Simulation main_sim = Simulation(params);

  main_sim.Run(Rcpp::as<int>(run_time));

  Rcpp::NumericMatrix ret_val = main_sim.GetCellsAsMatrix();

  delete params;
  
  return ret_val;

}

