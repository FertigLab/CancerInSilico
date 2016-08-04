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
	double timeIncrement

) {

	std::vector<double> gr_rates;
	for (unsigned int i = 0; i < growthRates.length(); ++i) {
	
		gr_rates.push_back(growthRates[i]);

	}	

    boost::unordered_map<double, std::vector<double> > drug_effect;	
    for (unsigned int i = 0; i < drugEffect.size(); ++i) {

        std::vector<double> dist = Rcpp::as< std::vector<double> >(drugEffect[i]);   

        double cycle_time = dist[0];
        dist[0] = dist.back();
        dist.pop_back();     

        drug_effect.insert(std::pair<double, std::vector<double> >(cycle_time, dist));

    }

    /*boost::unordered_map<double, std::vector<double> >::iterator iter;
    for (iter = drug_effect.begin(); iter != drug_effect.end(); ++iter) {

        std::cout << iter->first << std::endl;
        for (unsigned int i = 0; i < iter->second.size(); ++i) {

            std::cout << "\t" << iter->second[i] << std::endl;

        }

    }*/

    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(randSeed);

    Parameters *params = new Parameters(pow(2,0.5));

    params->SetMaxTranslation(maxTranslation);
    params->SetMaxDeform(maxDeform);
    params->SetMaxRotate(maxRotate);
    params->SetResistanceEPSILON(epsilon);
    params->SetCompressionDELTA(delta);
	params->StoreGrowthDistribution(gr_rates);
	params->SetInheritGrowth(inheritGrowth);
	params->SetNG(nG);
    params->StoreDrugEffect(drug_effect);

    Simulation main_sim = Simulation(params, initialNum, density);
    main_sim.Run(numMCSteps, outIncrement, timeIncrement);
    Rcpp::List ret_val = main_sim.GetCellsAsList();

    delete params;
    return ret_val;

}




