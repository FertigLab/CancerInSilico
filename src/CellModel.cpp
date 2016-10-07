#include "Parameters.h"
#include "Simulation.h"

#include <Rcpp.h>

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
    Rcpp::NumericVector growthRates,
    bool inheritGrowth,
    double nG,
    double timeIncrement

) {

    std::vector<double> gr_rates;
    for (unsigned int i = 0; i < growthRates.length(); ++i) {
    
        gr_rates.push_back(growthRates[i]);

    }        

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

    Simulation main_sim = Simulation(params, initialNum, density);
    main_sim.Run(numMCSteps, outIncrement, timeIncrement);
    Rcpp::List ret_val = main_sim.GetCellsAsList();

    delete params;
    return ret_val;

}




