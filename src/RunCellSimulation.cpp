#include <Rcpp.h>

#include "Parameters.h"
#include "CellBasedModel.h"

void createModel(Rcpp::List rParams, Parameters*& params,
CellBasedModel*& model)
{
    params = new ContinuumParameters(rParams);
    model = new ContinuumModel(params);
}

// [[Rcpp::export]]
Rcpp::List runCellSimulation(Rcpp::List rParams)
{
    Random::setSeed(rParams["randSeed"]);

    Parameters* modelParams; 
    CellBasedModel* model;

    createModel(rParams, modelParams, model);

    model->run();

    Rcpp::List modelOutput;
    modelOutput["cells"] = model->getCellsAsList();
    modelOutput["params"] = modelParams->getRparameters();
    
    delete modelParams;
    delete model;
    return modelOutput;
}




