#include <Rcpp.h>

#include "Core/Parameters.h"
#include "Core/CellBasedModel.h"
#include "CellModels/DrasdoHohmeModel.h"

void createModel(Rcpp::List rParams, Parameters*& params,
CellBasedModel*& model)
{
    params = new DrasdoHohmeParameters(rParams);
    model = new DrasdoHohmeModel(static_cast<DrasdoHohmeParameters*>
        (params));
    Random::setSeed(params->randSeed());
}

// [[Rcpp::export]]
Rcpp::List cppRunModel(Rcpp::List rParams)
{
    Parameters* modelParams; 
    CellBasedModel* model;

    createModel(rParams, modelParams, model);
    model->run();

    Rcpp::List modelOutput;
    modelOutput["cells"] = model->getCellsAsList();
    modelOutput["params"] = modelParams->getRParameters();
    
    delete modelParams;
    delete model;
    return modelOutput;
}




