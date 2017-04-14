#include <Rcpp.h>

#include "Core/Parameters.h"
#include "Core/CellBasedModel.h"
#include "CellModels/DrasdoHohmeModel.h"

void createModel(Rcpp::S4* rModel, CellBasedModel*& cModel,
std::string type)
{
    if (!type.compare("DrasdoHohme"))
    {
        cModel = new DrasdoHohmeModel(rModel);
    }
    else
    {
        throw std::invalid_argument("invalid model type");
    }
}

// [[Rcpp::export]]
Rcpp::S4 cppRunModel(Rcpp::S4 rModel, std::string type)
{
    Random::setSeed(rModel.slot("randSeed"));
    
    CellBasedModel* cModel;
    createModel(&rModel, cModel, type);
    cModel->run();

    rModel.slot("cells") = cModel->getCellsAsList();

    delete cModel;
    return rModel;
}




