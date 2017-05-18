#include <Rcpp.h>

#include "Core/CellBasedModel.h"
#include "CellModels/DrasdoHohmeModel.h"

// figure out which model implementation to use
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
    
    // create C model, run, and store in R model
    CellBasedModel* cModel;
    createModel(&rModel, cModel, type);
    
    try
    {
        cModel->run();
    }
    catch (const std::exception& e)
    {
        Rcpp::Rcout << e.what() << std::endl;
    }        
    rModel.slot("cells") = cModel->getCellRecord();

    delete cModel;
    return rModel;
}




