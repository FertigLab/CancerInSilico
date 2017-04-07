#include <Rcpp.h>

#include "Parameters.h"
#include "Random.h"
#include "CellType.h"
#include "Drug.h"

Parameters::Parameters(Rcpp::List rParams)
{
    mParams = rParams;

    Rcpp::List types = mParams["cellTypes"];
    for (unsigned i = 0; i < types.size(); ++i)
    {
        mCellTypes.push_back(CellType(i, types[i]));
    }
    
    Rcpp::List drugs = mParams["drugs"];
    for (unsigned i = 0; i < drugs.size(); ++i)
    {
        mDrugs.push_back(Drug(i, drugs[i]));
    }
}

const CellType& Parameters::randomCellType()
{
    Rcpp::NumericVector freq = mParams["cellTypeInitFreq"];

    double total = 0.0, u = Random::uniform(0,1);
    unsigned i = 0;
    while (u > total) {total += freq[i++];}

    return mCellTypes[i-1];
}
