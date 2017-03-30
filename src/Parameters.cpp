#include <Rcpp.h>

#include "Parameters.h"

Parameters::Parameters(Rcpp::List Rparams)
{
    mParams = Rparams;

    Rcpp::List types = mParams["cellTypes"];
    for (unsigned i = 0; i < types.size(); ++i)
    {
        mCellTypes.push_back(types[i]);
    }
}

Rcpp::S4* Parameters::randomCellType()
{
    Rcpp::NumericVector freq = mParams["cellTypeInitFreq"];

    double total = 0.0, u = Random::uniform(0,1);
    unsigned i = 0;
    while (u > total) { total += freq[i++];}

    return &mCellTypes[i-1];
}
