#include <Rcpp.h>

#include "Parameters.h"
#include "Random.h"
#include "CellType.h"
#include "Drug.h"

Parameters::Parameters(Rcpp::S4* rModel)
{
    mRModel = rModel;

    mInitialNum       = rModel->slot("initialNum");
    mRunTime          = rModel->slot("runTime");
    mDensity          = rModel->slot("density");
    mRandSeed         = rModel->slot("randSeed");
    mOutputIncrement  = rModel->slot("outputIncrement");
    mRecordIncrement  = rModel->slot("recordIncrement");
    mTimeIncrement    = rModel->slot("timeIncrement");
    mBoundary         = rModel->slot("boundary");
    mSyncCycles       = rModel->slot("syncCycles");

    Rcpp::List types = rModel->slot("cellTypes");
    for (unsigned i = 0; i < types.size(); ++i)
    {
        mCellTypes.push_back(CellType(i, types[i]));
    }
    
    Rcpp::List drugs = rModel->slot("drugs");
    for (unsigned i = 0; i < drugs.size(); ++i)
    {
        mDrugs.push_back(Drug(i, drugs[i]));
    }
}

CellType Parameters::randomCellType()
{
    Rcpp::NumericVector freq = mRModel->slot("cellTypeInitFreq");

    double total = 0.0, u = Random::uniform(0,1);
    unsigned i = 0;
    while (u >= total) {total += freq[i++];}
    if (i > mCellTypes.size() || i == 0)
        {throw std::runtime_error("cell type out of bounds");}

    return mCellTypes[i-1];
}

void Parameters::setBoundary(double b)
{
    mBoundary = b;
    mRModel->slot("boundary") = b;
}
