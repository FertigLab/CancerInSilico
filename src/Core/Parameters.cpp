#include <Rcpp.h>

#include "Parameters.h"
#include "Random.h"
#include "CellType.h"
#include "Drug.h"

// construct parameters from R model
Parameters::Parameters(Rcpp::S4* rModel)
{
    mRModel = rModel;

    // basic parameters
    mInitialNum       = rModel->slot("initialNum");
    mRunTime          = rModel->slot("runTime");
    mDensity          = rModel->slot("density");
    mRandSeed         = rModel->slot("randSeed");
    mOutputIncrement  = rModel->slot("outputIncrement");
    mRecordIncrement  = rModel->slot("recordIncrement");
    mTimeIncrement    = rModel->slot("timeIncrement");
    mBoundary         = rModel->slot("boundary");
    mSyncCycles       = rModel->slot("syncCycles");

    // create cell types
    Rcpp::List types = rModel->slot("cellTypes");
    for (unsigned i = 0; i < types.size(); ++i)
    {
        mCellTypes.push_back(CellType(i, types[i]));
    }
    
    // create drugs
    Rcpp::List drugs = rModel->slot("drugs");
    for (unsigned i = 0; i < drugs.size(); ++i)
    {
        mDrugs.push_back(Drug(i, drugs[i]));
    }
}

// get random cell type based on given frequency
CellType Parameters::randomCellType()
{
    Rcpp::NumericVector freq = mRModel->slot("cellTypeInitFreq");

    double total = 0.0, u = Random::uniform(0,1);
    unsigned i = 0;

    // increment until correct bin is found
    while (u >= total) {total += freq[i++];}
    if (i > mCellTypes.size() || i == 0)
        {throw std::runtime_error("cell type out of bounds");}

    // return cell type in this random bin
    return mCellTypes[i-1];
}

// set actual boundary value
void Parameters::setBoundary(double b)
{
    mBoundary = b;
    mRModel->slot("boundary") = b;
}
