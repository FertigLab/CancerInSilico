#include <Rcpp.h>
#include <cmath>

#include "CellBasedModel.h"
#include "Random.h"

CellBasedModel::CellBasedModel(Rcpp::S4* rModel)
{
    // store basic parameters
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
        mDrugs.push_back(drugs[i]);
    }
}

// set actual boundary value
void CellBasedModel::setBoundary(double b)
{
    mBoundary = b;
    mRModel->slot("boundary") = b;
}

// get random cell type based on given frequency
CellType CellBasedModel::randomCellType()
{
    // get distribution of cell types
    Rcpp::NumericVector freq = mRModel->slot("cellTypeInitFreq");
    double total = 0.0, u = Random::uniform(0,1);
    unsigned i = 0;

    // increment until correct bin is found
    while (u >= total && total < 1) {total += freq[i++];}

    // return cell type in this random bin
    return mCellTypes[i-1];
}

// run the entire model
void CellBasedModel::run()
{
    double time = 0.0;
    double recordTime = 0.0, outputTime = 0.0;

    Rprintf("\n");

    while (time <= runTime())
    {
        Rcpp::checkUserInterrupt();

        if (time >= recordTime) // record cell state
        {
            recordPopulation();
            recordTime = std::min(recordTime + recordIncrement(),runTime());
        }

        if (time >= outputTime) // output time and num cells
        {
            Rprintf("time = %.2f\n", floor(time));
            Rprintf("size = %d\n", size());

            outputTime = std::min(outputTime + outputIncrement(),runTime());
        }            

        oneTimeStep(time); // run the simulation for one time step
		time += timeIncrement();
    
        // ensures last time step happens at end of given runTime
        if (time > runTime() && time < runTime() + timeIncrement())
        {
            time = runTime();
        }
    }
}
