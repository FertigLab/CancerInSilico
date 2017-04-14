#ifndef CIS_PARAMETERS_H
#define CIS_PARAMETERS_H

#include <vector>
#include <Rcpp.h>

#include "Drug.h"
#include "CellType.h"

typedef std::vector<Drug>::iterator DrugIterator;

class Parameters
{
protected:

    // general model parameters
    double mInitialNum;
    double mRunTime;
    double mDensity;
    double mRandSeed;
    double mOutputIncrement;
    double mRecordIncrement;
    double mTimeIncrement;
    double mBoundary;
    bool mSyncCycles;

    // pointer to CellModel S4 object from R
    Rcpp::S4* mRModel;

public:

    // constructor
    Parameters(Rcpp::S4*);

    // drugs used in the simulation
    std::vector<Drug> mDrugs;

    // possible cell types
    std::vector<CellType> mCellTypes;

    // get parameter values
    double initialNum()         {return mInitialNum;}
    double runTime()            {return mRunTime;}
    double density()            {return mDensity;}
    double randSeed()           {return mRandSeed;}
    double outputIncrement()    {return mOutputIncrement;}
    double recordIncrement()    {return mRecordIncrement;}
    double timeIncrement()      {return mTimeIncrement;}
    double boundary()           {return mBoundary;}
    bool syncCycles()           {return mSyncCycles;}

    // iterator to drug list
    DrugIterator drugsBegin() {return mDrugs.begin();}
    DrugIterator drugsEnd()   {return mDrugs.end();}

    // set the boundary to a numeric value
    void setBoundary(double);

    // get a random cell type from initial distribution
    CellType randomCellType();
};

#endif
