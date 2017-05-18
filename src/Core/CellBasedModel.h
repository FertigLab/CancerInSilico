#ifndef CIS_CELL_BASED_MODEL_H
#define CIS_CELL_BASED_MODEL_H

// top-level abstract class for a cell model

#include <Rcpp.h>

#include "CellType.h"
#include "Drug.h"

class CellBasedModel
{
protected:

    // record of the cell states at specified times during the simulation
	std::vector< std::vector<double> > mPopulationRecord;

    // S4 class in R describing cell model
    Rcpp::S4* mRModel;

    // model parameters
    double mInitialNum;
    double mRunTime;
    double mDensity;
    double mRandSeed;
    double mOutputIncrement;
    double mRecordIncrement;
    double mTimeIncrement;
    double mBoundary;
    bool mSyncCycles;

    // drugs used in the simulation
    std::vector<Drug> mDrugs;

    // possible cell types
    std::vector<CellType> mCellTypes;

public:

    // constructor and virtual destructor
    CellBasedModel(Rcpp::S4*);
    virtual ~CellBasedModel() {}    

    // get parameter values
    double initialNum()         const {return mInitialNum;}
    double runTime()            const {return mRunTime;}
    double density()            const {return mDensity;}
    double randSeed()           const {return mRandSeed;}
    double outputIncrement()    const {return mOutputIncrement;}
    double recordIncrement()    const {return mRecordIncrement;}
    double timeIncrement()      const {return mTimeIncrement;}
    double boundary()           const {return mBoundary;}
    bool syncCycles()           const {return mSyncCycles;}

    // set actual boundary value (passed in as T/F)
    void setBoundary(double);

    // get a random cell type from initial distribution
    CellType randomCellType();

    // run the entire model
    void run();

    // update the R model
    Rcpp::List getCellRecord() {return Rcpp::wrap(mPopulationRecord);}

    // update the model for a single time step; must be implemented
    virtual void oneTimeStep(double time) = 0;

    // record the current state of the cell population
	virtual void recordPopulation() = 0;

    // number of cells in the model 
    virtual unsigned size() const = 0;
};

#endif
