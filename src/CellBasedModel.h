#ifndef CIS_CELL_BASED_MODEL_H
#define CIS_CELL_BASED_MODEL_H

#include <Rcpp.h>

#include "Parameters.h"
#include "CellType.h"
#include "Drug.h"

class CellBasedModel
{
protected:

	std::vector< std::vector<double> > mPopulationRecord;

    Parameters* mParams;
    std::vector<CellType> mCellTypes;
    std::vector<Drug> mDrugs;

public:

    CellBasedModel(Parameters* p) : mParams(p) {}

    /* run the entire model */
    void run();

    /* return model output as R list */
    Rcpp::List getCellsAsList() {return Rcpp::wrap(mPopulationRecord);}

    /* update the model for a single time step; must be implemented */
    virtual void oneTimeStep(double time) = 0;

    /* record the current state of the cell population */    
	virtual void recordPopulation() = 0;
};

#endif
