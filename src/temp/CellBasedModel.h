#ifndef CIS_CELL_BASED_MODEL_H
#define CIS_CELL_BASED_MODEL_H

#include <Rcpp.h>

#include "CellPopulation.h"
#include "Parameters.h"

class CellBasedModel
{
protected:

	std::vector< std::vector<double> > mPopulationRecord;

    CellPopulation* mCells;
    Parameters* mParams;

public:

    CellBasedModel(Parameters*);
    ~CellBasedModel();

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
