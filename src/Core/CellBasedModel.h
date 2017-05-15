#ifndef CIS_CELL_BASED_MODEL_H
#define CIS_CELL_BASED_MODEL_H

// top-level abstract class for a cell model

#include <Rcpp.h>

#include "Parameters.h"
#include "CellType.h"
#include "Drug.h"

class CellBasedModel
{
protected:

    // record of the cell states at specified times during the simulation
	std::vector< std::vector<double> > mPopulationRecord;

    // parameters object for this model
    Parameters* mParams;

public:

    // constructors
    CellBasedModel(Rcpp::S4* rM) {mParams = new Parameters(rM);}
    virtual ~CellBasedModel() {delete mParams;}

    // run the entire model
    void run();

    // update the R model
    virtual void updateRModel(Rcpp::S4*);

    // update the model for a single time step; must be implemented
    virtual void oneTimeStep(double time) = 0;

    // record the current state of the cell population
	virtual void recordPopulation() = 0;

    // number of cells in the model 
    virtual unsigned size() const = 0;
};

#endif
