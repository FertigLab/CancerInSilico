#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "CellPopulation.h"
#include "Parameters.h"

#include <Rcpp.h>

class CellBasedModel {

private:

    CellPopulation* mCells;
    Parameters* mParams;

public:

    CellBasedModel(Parameters*);
    ~CellBasedModel();

    void Run();
    virtual void OneTimeStep() = 0;
    Rcpp::List GetCellsAsList();

};

#endif
