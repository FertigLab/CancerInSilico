#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "CellPopulation.h"
#include "Parameters.h"

#include <Rcpp.h>

class Simulation {

  private:

    CellPopulation* mCells;
    Parameters* mParams;

  public:

    Simulation(Parameters*);
    ~Simulation();

    void Run();
    Rcpp::List GetCellsAsList();

};

#endif
