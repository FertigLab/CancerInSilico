#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "CellPopulation.h"
#include "Parameters.h"

#include <Rcpp.h>

class Simulation {

  private:

    CellPopulation m_cells;
    Parameters *m_param;

  public:

    Simulation(Parameters *);
    ~Simulation() {}

    void Run(int, int);
    Rcpp::NumericMatrix GetCellsAsMatrix();

};

#endif
