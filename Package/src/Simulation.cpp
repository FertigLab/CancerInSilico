#include "Simulation.hpp"
#include "Parameters.hpp"
#include <Rcpp.h>

Simulation::Simulation(Parameters* par) {

  m_param = par;  
  m_cells = CellPopulation(m_param, m_param->GetInitialNumCells(), m_param->GetInitialDensity());

}

void Simulation::Run(int dur) {

  m_cells.AddDrug();

  for (int i = 0; i < dur; i++) {
    
    Rcpp::checkUserInterrupt();

    m_cells.OneTimeStep();
    std::cout << i << std::endl;

  }  

}

Rcpp::NumericMatrix Simulation::GetCellsAsMatrix() {

  return m_cells.GetPopulationAsMatrix();

}
