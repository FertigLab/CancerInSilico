#include "Simulation.h"
#include "Parameters.h"
#include <Rcpp.h>

Simulation::Simulation(Parameters *par) {
    m_param = par;
    m_cells = CellPopulation(m_param, m_param->GetInitialNumCells(), m_param->GetInitialDensity());
}

void Simulation::Run(int dur, int out_incr) {
    m_cells.AddDrug();

    for (int i = 0; i < dur; i++) {
        Rcpp::checkUserInterrupt();

        if (i % out_incr == 0) {
            Rprintf("time = %d\n", i);
            Rprintf("size = %d\n", m_cells.size());
        }

        m_cells.OneTimeStep();
    }

    Rprintf("time = %d\n", dur);
    Rprintf("size = %d\n", m_cells.size());
}

Rcpp::NumericMatrix Simulation::GetCellsAsMatrix() {
    return m_cells.GetPopulationAsMatrix();
}
