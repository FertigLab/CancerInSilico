#include "Simulation.h"
#include "Parameters.h"
#include <Rcpp.h>

Simulation::Simulation(Parameters *par, int init_num, double den) {

    m_param = par;
    m_cells = new CellPopulation(m_param, init_num, den);

}

Simulation::~Simulation() {

	delete m_cells;

}

void Simulation::Run(int MCsteps, int out_incr, double time_incr, double rec_incr) {

	double time = 0.0;
    bool drug_added = false;
	m_cells->RecordPopulation();
	
    for (int i = 0; i < MCsteps; i++) {

        if (!drug_added && time > 5.0) {
    
            m_cells->AddDrug();
            drug_added = true;

        }
        
        Rcpp::checkUserInterrupt();

        if (i % out_incr == 0) {

            Rprintf("time = %.2f\n", time);
            Rprintf("size = %d\n", m_cells->size());

        }

        m_cells->OneTimeStep();
		time += time_incr;
		
        m_cells->RecordPopulation(); //TODO: don't record every MC step

    }

    Rprintf("time = %.2f\n", time);
    Rprintf("size = %d\n", m_cells->size());

}

Rcpp::List Simulation::GetCellsAsList() {

    return m_cells->GetPopulationAsList();

}
