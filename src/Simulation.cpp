#include "Simulation.h"
#include "Parameters.h"
#include <Rcpp.h>
#include <cmath>

Simulation::Simulation(Parameters *par, int init_num, double den) {

    m_param = par;
    m_cells = new CellPopulation(m_param, init_num, den);

}

Simulation::~Simulation() {

	delete m_cells;

}

void Simulation::Run(Rcpp::List params) {

	double time = 0.0;
    bool drug_added = false;
    double timeIncrement = params["timeIncrement"];
    int recordStep = (double) params["recordIncrement"] / timeIncrement;
    int outputStep = (double) params["outputIncrement"] / timeIncrement;
    int totalSteps = (double) params["runtime"] / timeIncrement;

    for (int i = 0; i < totalSteps; i++) {

        if (i % (recordStep ? recordStep : 1) == 0) {

            m_cells->RecordPopulation();
        
        }

        if (!drug_added && time > m_param->GetDrugTime()) {
    
            m_cells->AddDrug();
            drug_added = true;

        }
        
        Rcpp::checkUserInterrupt();

        if (i % outputStep == 0) {

            Rprintf("time = %.2f\n", ceil(time));
            Rprintf("size = %d\n", m_cells->size());

        }

        m_cells->OneTimeStep();
		time += timeIncrement;

    }

    Rprintf("time = %.2f\n", time);
    Rprintf("size = %d\n", m_cells->size());

}

Rcpp::List Simulation::GetCellsAsList() {

    return m_cells->GetPopulationAsList();

}
