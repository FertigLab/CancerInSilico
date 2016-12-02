#include "Simulation.h"
#include "Parameters.h"
#include <Rcpp.h>
#include <cmath>

Simulation::Simulation(Parameters *par) {

    mParams = par;
    mCells = new CellPopulation(mParams, mParams->initialNum(),
                                mParams->density());

}

Simulation::~Simulation() {

	delete mCells;

}

void Simulation::Run() {

	double time = 0.0;
    bool drug_added = false;
    int recordStep = mParams->recordIncrement() / mParams->timeIncrement();
    int outputStep = mParams->outputIncrement() / mParams->timeIncrement();
    int totalSteps = mParams->runTime() / mParams->timeIncrement();

    for (int i = 0; i < totalSteps; i++) {

        if (i % (recordStep ? recordStep : 1) == 0) {

            mCells->RecordPopulation();
        
        }

        if (!drug_added && time > mParams->drugTime()) {
    
            mCells->AddDrug();
            drug_added = true;

        }
        
        Rcpp::checkUserInterrupt();

        if (i % outputStep == 0) {

            Rprintf("time = %.2f\n", ceil(time));
            Rprintf("size = %d\n", mCells->size());

        }

        mCells->OneTimeStep();
		time += mParams->timeIncrement();

    }

    Rprintf("time = %.2f\n", time);
    Rprintf("size = %d\n", mCells->size());

}

Rcpp::List Simulation::GetCellsAsList() {

    return mCells->GetPopulationAsList();

}
