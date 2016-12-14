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
    double recordTime = 0.0, outputTime = 0.0;

    while (time <= mParams->runTime()) {

        Rcpp::checkUserInterrupt();

        if (time >= recordTime) {

            mCells->RecordPopulation();
            recordTime += mParams->recordIncrement();
        
        }

        if (time >= outputTime) {

            Rprintf("time = %.2f\n", outputTime);
            Rprintf("size = %d\n", mCells->size());

            outputTime += mParams->outputIncrement();

        }            

        if (!drug_added && time > mParams->drugTime()) {
    
            mCells->AddDrug();
            drug_added = true;

        }
        
        mCells->OneTimeStep();
		time += mParams->timeIncrement();

        if (time >= mParams->runTime()
                && time < mParams->runTime() + mParams->timeIncrement()) {

            time = mParams->runTime();
            outputTime = time;
            recordTime = time;

        }

    }

}

Rcpp::List Simulation::GetCellsAsList() {

    return mCells->GetPopulationAsList();

}
