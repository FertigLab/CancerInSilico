#include <Rcpp.h>
#include <cmath>

#include "CellBasedModel.h"
#include "Parameters.h"
#include "CellPopulation.h"

CellBasedModel::CellBasedModel(Parameters *par)
{
    mParams = par;
    mCells = new CellPopulation(mParams, mParams->initialNum(),
        mParams->density());
}

CellBasedModel::~CellBasedModel()
{
	delete mCells;
}

void CellBasedModel::run()
{
	double time = 0.0;
    double recordTime = 0.0, outputTime = 0.0;

    while (time <= mParams->runTime())
    {
        Rcpp::checkUserInterrupt();

        if (time >= recordTime)
        {
            mCells->RecordPopulation();
            recordTime += mParams->recordIncrement();
        }

        if (time >= outputTime)
        {
            Rprintf("time = %.2f\n", ceil(time));
            Rprintf("size = %d\n", mCells->size());

            outputTime += mParams->outputIncrement();
        }            

        OneTimeStep(time);
		time += mParams->timeIncrement();
    }

    Rprintf("final time = %.2f\n", ceil(time));
    Rprintf("final size = %d\n", mCells->size());
}


