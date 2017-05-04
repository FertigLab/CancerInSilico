#include <Rcpp.h>
#include <cmath>

#include "CellBasedModel.h"

void CellBasedModel::updateRModel(Rcpp::S4* rModel)
{
    rModel->slot("cells") = Rcpp::wrap(mPopulationRecord);
}

void CellBasedModel::run()
{
	double time = 0.0;
    double recordTime = 0.0, outputTime = 0.0;

    Rprintf("\n");

    while (time <= mParams->runTime())
    {
        Rcpp::checkUserInterrupt();

        if (time >= recordTime)
        {
            recordPopulation();
            recordTime = std::min(recordTime + mParams->recordIncrement(),
                mParams->runTime());
        }

        if (time >= outputTime)
        {
            Rprintf("time = %.2f\n", floor(time));
            Rprintf("size = %d\n", size());

            outputTime = std::min(outputTime + mParams->outputIncrement(),
                mParams->runTime());
        }            

        oneTimeStep(time);
		time += mParams->timeIncrement();
    
        if (time > mParams->runTime() && time < mParams->runTime() +
        mParams->timeIncrement())
        {
            time = mParams->runTime();
        }
    }
}
