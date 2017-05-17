#include <Rcpp.h>
#include <cmath>

#include "CellBasedModel.h"

// update the R model after changes to the C model
void CellBasedModel::updateRModel(Rcpp::S4* rModel)
{
    rModel->slot("cells") = Rcpp::wrap(mPopulationRecord);
}

// run the entire model
void CellBasedModel::run()
{
    double time = 0.0;
    double recordTime = 0.0, outputTime = 0.0;

    Rprintf("\n");

    while (time <= mParams->runTime())
    {
        Rcpp::checkUserInterrupt();

        if (time >= recordTime) // record cell state
        {
            recordPopulation();
            recordTime = std::min(recordTime + mParams->recordIncrement(),
                mParams->runTime());
        }

        if (time >= outputTime) // output time and num cells
        {
            Rprintf("time = %.2f\n", floor(time));
            Rprintf("size = %d\n", size());

            outputTime = std::min(outputTime + mParams->outputIncrement(),
                mParams->runTime());
        }            

        oneTimeStep(time); // run the simulation for one time step
		time += mParams->timeIncrement();
    
        // ensures last time step happens at end of given runTime
        if (time > mParams->runTime() && time < mParams->runTime() +
        mParams->timeIncrement())
        {
            time = mParams->runTime();
        }
    }
}
