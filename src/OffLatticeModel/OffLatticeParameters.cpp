#include <Rcpp.h>

#include "Continuum_Parameters.h"
#include "CellType.h"

void ContinuumParameters::CalculateTimeIncrement()
{
    double minCycle = std::numeric_limits<double>::max();
    std::vector<CellType>::iterator it = mCellTypes.begin();
    for (; it != mCellTypes.end(); ++it)
    {
        minCycle = std::min(minCycle, (*it).minCycleLength());
    }

    double t1 = delta() / (4 * nG() * (4 - pow(2, 0.5)));
    double t2 = delta() * (minCycle - 1) / (8 * nG() * (pow(2, 0.5) - 1));

    mParams["timeIncrement"] = std::min(t1, t2);
}
