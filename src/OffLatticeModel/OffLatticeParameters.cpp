#include <Rcpp.h>

#include "OffLatticeParameters.h"
#include "../Core/CellType.h"

double OffLatticeParameters::maxRadius()
{
    double maxRad = 0.0;
    std::vector<CellType>::iterator it = mCellTypes.begin();
    for (; it != mCellTypes.end(); ++it)
    {
        if (maxRad < it->size())
        {
            maxRad = it->size();
        }
    }
    return sqrt(2 * maxRad);
}
