#include <Rcpp.h>

#include "OffLatticeParameters.h"
#include "../Core/CellType.h"

// construct from equivalent R class
OffLatticeParameters::OffLatticeParameters(Rcpp::S4* rModel)
: Parameters(rModel)
{
    mMaxDeformation = rModel->slot("maxDeformation");
    mMaxTranslation = rModel->slot("maxTranslation");
    mMaxRotation    = rModel->slot("maxRotation");
}

// find largest possible radius across cell types
double OffLatticeParameters::maxRadius()
{
    double maxRad = 0.0;
    std::vector<CellType>::iterator it = mCellTypes.begin();
    for (; it != mCellTypes.end(); ++it)
    {
        if (maxRad < it->size()) {maxRad = it->size();}
    }
    return sqrt(2 * maxRad);
}
