#include <Rcpp.h>

#include "Cell.h"
#include "CellType.h"

// constructor for initial cells
Cell::Cell(const CellType& type)
{
    mType = &type;
    mCycleLength = mType->cycleLength();
    mDrugApplied = 0;
    mReadyToDivide = false;
    mPhase = INTERPHASE;
}

void Cell::applyDrug(const Drug& d)
{

}

