#include <Rcpp.h>

#include "Cell.h"
#include "CellType.h"

// constructor for initial cells
Cell::Cell(CellType type) : mType(type)
{
    mCycleLength = mType.cycleLength();
    mDrugApplied = 0;
    mReadyToDivide = false;
    mPhase = INTERPHASE;
}

void Cell::applyDrug(const Drug& drug)
{
    mCycleLength = drug.cycleLengthEffect(type(), cycleLength());
    mDrugApplied |= 1 << drug.id();
}

