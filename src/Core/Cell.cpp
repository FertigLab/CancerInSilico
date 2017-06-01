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

// apply the effect of a drug on the cell and record the action
double Cell::applyDrug(const Drug& drug)
{
    double oldLength = cycleLength();
    mCycleLength = drug.cycleLengthEffect(type(), cycleLength());
    mDrugApplied |= 1 << drug.id();
    return oldLength - mCycleLength;
}

