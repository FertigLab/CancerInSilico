#include <Rcpp.h>

#include "Drug.h"

Drug::Drug(unsigned id, Rcpp::S4 drug) : mID(id)
{
    mTimeAdded = drug.slot("timeAdded");
    mInheritedEffect = drug.slot("inheritedEffect");
    mCycleLengthEffect = drug.slot("cycleLengthEffect");
}

double cycleLengthEffect(const CellType* type, double cycleLength,
CellPhase phase) const
{
    return mCycleLengthEffect(&type, cycleLength, phase);
}
