#include <Rcpp.h>

#include "Drug.h"
#include "CellType.h"
#include "Cell.h"

Drug::Drug(unsigned id, const Rcpp::S4& drug)
{
    mID = id;
    mTimeAdded = Rcpp::as<double>(drug.slot("timeAdded"));
    mDrugClass = drug;
}

double Drug::cycleLengthEffect(const CellType& type, double cycleLength,
CellPhase phase) const
{
    Rcpp::Function effect = mDrugClass.slot("cycleLengthEffect");
    return Rcpp::as<double>(effect(type.name(), cycleLength, (int) phase));
}
