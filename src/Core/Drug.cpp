#include <Rcpp.h>
#include <iostream>

#include "Drug.h"
#include "CellType.h"
#include "Cell.h"

Drug::Drug(unsigned id, const Rcpp::S4& drug)
{
    mID = id;
    mTimeAdded = Rcpp::as<double>(drug.slot("timeAdded"));
    mDrugClass = drug;
}

double Drug::cycleLengthEffect(CellType type, double cycleLength) const
{ 
    Rcpp::Function effect = mDrugClass.slot("cycleLengthEffect");
    return Rcpp::as<double>(effect(Rcpp::wrap(type.name()), cycleLength));
}
