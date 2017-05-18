#include <Rcpp.h>
#include <iostream>

#include "Drug.h"
#include "CellType.h"
#include "Cell.h"

// construct drug with equivalent R class
Drug::Drug(const Rcpp::S4& drug)
{
    mTimeAdded = Rcpp::as<double>(drug.slot("timeAdded"));
    mDrugClass = drug;
}

// calculate drug effect on cycle length of a cell
double Drug::cycleLengthEffect(const CellType& type, double cycle) const
{ 
    // call R function
    Rcpp::Function effect = mDrugClass.slot("cycleLengthEffect");
    return Rcpp::as<double>(effect(Rcpp::wrap(type.name()), cycle));
}
