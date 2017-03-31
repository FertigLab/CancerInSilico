#include <Rcpp.h>

#include "Cell.h"

// constructor for initial cells
Cell::Cell(Rcpp::S4* type)
{
    mCellType = type;

    Rcpp::Function cycleLengthDist = type->slot("cycleLength");
    mProperties["CycleLength"] = Rcpp::as<double>(cycleLengthDist());
}

// constructor for daughter cell
Cell::Cell(const Cell* parent)
{
    divide(this, parent);
}

