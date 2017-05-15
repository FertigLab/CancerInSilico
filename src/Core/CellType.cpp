#include <Rcpp.h>
#include <algorithm>

#include "CellType.h"
#include "Random.h"

// construct cell type from equivalent class in R
CellType::CellType(unsigned id, const Rcpp::S4& type)
{
    mID = id;
    mName = Rcpp::as<std::string>(type.slot("name"));
    mSize = Rcpp::as<double>(type.slot("size"));
    mMinCycle = Rcpp::as<double>(type.slot("minCycle"));
    mCellTypeClass = type;
}

// get a random cycle length
double CellType::cycleLength() const
{
    Rcpp::Function cl = mCellTypeClass.slot("cycleLength");
    return Rcpp::as<double>(cl());
}
