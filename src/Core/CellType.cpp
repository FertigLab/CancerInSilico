#include <Rcpp.h>
#include <algorithm>

#include "CellType.h"
#include "Random.h"

#define DIST_SIZE 1000

CellType::CellType(unsigned id, const Rcpp::S4& type)
{
    mID = id;
    mName = Rcpp::as<std::string>(type.slot("name"));
    mSize = Rcpp::as<double>(type.slot("size"));
    mMinCycle = Rcpp::as<double>(type.slot("minCycle"));
    mCellTypeClass = type;
}

double CellType::cycleLength() const
{
    Rcpp::Function cl = mCellTypeClass.slot("cycleLength");
    return Rcpp::as<double>(cl());
}
