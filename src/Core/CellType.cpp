#include <Rcpp.h>
#include <algorithm>

#include "CellType.h"
#include "Random.h"

CellType::CellType(unsigned id, const Rcpp::S4& type)
{
    mID = id;
    mSize = Rcpp::as<double>(type.slot("size"));
    mName = Rcpp::as<std::string>(type.slot("name"));
    Rcpp::Function cl = type.slot("cycleLength");

    for (unsigned i = 0; i < DIST_SIZE; ++i)
    {
        mCycleLength[i] = Rcpp::as<double>(cl());
    }
}

double CellType::cycleLength() const
{
    return mCycleLength[Random::uniformInt(0, DIST_SIZE - 1)];
}

double CellType::minCycleLength() const
{
    return *std::min_element(mCycleLength, mCycleLength + DIST_SIZE);
}
