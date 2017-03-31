#include "CellType.h"
#include "Random.h"

CellType::CellType(unsigned id, Rcpp::S4 type) : mID(id)
{
    mSize = type.slot("size");
    mInheritCycle = type.slot("inheritCycleLength");
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
