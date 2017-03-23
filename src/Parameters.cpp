#if 0

#include <cmath>
#include <Rcpp.h>
#include <algorithm>
#include <exception>

#include "Parameters.h"

Parameters::Parameters(Rcpp::List Rparams)
{
    mParams = Rparams;
    CalculateTimeIncrement();
}

void Parameters::StoreTimeIncrement()
{
    Rcpp::NumericVector cycleDist = mParams["cycleLengthDist"];
    double minCycle = Rcpp::min(cycleDist);

    double t1 = delta() / (4 * nG() * (4 - pow(2, 0.5)));
    double t2 = delta() * (minCycle - 1) / (8 * nG() * (pow(2, 0.5) - 1));

    mParams["timeIncrement"] = std::min(t1, t2);
}

void Parameters::StoreDrasdoParameters()
{
    mParams["maxDeform"] = 2 * timeIncrement() * nG() * (4 - pow(2, 0.5));
    mParams["maxTranslation"] = delta() / 2;
    mParams["maxRotate"] = acos((16 + pow(delta(), 2) - 4 * delta()) / 16);
}

Rcpp::S4 Parameters::getRandomCellType()
{

}
#endif
