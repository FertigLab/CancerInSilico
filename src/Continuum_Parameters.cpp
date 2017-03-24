#include <Rcpp.h>

#include "Continuum_Parameters.h"

void ContinuumParameters::CalculateTimeIncrement()
{
    Rcpp::NumericVector cycleDist = mParams["cycleLengthDist"];
    double minCycle = Rcpp::min(cycleDist);

    double t1 = delta() / (4 * nG() * (4 - pow(2, 0.5)));
    double t2 = delta() * (minCycle - 1) / (8 * nG() * (pow(2, 0.5) - 1));

    mParams["timeIncrement"] = std::min(t1, t2);
}
