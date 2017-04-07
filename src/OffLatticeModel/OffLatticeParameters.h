#ifndef CIS_OFF_LATTICE_PARAMETERS_H
#define CIS_OFF_LATTICE_PARAMETERS_H

#include <vector>
#include <Rcpp.h>

#include "../Core/Parameters.h"
#include "OffLatticeRadiusSolver.h"

class OffLatticeParameters : public Parameters
{
public:

    // constructors
    OffLatticeParameters(Rcpp::List p) : Parameters(p) {}

    // basic off lattice parameters
    double maxDeform()      {return mParams["maxDeform"];}
    double maxTranslation() {return mParams["maxTranslation"];}
    double maxRotate()      {return mParams["maxRotate"];}
    double maxGrowth()      {return mParams["maxGrowth"];}
    
    // get the largest possible radius
    double maxRadius();
};

#endif
