#ifndef CIS_OFF_LATTICE_PARAMETERS_H
#define CIS_OFF_LATTICE_PARAMETERS_H

// stores parameters of off lattice model

#include <vector>
#include <Rcpp.h>

#include "../Core/Parameters.h"
#include "OffLatticeRadiusSolver.h"

class OffLatticeParameters : public Parameters
{
private:

    // geometric parameters
    double mMaxDeformation;
    double mMaxTranslation;
    double mMaxRotation;

public:

    // constructor
    OffLatticeParameters(Rcpp::S4*);

    // basic off lattice parameters
    double maxDeformation() {return mMaxDeformation;}
    double maxTranslation() {return mMaxTranslation;}
    double maxRotation()    {return mMaxRotation;}
    
    // get the largest possible radius
    double maxRadius();
};

#endif
