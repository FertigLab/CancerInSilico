#ifndef CIS_DRUG_H
#define CIS_DRUG_H

#include <Rcpp.h>

#include "CellType.h"
#include "CellPhase.h"

class Drug
{
private:

    // numerical id of this type
    unsigned mID;

    // time this drug is added to the simulation
    double mTimeAdded;

    // needed for drug effect on cycle length
    Rcpp::S4 mDrugClass;

public:

    Drug() {}
    Drug(unsigned, const Rcpp::S4&);

    unsigned id() const {return mID;}
    double timeAdded() const {return mTimeAdded;}
    
    double cycleLengthEffect(CellType, double) const;
};

#endif
