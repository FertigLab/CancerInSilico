#ifndef CIS_DRUG_H
#define CIS_DRUG_H

// contains all info about a drug and it's interactions with cells

#include <Rcpp.h>

#include "CellType.h"
#include "CellPhase.h"

class Drug
{
private:

    // used to differentiate between drugs
    unsigned mID;

    // time this drug is added to the simulation
    double mTimeAdded;

    // needed for drug effect on cycle length
    Rcpp::S4 mDrugClass;

public:

    Drug() {}
    Drug(unsigned, const Rcpp::S4&);

    // getters  
    unsigned id()      const {return mID;}
    double timeAdded() const {return mTimeAdded;}
    
    // calculate effect this drug has on the cycle length of a certain cell
    double cycleLengthEffect(const CellType&, double) const;
};

#endif
