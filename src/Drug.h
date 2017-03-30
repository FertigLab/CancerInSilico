#ifndef CIS_DRUG_H
#define CIS_DRUG_H

class Drug
{
private:

private:

    // numerical id of this type
    unsigned mID;

    // time this drug is added to the simulation
    double mTimeAdded;

    // whether or not the same effect is inherited through generations
    bool mInheritedEffect;

    // drug effect on cycle length
    Rcpp::Function mCycleLengthEffect;

public:

    Drug(unsigned, Rcpp::S4);

    unsigned id() const {return mID;}
    double timeAdded() const {return mSize;}
    bool inheritEffect() const {return mInheritCycle;}
    
    double cycleLengthEffect(const CellType*, double, CellPhase) const;
};

#endif
