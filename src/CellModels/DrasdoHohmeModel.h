#ifndef CIS_DRASDO_HOHME_MODEL_H
#define CIS_DRASDO_HOHME_MODEL_H

// Off-Lattice Monte-Carlo cell-based model based
// on the work of Drasdo, Hohme (2003)

#include "../OffLatticeModel/OffLatticeCellBasedModel.h"
#include "../OffLatticeModel/OffLatticeParameters.h"

#define DH_PARAMS   static_cast<DrasdoHohmeParameters*>(mParams)

class DrasdoHohmeParameters : public OffLatticeParameters
{
private:

    // parameters defined in paper - control cell mechanics
    double mNG;
    double mEpsilon;
    double mDelta;

public:

    // constructor
    DrasdoHohmeParameters(Rcpp::S4*);

    // get parameters
    double nG()         {return mNG;}
    double epsilon()    {return mEpsilon;}
    double delta()      {return mDelta;}
};

class DrasdoHohmeModel : public OffLatticeCellBasedModel
{
public:

    DrasdoHohmeModel(Rcpp::S4*);

    double maxGrowth(OffLatticeCell&) const;

    // functions for attempting/accepting monte carlo trials
    bool attemptTrial(OffLatticeCell&);
    bool acceptTrial(Energy, Energy, unsigned, unsigned) const;
    Energy calculateHamiltonian(const OffLatticeCell&);
    unsigned numNeighbors(const OffLatticeCell&);
};

#endif
