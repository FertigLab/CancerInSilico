#ifndef CIS_DRASDO_HOHME_MODEL_H
#define CIS_DRASDO_HOHME_MODEL_H

#include "../OffLatticeModel/OffLatticeCellBasedModel.h"
#include "../OffLatticeModel/OffLatticeParameters.h"

#define DH_PARAMS   static_cast<DrasdoHohmeParameters*>(mParams)

class DrasdoHohmeParameters : public OffLatticeParameters
{
private:

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

    bool attemptTrial(OffLatticeCell&);
    bool acceptTrial(Energy, Energy, unsigned, unsigned) const;
    Energy calculateHamiltonian(const OffLatticeCell&);
    unsigned numNeighbors(const OffLatticeCell&);
};

#endif
