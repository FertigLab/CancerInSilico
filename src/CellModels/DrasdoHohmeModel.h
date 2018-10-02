#ifndef CIS_DRASDO_HOHME_MODEL_H
#define CIS_DRASDO_HOHME_MODEL_H

// Off-Lattice Monte-Carlo cell-based model based
// on the work of Drasdo, Hohme (2003)

#include "../OffLatticeModel/OffLatticeCellBasedModel.h"

class DrasdoHohmeModel : public OffLatticeCellBasedModel
{
protected:

    // parameters defined in paper - control cell mechanics
    double mNG;
    double mEpsilon;
    double mDelta;

public:

    // constructor from R class
    DrasdoHohmeModel(Rcpp::S4*);

    // get parameters
    double nG()      const  {return mNG;}
    double epsilon() const  {return mEpsilon;}
    double delta()   const  {return mDelta;}

    // largest amount radius can grow in single step
    double maxGrowth(OffLatticeCell&) const;
    double maxDeformation(OffLatticeCell&) const;

    // functions for attempting/accepting monte carlo trials
    bool attemptTrial(OffLatticeCell&);
    bool acceptTrial(Energy, Energy, unsigned, unsigned) const;
    NeighborInfo getNeighborInfo(const OffLatticeCell&);
};

#endif
