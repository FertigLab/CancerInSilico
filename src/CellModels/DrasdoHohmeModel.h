#ifndef CIS_CONTINUUM_MODEL_H
#define CIS_CONTINUUM_MODEL_H

#include "../OffLatticeModel/OffLatticeCellBasedModel.h"
#include "../OffLatticeModel/OffLatticeParameters.h"

#define DH_PARAMS   static_cast<DrasdoHohmeParameters*>(mParams)

class DrasdoHohmeParameters : public OffLatticeParameters
{
private:

    // convert drasdo parameters to general off lattice model parameters
    void processParameters() {}

public:

    // constructor
    DrasdoHohmeParameters(Rcpp::List rP) : OffLatticeParameters(rP)
    {
        processParameters();
    }

    // get parameters
    double nG()             {return mParams["nG"];}
    double epsilon()        {return mParams["epsilon"];}
    double delta()          {return mParams["delta"];}
};

class DrasdoHohmeModel : public OffLatticeCellBasedModel
{
public:

    DrasdoHohmeModel(DrasdoHohmeParameters* p)
        : OffLatticeCellBasedModel(p) {}

    void attemptTrial(OffLatticeCell&);
    bool acceptTrial(Energy, Energy, unsigned, unsigned) const;
    Energy calculateHamiltonian(const OffLatticeCell&);
    unsigned numNeighbors(const OffLatticeCell&);
};

#endif
