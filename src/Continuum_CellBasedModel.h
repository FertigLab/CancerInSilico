#ifndef CIS_CONTINUUM_MODEL_H
#define CIS_CONTINUUM_MODEL_H

#include "CellBasedModel.h"
#include "SquareLattice.h"
#include "Parameters.h"

class ContinuumModel : public CellBasedModel
{
private:

    SquareLattice<Cell> mCellPopulation;
    
public:

    ContinuumModel(ContinuumParameters*);

    void oneTimeStep(double time);
    void updateDrugs(double time);
    void recordPopulation();
};

#endif
