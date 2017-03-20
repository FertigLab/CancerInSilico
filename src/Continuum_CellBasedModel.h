#ifndef CIS_CONTINUUM_MODEL_H
#define CIS_CONTINUUM_MODEL_H

#include "CellBasedModel.h"
#include "SpatialHash.h"
#include "Parameters.h"

class ContinuumModel : public CellBasedModel
{
private:

    SpatialHash<Cell> mCellPopulation;

public:

    ContinuumModel(Parameters*);

    void oneTimeStep(double time);
    void updateDrugs(double time);
};

#endif
