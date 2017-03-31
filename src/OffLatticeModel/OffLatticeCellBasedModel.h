#ifndef CIS_CONTINUUM_MODEL_H
#define CIS_CONTINUUM_MODEL_H

#include <utility>

#include "CellBasedModel.h"
#include "SquareLattice.h"
#include "Parameters.h"

typedef std::pair<double, bool> Energy;

class OffLatticeCellBasedModel : public CellBasedModel
{
private:

    SquareLattice<Cell> mCellPopulation;
    
public:

    OffLatticeCellBasedModel(OffLatticeParameters*);

    void oneTimeStep(double time);
    void oneMonteCarloStep(double time);

    void doTrial(OffLatticeCell&);

    virtual void attemptTrial(OffLatticeCell&) = 0;
    virtual bool acceptTrial(Energy, Energy) const = 0;
    virtual double calculateHamiltonian(const OffLatticeCell&) const = 0;

    void checkMitosis(Cell&);
    void numNeighbors(const Cell&) const;
    void updateDrugs(double time);
    void recordPopulation() const;

    bool checkOverlap(const Cell&) const;
    bool checkBoundary(const Cell&) const;
};

#endif
