#include "Continuum_CellBasedModel.h"

ContinuumModel::ContinuumModel(ContinuumParameters* par)
: CellBasedModel(par)
{


}

void ContinuumModel::oneTimeStep(double time)
{
    // update all drugs in the system
    updateDrugs();
    
    // do N monte carlo steps
    unsigned sz = mCellPopulation.size();
    for (unsigned i = 0; i < sz; ++i)
    {
        oneMonteCarlostep(time);
    }
}

void ContinuumModel::oneMCstep(double time)
{
    OffLatticeCell& cell = mCellPopulation.randomValue();
    doTrial(cell);
    checkMitosis(cell);
}

void ContinuumModel::doTrial(OffLatticeCell& cell)
{
    Cell orig = cell;
    Energy pre = calculateHamiltonian(cell);
    if (pre.second) {throw std::runtime_error("infinite hamiltonian");}

    attemptTrial(cell);
    Energy post = calculateHamiltonian(cell);

    if (checkOverlap(cell) || checkBoundary(cell))
    {
        cell = orig;
    }
    else
    {
        mCellPopulation.update(orig.coord(), cell.coord());
        if (!acceptTrial(pre, post))
        {
            mCellPopulation.update(cell.coord(), orig.coord());
            cell = orig;
        }
    }
}

bool checkOverlap(const Cell& cell)
{
    double maxSearch = std::max(mParams->maxTranslation(),
        mParams->maxRadius());
    SquareLattice<Cell>::local_iterator it = mCellPopulation.begin(
        cell.coord(), maxSearch);
    SquareLattice<Cell>::local_iterator endIt = mCellPopulation.end(
        cell.coord(), maxSearch);

    for (; it != endIt; ++it)
    {
        if (cell.distance(*it) < 0)
        {
            return true;
        }
    }
    return false;
}

bool checkBoundary(const ContinuumCell& cell)
{
    Point<double> origin(0,0);
    double b = mParams->boundary();

    /* return true if cell is farther from center than the boundary line */
    return (cell.centers().first.dist(origin) + cell.radius() > b 
        || cell.centers().second.dist(origin) + cell.radius() > b);
}

unsigned ContinuumModel::numNeighbors(const Cell& cell)
{
    unsigned neighbors = 0;
    SquareLattice<Cell>::local_iterator it = mCellPopulation.begin(
        cell.coord(), mParams->delta() + 1);
    SquareLattice<Cell>::local_iterator endIt = mCellPopulation.end(
        cell.coord(), mParams->delta() + 1);

    for (; it != endIt; ++it)
    {
        if (cell.distance(*it) <= mParams->delta())
        {
            neighbors++;
        }
    }
    return neighbors;
}

void ContinuunModel::checkMitosis(Cell& cell)
{


}

void ContinuumModel::updateDrugs(double time)
{
    // cell population iterator
    SquareLattice<Cell>::full_iterator cellIt = mCellPopulation.begin();
    for (; cellIt != mCellPopulation.end(); ++cellIt)
    {
        // drug iterator
        std::vector<Drug>::iterator drugIt = mParams->drugs().begin();
        for (; drugIt != mParams->drugs().end(); ++drugIt)
        {
            // check if drug has been applied
            if (!(cellIt).drugApplied((*drugIt).id()))
            {
                (*cellIt).applyDrug(*drugIt);
            }
        }
    }
}

void ContinuumModel::recordPopulation()
{
    // the vector to hold the current population
    std::vector<double> current;
    
    // cell population iterator
    SquareLattice<Cell>::full_iterator it = mCellPopulation.begin();
  
    // loop through each cell, store info
    for (; it != mCellPopulation.end(); ++it)
    {
        current_pop.push_back((*it).coord().x);
        current_pop.push_back((*it).coord().y);
        current_pop.push_back((*it).radius());
        current_pop.push_back((*it).axisLength());
        current_pop.push_back((*it).axisAngle());
        current_pop.push_back((*it).GetGrowth());
        current_pop.push_back((*it).cellType());
    }

    // add current population to record
    mPopulationRecord.push_back(current);
}

