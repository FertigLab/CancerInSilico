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
        oneMCstep(time);
    }
}

void ContinuumModel::oneMCstep(double time)
{
    Cell cell = mCellPopulation.randomValue();
    doTrial(cell);
    checkMitosis(cell);
}

void ContinuumModel::doTrial(Cell& cell)
{
    double energy = totalEnergy(cell);
    unsigned numNeighbors = numNeighbors(cell);

    Cell orig = cell;
    bool growth = attemptTrial(cell);

    if (checkOverlap(cell) || checkBoundary(cell))
    {
        cell = orig;
    }
    else if (!growth)
    {
        mCellPopulation.update(orig.coord(), cell.coord());
        if (!acceptTrial(energy, numNeighbors, cell))
        {
            mCellPopulation.update(cell.coord(), orig.coord());
            cell = orig;
        }
    }
}

void ContinuumModel::attemptTrial(ContinuumCell& cell)
{
    double unif = Random::uniform(0,1);
    bool growthTrial = unif <= 1.0 / (mParams->nG() + 1.0);

    if (cell.phase() != MITOSIS)
    {
        growthTrial ? Growth(cell) : Translation(cell);
    }
    else
    {
        if (growthTrial)
        {
            Deformation(cell);
        }
        else if ((mParams->nG() + 1.0) * unif <= 1.0 + mParams->nG() / 2.0))
        {
            Rotation(cell);
        }
        else if (!cell.readyToDivide())
        {
            Translation(cell);
        }
    }
    return growthTrial;            
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

bool ContinuumModel::acceptTrial(double energy, unsigned prevNeighbors,
const Cell& cell)
{
    double deltaEnergy = CalculateTotalEnergy(cell) - energy;

    if (CalculateNumNeighbors(cell) < prevNeighbors)
    {
        return false;
    }
    else if (deltaEnergy <= 0.0)
    {
        return true;
    }
    else
    {
        return Random::uniform(0,1) < exp(-1 * deltaEnergy);
    }
}

double ContinuumModel::CalculateTotalEnergy(Cell& cell)
{
    double sum = 0.0;
    SquareLattice<Cell>::local_iterator it = mCellPopulation.begin(
        cell.coord(), mParams->delta() + 1);
    SquareLattice<Cell>::local_iterator endIt = mCellPopulation.end(
        cell.coord(), mParams->delta() + 1);

    for (; it != endIt; ++it)
    {
        double dist = cell.distance(*it);
        sum +=  dist > mParams->delta() ? 0.0 : mParams->epsilon()
            * (pow(2 * dist / mParams->delta() - 1, 2) - 1);
    }
    return sum;
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

