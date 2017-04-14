#include "DrasdoHohmeModel.h"
#include "../Core/Random.h"

DrasdoHohmeParameters::DrasdoHohmeParameters(Rcpp::S4* rModel)
: OffLatticeParameters(rModel)
{
    mNG      = rModel->slot("nG");
    mEpsilon = rModel->slot("epsilon");
    mDelta   = rModel->slot("delta");
}

DrasdoHohmeModel::DrasdoHohmeModel(Rcpp::S4* rModel)
: OffLatticeCellBasedModel(rModel)
{
    Parameters* temp = new DrasdoHohmeParameters(rModel);
    delete mParams;
    mParams = temp;
}

double DrasdoHohmeModel::growthRate(OffLatticeCell& cell) const
{
    return sqrt(cell.type().size())*(sqrt(2)-1) * mParams->timeIncrement()
        * (DH_PARAMS->nG() + 1) / (2 * cell.cycleLength());
}

void DrasdoHohmeModel::attemptTrial(OffLatticeCell& cell)
{
    double unif = Random::uniform(0,1);
    if (cell.phase() == INTERPHASE)
    {   
        if (unif <= (1 / (DH_PARAMS->nG() + 1)))
        {
            growth(cell);
        }
        else
        {
            translation(cell);
        }
    }
    else if (cell.phase() == MITOSIS)
    {
        if (unif <= (1 / (DH_PARAMS->nG() + 1)))
        {
            deformation(cell);
        }
        else if ((DH_PARAMS->nG() + 1) * unif <= 1 + DH_PARAMS->nG() / 2)
        {
            rotation(cell);
        }
        else if (!cell.readyToDivide())
        {
            translation(cell);
        }
    }
    else
    {
        throw std::runtime_error("invalid cell phase");
    }
}

bool DrasdoHohmeModel::acceptTrial(Energy preE, Energy postE,
unsigned preN, unsigned postN) const
{
    if (postN < preN)
    {
        return false;
    }
    else if (postE.first < preE.first)
    {
        return true;
    }
    else
    {
        double prob = exp(-DH_PARAMS->epsilon() * (postE.first -
            preE.first));
        return Random::uniform(0,1) < prob;
    }
}

Energy DrasdoHohmeModel::calculateHamiltonian(const OffLatticeCell& cell)
{
    double sum = 0.0;
    double maxSearch = cell.radius() + DH_PARAMS->maxRadius()
        + DH_PARAMS->delta();
    
    LocalCellIterator it = mCellPopulation.lbegin(
        cell.coordinates(), maxSearch);
    LocalCellIterator endIt = mCellPopulation.lend(
        cell.coordinates(), maxSearch);

    for (; it != endIt; ++it)
    {
        double dist = cell.distance(*it);
        if (cell != *it && dist <= DH_PARAMS->delta())
        {
            sum += pow(2 * dist / DH_PARAMS->delta() - 1, 2) - 1;
        }
    }
    return std::pair<double, bool> (sum, false);        
}

unsigned DrasdoHohmeModel::numNeighbors(const OffLatticeCell& cell)
{
    unsigned neighbors = 0;
    double maxSearch = cell.radius() + DH_PARAMS->maxRadius()
        + DH_PARAMS->delta();
    
    LocalCellIterator it = mCellPopulation.lbegin(
        cell.coordinates(), maxSearch);
    LocalCellIterator endIt = mCellPopulation.lend(
        cell.coordinates(), maxSearch);

    for (; it != endIt; ++it)
    {
        if (cell != *it && cell.distance(*it) <= DH_PARAMS->delta())
        {
            neighbors++;
        }
    }
    return neighbors;
}
