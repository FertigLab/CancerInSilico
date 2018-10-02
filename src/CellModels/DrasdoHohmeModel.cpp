#include "DrasdoHohmeModel.h"
#include "../Core/Random.h"

// constructor
DrasdoHohmeModel::DrasdoHohmeModel(Rcpp::S4* rModel)
: OffLatticeCellBasedModel(rModel)
{
    // store parameters
    mNG      = rModel->slot("nG");
    mEpsilon = rModel->slot("epsilon");
    mDelta   = rModel->slot("delta");
}

// calculate the maximum growth (radius increase) possible in one time step
double DrasdoHohmeModel::maxGrowth(OffLatticeCell& cell) const
{
    return 2 * timeIncrement() * nG() * sqrt(cell.type().size())
        * (sqrt(2) - 1) / (cell.cycleLength() - 2);
}

// calculate the maximum deformation (axis increase) possible in one time step
double DrasdoHohmeModel::maxDeformation(OffLatticeCell& cell) const
{
    return 2 * timeIncrement() * nG() * sqrt(cell.type().size())
        * (2 - sqrt(2));
}

// attempt a single monte carlo update for this cell
bool DrasdoHohmeModel::attemptTrial(OffLatticeCell& cell)
{
    // nG determines how likely growth (deform) trial is
    double unif = Random::uniform(0,1);
    if (cell.phase() == INTERPHASE)
    {   
        if (unif <= 1 / nG())
        {
            growth(cell);
            return true;
        }
        else
        {
            translation(cell);
        }
    }
    else if (cell.phase() == MITOSIS)
    {
        if (unif <= 1 / nG())
        {
            deformation(cell);
            return true;
        }
        else if (nG() * unif <= 1 + nG() / 2)
        {
            rotation(cell);
        }
        else
        {
            translation(cell);
        }
    }
    else
    {
        throw std::runtime_error("invalid cell phase");
    }
    return false;
}

// accept/reject trial based on change in hamiltonian (see paper)
bool DrasdoHohmeModel::acceptTrial(Energy preE, Energy postE,
unsigned preN, unsigned postN) const
{
    if (postN < preN) // lost a neighbor
    {
        return false;
    }
    else if (postE.first < preE.first) // hamiltonian decreased
    {
        return true;
    }
    else
    {
        // prob of acceptance, when hamiltonian increases
        double prob = exp(epsilon() * (postE.first - preE.first));
        return Random::uniform(0,1) < prob;
    }
}

// calculate hamiltonian locally around this cell
NeighborInfo DrasdoHohmeModel::getNeighborInfo(const OffLatticeCell& cell)
{
    double maxSearch = 2 * cell.radius() + 2 * maxRadius() + delta();
    double en = 0.0;
    unsigned neighbors = 0;

    // iterate locally around this cell - within the search radius
    LocalCellIterator it = mCellPopulation.lbegin(
        cell.coordinates(), maxSearch);
    LocalCellIterator endIt = mCellPopulation.lend(
        cell.coordinates(), maxSearch);

    for (; it != endIt; ++it)
    {
        if (cell != *it)
        {
            double dist = cell.distance(*it);
            if (dist <= delta())
            {
                neighbors++;

                // formula for hamiltonian taken from Drasdo, Hohme paper
                en += pow(2 * dist / delta() - 1, 2) - 1;
            }
        }
    }
    return NeighborInfo(en, false, neighbors);
}
