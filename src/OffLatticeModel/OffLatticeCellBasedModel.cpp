#include "OffLatticeCellBasedModel.h"

OffLatticeCellBasedModel::OffLatticeCellBasedModel(OffLatticeParameters* p)
{
    mParams = p;
    mCellPopulation.setWidth(1.0);
}

void OffLatticeCellBasedModel::oneTimeStep(double time)
{
    // update all drugs in the system
    updateDrugs(time);
    
    // do N monte carlo steps
    unsigned sz = mCellPopulation.size();
    for (unsigned i = 0; i < sz; ++i)
    {
        oneMCStep();
    }
}

void OffLatticeCellBasedModel::oneMCStep()
{
    OffLatticeCell& cell = mCellPopulation.randomValue();
    doTrial(cell);
    checkMitosis(cell);
}

void OffLatticeCellBasedModel::updateDrugs(double time)
{
    // cell population iterator
    CellIterator cellIt = mCellPopulation.begin();

    for (; cellIt != mCellPopulation.end(); ++cellIt)
    {
        // drug iterator
        std::vector<Drug>::iterator drugIt = mParams->drugsBegin();
        for (; drugIt != mParams->drugsEnd(); ++drugIt)
        {
            // check if drug has been applied
            if (!(*cellIt).drugApplied((*drugIt).id())
                && time >= (*drugIt).timeAdded())
            {
                (*cellIt).applyDrug(*drugIt);
            }
        }
    }
}

void OffLatticeCellBasedModel::doTrial(OffLatticeCell& cell)
{
    OffLatticeCell orig = cell;
    Energy preE = calculateHamiltonian(cell);
    unsigned preN = numNeighbors(cell);
    if (preE.second) {throw std::runtime_error("infinite hamiltonian");}

    attemptTrial(cell);
    Energy postE = calculateHamiltonian(cell);
    unsigned postN = numNeighbors(cell);

    if (checkOverlap(cell) || checkBoundary(cell))
    {
        cell = orig;
    }
    else
    {
        mCellPopulation.update(orig.coordinates(), cell.coordinates());
        if (!acceptTrial(preE, postE, preN, postN))
        {
            mCellPopulation.update(cell.coordinates(), orig.coordinates());
            cell = orig;
        }
    }
}

void OffLatticeCellBasedModel::checkMitosis(OffLatticeCell& cell)
{
    double fullRad = sqrt(2 / cell.type().size());
    if (cell.phase() == MITOSIS && cell.radius() == fullRad)
    {
        OffLatticeCell daughter(cell.type());
        cell.divide(daughter);
    }
}

bool OffLatticeCellBasedModel::checkOverlap(const OffLatticeCell& cell)
{
    double maxSearch = std::max(OL_PARAMS->maxTranslation(),
        OL_PARAMS->maxRadius());
    LocalCellIterator it =
        mCellPopulation.lbegin(cell.coordinates(), maxSearch);
    LocalCellIterator endIt =
        mCellPopulation.lend(cell.coordinates(), maxSearch);

    for (; it != endIt; ++it)
    {
        if (cell.distance(*it) < 0)
        {
            return true;
        }
    }
    return false;
}

bool OffLatticeCellBasedModel::checkBoundary(const OffLatticeCell& cell)
{
    Point<double> origin(0,0);
    double b = mParams->boundary();

    // return true if cell is farther from center than the boundary line
    return (cell.centers().first.distance(origin) + cell.radius() > b 
        || cell.centers().second.distance(origin) + cell.radius() > b);
}

void OffLatticeCellBasedModel::growth(OffLatticeCell& cell)
{
    double growthRate = 0;
    double growth = Random::uniform(0, growthRate);
    cell.setRadius(std::min(cell.type().size() * pow(2,0.5), 
        cell.radius() + growth));
    
    if (cell.radius() == sqrt(2 / cell.type().size()))
    {
        cell.setPhase(MITOSIS);
    }
}

void OffLatticeCellBasedModel::translation(OffLatticeCell& cell)
{
    double len = OL_PARAMS->maxTranslation() * sqrt(Random::uniform(0,1));
    double dir = Random::uniform(0, 2 * M_PI);
    cell.setCoordinates(Point<double>(len * cos(dir), len * sin(dir)));
}

void OffLatticeCellBasedModel::deformation(OffLatticeCell& cell)
{
    double deform = Random::uniform(0, OL_PARAMS->maxDeform());
    cell.setAxisLength(std::min(4.0, cell.axisLength() + deform));

    if (cell.axisLength() == 4.0 * cell.type().size())
    {
        cell.setReadyToDivide(true);
    }
}

void OffLatticeCellBasedModel::rotation(OffLatticeCell& cell)
{
    cell.setAxisAngle(Random::uniform(-OL_PARAMS->maxRotate(),
        OL_PARAMS->maxRotate()));
}

void OffLatticeCellBasedModel::recordPopulation()
{
    // the vector to hold the current population
    std::vector<double> current;
    
    // cell population iterator
    CellIterator it = mCellPopulation.begin();
  
    // loop through each cell, store info
    for (; it != mCellPopulation.end(); ++it)
    {
        current.push_back((*it).coordinates().x);
        current.push_back((*it).coordinates().y);
        current.push_back((*it).radius());
        current.push_back((*it).axisLength());
        current.push_back((*it).axisAngle());
        current.push_back((*it).cycleLength());
        current.push_back((*it).phase());
        current.push_back((*it).type().id());
    }

    // add current population to record
    mPopulationRecord.push_back(current);
}

