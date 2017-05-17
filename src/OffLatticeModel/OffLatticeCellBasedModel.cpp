#include "OffLatticeCellBasedModel.h"

// private helper function - get random point in a circle
static Point<double> getRandomPoint(double radius)
{
    double dist = Random::uniform(0,1);
    double ang = Random::uniform(0, 2 * M_PI);
    double x = radius * sqrt(dist) * cos(ang);
    double y = radius * sqrt(dist) * sin(ang);
    return Point<double>(x,y);
}

// TODO: add burn in period without any growth steps
// constructor
OffLatticeCellBasedModel::OffLatticeCellBasedModel(Rcpp::S4* rModel)
: CellBasedModel(rModel)
{
    // setup parameters and lattice structure
    Parameters* temp = new OffLatticeParameters(rModel);
    delete mParams;
    mParams = temp;
    mCellPopulation.setWidth(sqrt(2) - 0.001);

    // seed default cells
    std::vector<OffLatticeCell> defaultCells;
    double area = 0.0;
    for (unsigned i = 0; i < mParams->initialNum(); ++i)
    {
        OffLatticeCell cell (mParams->randomCellType());
        if (!mParams->syncCycles()) {cell.gotoRandomCyclePoint();}
        defaultCells.push_back(cell);
        area += cell.area();
    }

    // calculate boundary
    double seedBoundary = sqrt(area / (M_PI * mParams->density()));
    if (mParams->boundary() > 0) {mParams->setBoundary(seedBoundary);}
    
    // place cells randomly
    std::vector<OffLatticeCell>::iterator it = defaultCells.begin();
    for (; it != defaultCells.end(); ++it)
    {
        do
        {
            (*it).setCoordinates(getRandomPoint(seedBoundary));
            Rcpp::checkUserInterrupt();

        } while (checkOverlap(*it) || checkBoundary(*it));
        
        mCellPopulation.insert((*it).coordinates(), *it);
    }            
}

// run model for one time step
void OffLatticeCellBasedModel::oneTimeStep(double time)
{
    // update all drugs in the system
    updateDrugs(time);

    // do N monte carlo steps
    unsigned sz = size();
    for (unsigned i = 0; i < sz; ++i)
    {
        oneMCStep();
    }
}

// execute a single monte carlo step
void OffLatticeCellBasedModel::oneMCStep()
{
    OffLatticeCell& cell = mCellPopulation.randomValue();
    doTrial(cell);
    checkMitosis(cell);
}

// update drugs in the system
void OffLatticeCellBasedModel::updateDrugs(double time)
{
    CellIterator cellIt = mCellPopulation.begin();
    for (; cellIt != mCellPopulation.end(); ++cellIt)
    {
        std::vector<Drug>::iterator drugIt = mParams->drugsBegin();
        for (; drugIt != mParams->drugsEnd(); ++drugIt)
        {
            // check if drug hasn't been applied and it's time to apply
            if (!(*cellIt).drugApplied((*drugIt).id())
                && time >= (*drugIt).timeAdded())
            {
                (*cellIt).applyDrug(*drugIt);
            }
        }
    }
}

// attempt a single monte carlo trial on this cell
void OffLatticeCellBasedModel::doTrial(OffLatticeCell& cell)
{
    // store current state of energy/cell/num neighbors
    OffLatticeCell orig = cell;
    Energy preE = calculateHamiltonian(cell);
    unsigned preN = numNeighbors(cell);

    bool accepted, growth = attemptTrial(cell); // attempt the trial

    // auto reject if overlap or boundary violated
    if (checkOverlap(cell) || checkBoundary(cell))
    {
        cell = orig;
        accepted = false;
    }
    else
    {
        mCellPopulation.update(orig.coordinates(),
            cell.coordinates()); //update lattice with new cell position
        Energy postE = calculateHamiltonian(cell); // energy after update
        unsigned postN = numNeighbors(cell); // num neighbors after update

        // auto accept growth, otherwise accept based on energy change
        accepted = growth || acceptTrial(preE, postE, preN, postN);
        if (!accepted)
        {
            mCellPopulation.update(cell.coordinates(),
                orig.coordinates()); //revert cell position in lattice
            cell = orig;
        }            
    }
    if (growth) {cell.addToTrialRecord(accepted);} //record success/failure
}

// check if cell is ready to divide - execute division if it is
void OffLatticeCellBasedModel::checkMitosis(OffLatticeCell& cell)
{
    if (cell.readyToDivide())
    {
        OffLatticeCell daughter(cell.type());
        Point<double> old = cell.coordinates();
        cell.divide(daughter);
        mCellPopulation.update(old, cell.coordinates());
        mCellPopulation.insert(daughter.coordinates(), daughter);
    }
}

// check if cell overlaps any neighbors
bool OffLatticeCellBasedModel::checkOverlap(const OffLatticeCell& cell)
{
    double maxSearch = 4 * OL_PARAMS->maxRadius()
        + OL_PARAMS->maxTranslation(); // search radius
    LocalCellIterator it = // iterator around cell within radius
        mCellPopulation.lbegin(cell.coordinates(), maxSearch);
    LocalCellIterator endIt =
        mCellPopulation.lend(cell.coordinates(), maxSearch);*/

    for (; it != endIt; ++it)
    {
        if (cell != *it && cell.distance(*it) < 0) {return true;}
    }
    return false;
}

// check if cell extends past boundary
bool OffLatticeCellBasedModel::checkBoundary(const OffLatticeCell& cell)
{
    Point<double> origin(0,0);
    double b = mParams->boundary();

    // return true if cell is farther from center than the boundary line
    return (b > 0 &&
        (cell.centers().first.distance(origin) + cell.radius() > b 
        || cell.centers().second.distance(origin) + cell.radius() > b));
}

// execute growth trial - if max radius hit, enter mitosis
void OffLatticeCellBasedModel::growth(OffLatticeCell& cell)
{
    double growth = Random::uniform(0, maxGrowth(cell));
    double maxRadius = sqrt(2 * cell.type().size());

    cell.setRadius(std::min(maxRadius, cell.radius() + growth));
    if (cell.radius() == maxRadius) {cell.setPhase(MITOSIS);}
}

// execture translation trial (move cell)
void OffLatticeCellBasedModel::translation(OffLatticeCell& cell)
{
    double len = OL_PARAMS->maxTranslation() * sqrt(Random::uniform(0,1));
    double dir = Random::uniform(0, 2 * M_PI);
    cell.setCoordinates(Point<double>(cell.coordinates().x + len * cos(dir),
        cell.coordinates().y + len * sin(dir)));
}

// execute deformation trial, if max axis hit, get ready to divide
void OffLatticeCellBasedModel::deformation(OffLatticeCell& cell)
{
    double deform = Random::uniform(0, sqrt(cell.type().size())
        * OL_PARAMS->maxDeformation());
    double maxAxis = sqrt(16 * cell.type().size());

    cell.setAxisLength(std::min(maxAxis, cell.axisLength() + deform));

    if (cell.axisLength() == maxAxis)
    {
        cell.setReadyToDivide(true);
    }
}

// execute rotation trial
void OffLatticeCellBasedModel::rotation(OffLatticeCell& cell)
{
    double change = Random::uniform(-OL_PARAMS->maxRotation(),
        OL_PARAMS->maxRotation());
    cell.setAxisAngle(cell.axisAngle() + change / sqrt(cell.type().size()));
}

// record state of all cells
void OffLatticeCellBasedModel::recordPopulation()
{
    std::vector<double> current; // holds current population
    
    // loop through each cell, store info in current population
    CellIterator it = mCellPopulation.begin();
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
        current.push_back((*it).getTrialRecord());
    }

    mPopulationRecord.push_back(current); // add current pop to record
}

