#include "OffLatticeCellBasedModel.h"

void OffLatticeCellBasedModel::updateRModel(Rcpp::S4* rModel)
{
    CellBasedModel::updateRModel(rModel);
    rModel->slot("acceptRecord") = Rcpp::wrap(mAcceptRecord);
}

// private helper function
static Point<double> getRandomPoint(double radius)
{
    double dist = Random::uniform(0,1);
    double ang = Random::uniform(0, 2 * M_PI);
    double x = radius * sqrt(dist) * cos(ang);
    double y = radius * sqrt(dist) * sin(ang);
    return Point<double>(x,y);
}

OffLatticeCellBasedModel::OffLatticeCellBasedModel(Rcpp::S4* rModel)
: CellBasedModel(rModel), mAcceptRecord(100, 0), mTrials(100, 0)
{
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
    double seedBoundary = sqrt(2 * area / M_PI);
    mDensity = 0.5;
    if (mParams->boundary() > 0)
    {
        mParams->setBoundary(sqrt(area / (M_PI * mParams->density())));
        seedBoundary = mParams->boundary();
        mDensity = mParams->density();
    }
    
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
    if (preE.second) {throw std::runtime_error("infinite hamiltonian");}
    unsigned preN = numNeighbors(cell);

    bool rejected, record = attemptTrial(cell);
    mDensity += (pow(cell.radius(), 2) - pow(orig.radius(), 2)) / 
        pow(mParams->boundary(), 2);
    Energy postE = calculateHamiltonian(cell);
    unsigned postN = numNeighbors(cell);

    if (checkOverlap(cell) || checkBoundary(cell) ||
    !acceptTrial(preE, postE, preN, postN))
    {
        mDensity -= (pow(cell.radius(), 2) - pow(orig.radius(), 2)) / 
            pow(mParams->boundary(), 2);
        cell = orig;
        rejected = true;
    }
    else
    {
        mCellPopulation.update(orig.coordinates(), cell.coordinates());
        rejected = false;
    }

    if (record && mParams->boundary() > 0)
    {
        int index = density() * 100;
        if (index < 0 || index > 99)
            {throw std::runtime_error("invalid index from density");}
        int n = ++mTrials[index];
        mAcceptRecord[index] *= (n - 1) / n;
        mAcceptRecord[index] += (rejected ? 0.0 : 1.0) / n;
    }
}

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

bool OffLatticeCellBasedModel::checkOverlap(const OffLatticeCell& cell)
{
    double maxSearch = 4 * OL_PARAMS->maxRadius();
    LocalCellIterator it =
        mCellPopulation.lbegin(cell.coordinates(), maxSearch);
    LocalCellIterator endIt =
        mCellPopulation.lend(cell.coordinates(), maxSearch);

    for (; it != endIt; ++it)
    {
        if (cell != *it && cell.distance(*it) < 0)
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
    return (b > 0 &&
        (cell.centers().first.distance(origin) + cell.radius() > b 
        || cell.centers().second.distance(origin) + cell.radius() > b));
}

void OffLatticeCellBasedModel::growth(OffLatticeCell& cell)
{
    double growth = Random::uniform(0, maxGrowth(cell));
    double sz = cell.type().size();
    cell.setRadius(std::min(sqrt(2 * sz), cell.radius() + growth));
    if (cell.radius() == sqrt(2 * sz)) {cell.setPhase(MITOSIS);}
}

void OffLatticeCellBasedModel::translation(OffLatticeCell& cell)
{
    double len = OL_PARAMS->maxTranslation() * sqrt(Random::uniform(0,1));
    double dir = Random::uniform(0, 2 * M_PI);
    cell.setCoordinates(Point<double>(cell.coordinates().x + len * cos(dir),
        cell.coordinates().y + len * sin(dir)));
}

void OffLatticeCellBasedModel::deformation(OffLatticeCell& cell)
{
    double deform = Random::uniform(0, sqrt(cell.type().size())
        * OL_PARAMS->maxDeformation());

    cell.setAxisLength(std::min(sqrt(16 * cell.type().size()),
        cell.axisLength() + deform));

    if (cell.axisLength() == sqrt(16 * cell.type().size()))
    {
        cell.setReadyToDivide(true);
    }
}

void OffLatticeCellBasedModel::rotation(OffLatticeCell& cell)
{
    double change = Random::uniform(-OL_PARAMS->maxRotation(),
        OL_PARAMS->maxRotation());
    cell.setAxisAngle(cell.axisAngle() + change / sqrt(cell.type().size()));
}

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
    }

    mPopulationRecord.push_back(current); // add current pop to record
}

