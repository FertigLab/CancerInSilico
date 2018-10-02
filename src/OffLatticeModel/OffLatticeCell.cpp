#include <cmath>
#include <Rcpp.h>

#include "OffLatticeCell.h"
#include "../Core/Random.h"

OffLatticeRadiusSolver OffLatticeCell::mSolver = OffLatticeRadiusSolver();

double sinTable[65536] = {0};
initTable initTableInstance;

// constructor
OffLatticeCell::OffLatticeCell(CellType type) : Cell(type)
{
    clearTrialRecord();
    mCoordinates = Point<double>(0.0, 0.0);
    setRadius(sqrt(mType.size()));
    mAxisAngle = Random::uniform(0, TWO_PI);
}

// set axis angle of cell
void OffLatticeCell::setAxisAngle(double angle)
{
    mAxisAngle = fmod(angle, TWO_PI);
    if (mAxisAngle < 0.0)
        mAxisAngle += TWO_PI;
}

// set radius of cell - adjust axis length accordingly
void OffLatticeCell::setRadius(double rad)
{
    mRadius = rad;
    mAxisLength = 2 * mRadius;
}

// set axis length of cell - adjust radius accordingly
void OffLatticeCell::setAxisLength(double len)
{
    if (len < sqrt(8 * mType.size()))
    {
        throw std::invalid_argument("adjusting axis on interphase cell");
    }
    mAxisLength = len;
    double adjusted = len / sqrt(mType.size());
    mRadius = sqrt(mType.size()) * mSolver.radius(adjusted);
}

// undergo cell division, return daughter cell
void OffLatticeCell::divide(OffLatticeCell& daughter)
{
    // update coordinates
    updateCenters();
    daughter.setCoordinates(centers().first);
    setCoordinates(centers().second);

    // update state of parent
    setRadius(sqrt(mType.size()));
    mAxisAngle = Random::uniform(0, TWO_PI);
    mPhase = INTERPHASE;
    mReadyToDivide = false;
    clearTrialRecord();
}

// go to random point in the cell cycle
void OffLatticeCell::gotoRandomCyclePoint()
{
    double chance = 1 - 2 / (mCycleLength + 2);
    if (Random::uniform(0,1) < chance) // random point of interphase
    {
        mPhase = INTERPHASE;
        setRadius(Random::uniform(sqrt(mType.size()),
            sqrt(2 * mType.size())));
    }
    else // random point of mitosis
    {
        mPhase = MITOSIS;
        setAxisLength(Random::uniform(sqrt(8 * mType.size()),
            sqrt(16 * mType.size())));
    }
}

// calculate the center of each side of the 'dumbell' shape
void OffLatticeCell::updateCenters() const
{
    // offset from true center in both dimensions
    double xOffset = ((0.5 * mAxisLength) - mRadius) * fastCos(mAxisAngle);
    double yOffset = ((0.5 * mAxisLength) - mRadius) * fastSin(mAxisAngle);

    // apply offset
    mCenters.first = Point<double>(mCoordinates.x + xOffset,
        mCoordinates.y + yOffset);
    mCenters.second = Point<double>(mCoordinates.x - xOffset,
        mCoordinates.y - yOffset);
}

// calculate distance between two cells (distance between edges)
double OffLatticeCell::distance(const OffLatticeCell& b) const
{
    // find smallest between two centers 
    updateCenters();
    b.updateCenters();

    // squared distances
    double AA = centers().first.distance2(b.centers().first);
    double AB = centers().first.distance2(b.centers().second);
    double BA = centers().second.distance2(b.centers().first);
    double BB = centers().second.distance2(b.centers().second);
    double minD2 = std::min(std::min(AA, AB), std::min(BA, BB));

    // return distance (between centers) minus the radii
    return sqrt(minD2) - radius() - b.radius();
}

bool OffLatticeCell::operator!=(const OffLatticeCell& other) const
{
    return coordinates() != other.coordinates();
}

bool OffLatticeCell::operator==(const OffLatticeCell& other) const
{
    return coordinates() == other.coordinates();
}

// clear accept/reject trail record
void OffLatticeCell::clearTrialRecord()
{
    mAcceptedTrials = 0;
    mTotalTrials = 0;
}

// add trial to record
void OffLatticeCell::addToTrialRecord(bool result)
{
    mTotalTrials += 1;
    if (result) {mAcceptedTrials += 1;}
}

// get trial record, return 1 if low number of trials done
double OffLatticeCell::getTrialRecord()
{
    return mTotalTrials < 5 ? 1 : mAcceptedTrials / mTotalTrials;
}

