#include <cmath>
#include <Rcpp.h>

#include "OffLatticeCell.h"
#include "../Core/Random.h"

OffLatticeRadiusSolver OffLatticeCell::mSolver = OffLatticeRadiusSolver();

// constructor
OffLatticeCell::OffLatticeCell(const CellType& type) : Cell(type)
{
    mCoordinates = Point<double>(0.0, 0.0);
    mRadius = sqrt(mType->size());
    mAxisLength = 2 * mRadius;
    mAxisAngle = Random::uniform(0, 2 * M_PI);
}

// undergo cell division, return daughter cell
void OffLatticeCell::divide(OffLatticeCell& daughter)
{
    // update coordinates
    daughter.setCoordinates(Point<double>(mCoordinates.x
        - cos(daughter.axisAngle()), mCoordinates.y
        - sin(daughter.axisAngle())));

    setCoordinates(Point<double>(mCoordinates.x
        + cos(axisAngle()), mCoordinates.y
        + sin(axisAngle())));

    // update stats of parent
    mRadius = sqrt(mType->size());
    mAxisLength = 2 * mRadius;
    mAxisAngle = R::runif(0, 2 * M_PI);
    mPhase = INTERPHASE;
}

// go to random point in the cell cycle
void OffLatticeCell::gotoRandomCyclePoint()
{
    // random point of interphase
    if (R::runif(0,1) < 0.75)
    {
        mPhase = INTERPHASE;
        mRadius = Random::uniform(sqrt(mType->size()),
            sqrt(2 * mType->size()));
        mAxisLength = 2 * mRadius;
    }
    // random point of mitosis
    else
    {
        mPhase = MITOSIS;
        mAxisLength = Random::uniform(2 * mRadius * sqrt(2), 4 * mRadius);
        mRadius = mSolver.radius(mAxisLength);
    }
}

std::pair< Point<double>, Point<double> > OffLatticeCell::centers() const
{
    std::pair< Point<double>, Point<double> > centers;
    
    double xOffset = 0.5 * (mAxisLength - mRadius) * cos(mAxisAngle);
    double yOffset = 0.5 * (mAxisLength - mRadius) * sin(mAxisAngle);

    centers.first = Point<double>(mCoordinates.x + xOffset,
        mCoordinates.y + yOffset);
    centers.second = Point<double>(mCoordinates.x - xOffset,
        mCoordinates.y - yOffset);

    return centers;
}

// calculate distance between two cells (distance between edges)
double OffLatticeCell::distance(const OffLatticeCell& b) const
{
    // find smallest between two centers 
    double minD = centers().first.distance(b.centers().first);
    minD = std::min(minD, centers().first.distance(b.centers().second));
    minD = std::min(minD, centers().second.distance(b.centers().first));
    minD = std::min(minD, centers().second.distance(b.centers().second));

    // return distance (between centers) minus the radii
    return minD - radius() - b.radius();
}

bool OffLatticeCell::operator!=(const OffLatticeCell& other) const
{
    return coordinates() != other.coordinates();
}

bool OffLatticeCell::operator==(const OffLatticeCell& other) const
{
    return coordinates() == other.coordinates();
}

