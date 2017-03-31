#include <cmath>
#include <Rcpp.h>

#include "Cell.h"
#include "Random.h"

// constructor for initial cells
Cell::Cell(Point coords, Rcpp::S4* type) : Cell(type)
{
    /* store variables */
    mCoordinates = coords;

    /* calculate geometric properties */
    mRadius = pow(type->slot("size"), 0.5);
    mAxisLength = 2 * mRadius;    
    mAxisAngle = Random::uniform(0, 2 * M_PI);
    mPhase = I;

    // go to random point in cycle if not synced
    if (!params->syncCellCycles())
    {
        gotoRandomCyclePoint();
    }

    /* get cycle length time */
    Rcpp::Function cycleLengthDist = type->slot("cycleLength");
    mCycleLength = Rcpp::as<double>(cycleLengthDist());
}

/* constructor for daughter cell, pass reference to parent */
Cell::Cell(Point coords, Cell& parent)
{
    /* store variables */
    mCoordinates = coords;
    mParameters = parent.parameters();
    mCellType = parent.cellType();

    /* calculate geometric properties */
    mRadius = pow(mCellType->slot("size"), 0.5);
    mAxisLength = 2 * mRadius;
    mAxisAngle = R::runif(0, 2 * M_PI);
    mPhase = I;
    
    /* check if cycle length is inherited */
    if (mCellType->slot("inheritCycleLength"))
    {
        mCycleLength = parent.cycleLength();
    }
    else
    {
        Rcpp::Function cycleLengthDist = type->slot("cycleLength");
        mCycleLength = Rcpp::as<double>(cycleLengthDist());
    }
}

// undergo cell division, return daughter cell
Cell Cell::divide()
{
    // store coordinates of daughter cell
    Point daughterCoords = Point(coordinates().x - cos(axisAngle()),
        coordinates().y - sin(axisAngle()));

    // update coordinates of parent cell
    mCoordinates = Point(coordinates().x + cos(axisAngle()),
        coordinates().y + sin(axisAngle())));

    // reset properties of parent
    mRadius = pow(mCellType->slot("size"), 0.5);
    mAxisLength = 2 * mRadius;
    mAxisAngle = R::runif(0, 2 * M_PI);
    mPhase = I;

    // return daughter cell
    return Cell(daughterCoords, *this);
}

// go to random point in the cell cycle
void Cell::gotoRandomCyclePoint()
{
    // random point of interphase
    if (R::runif(0,1) < 0.75)
    {
        mRadius = R::runif(mRadius, mRadius * pow(2, 0.5));
        mAxisLength = 2 * mRadius;
    }
    // random point of mitosis
    else
    {
        mPhase = M;
        mAxisLength = R::runif(2 * mRadius * pow(2, 0.5), 4 * mRadius);
        mRadius = mParams->GetRadius(mAxisLength);
    }
}

// calculate distance between two cells (distance between edges)
double ContinuumCell::distance(const ContinuumCell& b) const
{
    // centers of the dumbells for each cell 
    Point aCenters[2];
    Point bCenters[2];

    // offsets for each coordinate
    double aX = 0.5 * (axisLength() - radius()) * cos(axisAngle());
    double aY = 0.5 * (axisLength() - radius()) * sin(axisAngle());
    double bX = 0.5 * (b.axisLength() - b.radius()) * cos(b.axisAngle());
    double bY = 0.5 * (b.axisLength() - b.radius()) * sin(b.axisAngle());

    // get centers by offsetting coordinates in each direction
    aCenters[0] = Point(coordinates().x + aX, coordinates().y + aY);
    aCenters[1] = Point(coordinates().x - aX, coordinates().y - aY);
    bCenters[0] = Point(b.coordinates().x + bX, b.coordinates().y + bY);
    bCenters[1] = Point(b.coordinates().x - bX, b.coordinates().y - bY);

    // find smallest between two centers 
    double minDist = std::numeric_limits<double>::max();
    minDist = std::min(minDist, aCenters[0].distance(bCenters[0]));
    minDist = std::min(minDist, aCenters[0].distance(bCenters[1]));
    minDist = std::min(minDist, aCenters[1].distance(bCenters[0]));
    minDist = std::min(minDist, aCenters[1].distance(bCenters[1]));

    // return distance (between centers) minus the radii
    return minDist - radius() - b.radius();
}

bool ContinuumCell::operator!=(const ContinuumCell& other) const
{
    return coordinates() != other.coordinates();
}

bool ContinuumCell::operator==(const ContinuumCell& other) const
{
    return coordinates() == other.coordinates();
}

