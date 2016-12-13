#include <cmath>
#include <Rcpp.h>

#include "Cell.h"

Cell::Cell(Point coords, Parameters* params, Rcpp::S4* type)
{
    mCoordinates = coords;
    mParams = params;
    mCellType = type;

    mRadius = pow(type.slot("size"), 0.5);
    mAxisLength = 2 * mRadius;    
    mAxisAngle = R::runif(0, 2 * M_PI);

    Rcpp::Function cycleLength = type.slot("cycleLength");
    mGrowthRate = 
}

/* constructor for daughter cell, pass reference to parent */
Cell::Cell(Point coords, Cell& parent)
{
    mCoordinates = coords;
    mParams = parent.parameters();
    mCellType = parent.cellType();
    mAxisAngle = R::runif(0, 2 * M_PI);
}

Cell Cell::Divide()
{
    /* store coordinates of daughter cell */
    Point daughterCoords = Point(coordinates().x - cos(axisAngle()),
                                 coordinates().y - sin(axisAngle()));

    /* update coordinates of parent cell */
    setCoordinates(Point(coordinates().x + cos(axisAngle()),
                         coordinates().y + sin(axisAngle()));

    /* reset properties of parent */
    mAxisAngle = R::runif(0, 2 * M_PI);

    return Cell(daughterCoords, *this);
}

double Cell::distance(const Cell& b) const
{
    Point aCenters[2];
    Point bCenters[2];

    double aX = 0.5 * (axisLength() - radius()) * cos(axisAngle());
    double aY = 0.5 * (axisLength() - radius()) * sin(axisAngle());

    double bX = 0.5 * (b.axisLength() - b.radius()) * cos(b.axisAngle());
    double bY = 0.5 * (b.axisLength() - b.radius()) * sin(b.axisAngle());

    aCenters[0] = Point(coordinates().x + aX, coordinates().y + aY);
    aCenters[1] = Point(coordinates().x - aX, coordinates().y - aY);

    bCenters[0] = Point(b.coordinates().x + bX, b.coordinates().y + bY);
    bCenters[1] = Point(b.coordinates().x - bX, b.coordinates().y - bY);

    double minDist = std::numeric_limits<double>::max();

    for (unsigned int i = 0; i < 2; ++i)
    {
        for (unsigned int j = 0; j < 2; ++j)
        {
            minDist = std::min(minDist, aCenters[i].distance(bCenters[j]));
        }
    }

    return minDist - radius() - b.radius();

}

bool Cell::operator!=(const Cell& other) const
{
    return coordinates() != other.coordinates();
}

bool Cell::operator==(const Cell& other) const
{
    return coordinates() == other.coordinates();
}

