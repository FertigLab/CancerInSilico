#ifndef CIS_SQUARE_LATTICE_H
#define CIS_SQUARE_LATTICE_H

#include <cmath>

#include "Lattice.h"

template <class T>
class SquareLattice : public Lattice<T>
{
private:
    // space between grid lines
    double mGridWidth;

    // hash a point in 2D space to a square grid
    GridPoint hash(const Point& pt) const
    {
        GridPoint hashedPt;
        hashedPt.x = ceil((fabs(pt.x) - mGridWidth / 2) / mGridWidth);
        hashedPt.y = ceil((fabs(pt.y) - mGridWidth / 2) / mGridWidth);

        hashedPt.x *= pt.x < 0 ? -1.0 : 1.0;
        hashedPt.y *= pt.y < 0 ? -1.0 : 1.0;

        return hashedPt;
    }

public:
    
    SquareLattice(double width) {mGridWidth = width;}
};

#endif
