#ifndef CIS_OFF_LATTICE_CELL_HPP
#define CIS_OFF_LATTICE_CELL_HPP

#include <Rcpp.h>

#include "../Core/Point.h"
#include "../Core/Cell.h"
#include "OffLatticeRadiusSolver.h"

class OffLatticeCell : public Cell
{
private:

    // radius solver for cell geometry calculations
    static OffLatticeRadiusSolver mSolver;

    // geometric properties
    Point<double> mCoordinates;
    double mRadius;
    double mAxisLength, mAxisAngle;

public:

    // constructors
    OffLatticeCell(const CellType&);

    // getters
    Point<double> coordinates() const {return mCoordinates;}
    double radius() const {return mRadius;}
    double axisLength() const {return mAxisLength;}
    double axisAngle() const {return mAxisAngle;}
    double area() const  {return M_PI * pow(mRadius, 2);}

    // get the two centers of the cell
    std::pair< Point<double>, Point<double> > centers() const;    

    // setters
    void setCoordinates(Point<double> coords) {mCoordinates = coords;}
    void setRadius(double radius) {mRadius = radius;}
    void setAxisLength(double length) {mAxisLength = length;}
    void setAxisAngle(double angle) {mAxisAngle = angle;}

    // operators
    bool operator!=(const OffLatticeCell&) const;
    bool operator==(const OffLatticeCell&) const;
    double distance(const OffLatticeCell&) const;

    // process cell division
 	void divide(OffLatticeCell&);

    // go to random point in the cell cycle
    void gotoRandomCyclePoint();
};

#endif
