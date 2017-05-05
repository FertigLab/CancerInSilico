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

    // monte carlo acceptance record
    std::vector<double> mAcceptRecord;
    std::vector<double> mTrialTimePoints;

public:

    // constructors
    OffLatticeCell(CellType);

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
    void setAxisAngle(double angle) {mAxisAngle = angle;}
    void setRadius(double);
    void setAxisLength(double);

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
