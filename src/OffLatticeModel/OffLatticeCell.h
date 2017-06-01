#ifndef CIS_OFF_LATTICE_CELL_HPP
#define CIS_OFF_LATTICE_CELL_HPP

// used for fast sin calculation

#define TWO_PI 6.28318530718
#define HALF_PI 1.57079632679
#define FIVE_HALF_PI 7.85398163397

extern double sinTable[65536];
struct initTable
{
    initTable()
    {
        for (int i = 0; i < 65536; ++i)
        {
            sinTable[i] = sin(i * TWO_PI / 65536);
        }
    }
};

inline double fastSin(double x)
{
    return sinTable[(int) (x * 65535 / TWO_PI)];
}

inline double fastCos(double x)
{
    return x < HALF_PI ? fastSin(HALF_PI - x) : fastSin(FIVE_HALF_PI - x);
}

// cell implementation for off lattice models

#include <Rcpp.h>

#include "../Core/Point.h"
#include "../Core/Cell.h"
#include "OffLatticeRadiusSolver.h"

typedef std::pair< Point<double>, Point<double> > PointPair;

class OffLatticeCell : public Cell
{
private:

    // radius solver for cell geometry calculations
    static OffLatticeRadiusSolver mSolver;

    // geometric properties
    Point<double> mCoordinates;
    double mRadius;
    double mAxisLength, mAxisAngle;

    // percent of trials accepted
    double mAcceptedTrials, mTotalTrials;
    mutable PointPair mCenters;

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
    void updateCenters() const;
    PointPair centers() const {return mCenters;}

    // setters
    void setCoordinates(Point<double> coords) {mCoordinates = coords;}
    void setAxisAngle(double);
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

    // update trial record
    void clearTrialRecord();
    void addToTrialRecord(bool);
    double getTrialRecord();
};

#endif
