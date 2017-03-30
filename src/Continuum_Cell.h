#ifndef CIS_CELL_HPP
#define CIS_CELL_HPP

#include <Rcpp.h>

#include "Point.h"
#include "Cell.h"

class ContinuumCell : public Cell
{
private:

    // coordinates of cell center
    Point<double> mCoordinates;
    double mAxisLength, mAxisAngle;
    double mRadius;

public:

    // initial constructor
    ContinuumCell(Rcpp::S4*);

    // NOT COPY CONSTRUCTOR - used for daughter cell, pass pointer to parent
    ContinuumCell(const Cell* p) : Cell(p) {}

    /* getters */
    Point coordinates() const {return mCoordinates;}
    double radius() const {return mRadius;}
    double axisLength() const {return mAxisLength;}
    double axisAngle() const {return mAxisAngle;}
    double area() const  {return M_PI * pow(mRadius, 2);}
    
    /* setters */   
    void setCoordinates(Point coords) {mCoordinates = coords;}
    void setRadius(double radius) {mRadius = radius;}
    void setAxisLength(double length) {mAxisLength = length;}
    void setAxisAngle(double angle) {mAxisAngle = angle;}

    /* operators */
    bool operator!=(const ContinuumCell&) const;
    bool operator==(const ContinuumCell&) const;
    double distance(const ContinuumCell&) const;

    /* undergo cell division, return daughter cell */
 	Cell divide();

    /* go to random point in the cell cycle */
    void gotoRandomCyclePoint();
};

#endif
