#ifndef CELL_HPP
#define CELL_HPP

#include "Point.h"
#include "Parameters.h"

#include <Rcpp.h>

class Cell
{

private:

    /* copy of parameters object */
    Parameters* mParams;

    /* R object with all info about the cell's type */
    Rcpp::S4* mCellType;    

    /* geometric properties */
    Point mCoordinates;
    double mRadius;
	double mAxisLength, mAxisAngle;

    /* growth rate */
    double mGrowthRate;
    
    /* phase in the cell cycle */
    enum CellPhase {I, M, G0, G1, S, G2};
    CellPhase mPhase;

public:

    /* constructor used for initial cell population */
    Cell(Point, Parameters*, Rcpp::S4*);

    /* constructor for daughter cell, pass reference to parent */
    Cell(Point, Cell&);

    /* getters */
    Point coord() { return mCoordinates;}
    double radius() { return mRadius;}
    double axisLength() { return mAxisLength;}
    double axisAngle() { return mAxisAngle;}
    double growthRate() { return mGrowthRate;}
    CellPhase phase() { return mPhase;}
    Rcpp::S4* cellType() { return mCellType;}
    Parameters* parameters() {return mParameters;}
    double area() { return M_PI * pow(mRadius, 2);}
    
    /* setters */   
    void setCoordinates(Point c) { mCoordinates = c;}
    void setRadius(double r) { mRadius = r;}
    void setAxisLength(double al) { mAxisLength = al;}
    void setAxisAngle(double aa) { mAxisAngle = aa;}
    void setGrowthRate(double gr) { mGrowthRate = gr;}
    void setPhase(CellPhase p) { mPhase = p;}

    /* undergo cell division, return daughter cell */
 	Cell Divide();

    /* operator definitions */
    bool operator!=(const Cell&) const;
    bool operator==(const Cell&) const;
    double distance(const Cell&) const;

};

#endif
