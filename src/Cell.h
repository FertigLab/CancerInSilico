#ifndef CELL_HPP
#define CELL_HPP

#include "Point.h"
#include "Parameters.h"

#include <Rcpp.h>

class Cell
{

private:

    /* copy of parameters object */
    Parameters* mParameters;

    /* R object with all info about the cell's type */
    Rcpp::S4* mCellType;    

    /* geometric properties */
    Point mCoordinates;
    double mRadius;
	double mAxisLength, mAxisAngle;

    /* cell properties that can be effected by drug - must be R
       object so that R function can handl information, includes;
       - cycle length of cell */
    Rcpp::List mProperties;

    /* for each drug in the simulation, record the seed used to 
       calcualte effect (if random), that way even random effects
       can be inherited */
    Rcpp::List mDrugSeeds;

    /* phase in the cell cycle */
    enum CellPhase {I, M};
    CellPhase mPhase;

public:

    /* constructor used for initial cell population */
    Cell(Point, Parameters*, Rcpp::S4*);

    /* constructor for daughter cell, pass reference to parent */
    Cell(Point, const Cell&);

    /* getters */
    Point coordinates() const { return mCoordinates;}
    double radius() const { return mRadius;}
    double axisLength() const { return mAxisLength;}
    double axisAngle() const { return mAxisAngle;}
    double cycleLength() const { return mCycleLength;}
    double drugEffectSeed() const { return mDrugEffectSeed;}
    CellPhase phase() const { return mPhase;}
    Rcpp::S4* cellType() const { return mCellType;}
    Parameters* parameters() const {return mParameters;}
    double area() const { return M_PI * pow(mRadius, 2);}
    
    /* setters */   
    void setCoordinates(Point coords) { mCoordinates = coords;}
    void setRadius(double radius) { mRadius = radius;}
    void setAxisLength(double length) { mAxisLength = length;}
    void setAxisAngle(double angle) { mAxisAngle = angle;}
    void setCycleLength(double length) { mCycleLength = length;}
    void drugEffectSeed(double seed) { mDrugEffectSeed = seed;}
    void setPhase(CellPhase phase) { mPhase = phase;}

    /* operators */
    bool operator!=(const Cell&) const;
    bool operator==(const Cell&) const;
    double distance(const Cell&) const;

    /* undergo cell division, return daughter cell */
 	Cell divide();

    /* go to random point in the cell cycle */
    void gotoRandomCyclePoint();
};

#endif
