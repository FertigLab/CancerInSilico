#ifndef CIS_CELL_H
#define CIS_CELL_H

#include <Rcpp.h>

enum CellPhase {I, M, G0, G1, S, G2};

class Cell
{
private:

    // R object with all info about the cell's type
    Rcpp::S4* mCellType;    

    // must be R object so that R functions can manipulate it
    double mCycleLength;

    // phase in the cell cycle
    CellPhase mPhase;

    // calculate post-division properties
    virtual void divide(Cell*, const Cell*) = 0;

public:

    // constructor for initial cells
    Cell(Rcpp::S4*);

    // NOT COPY CONSTRUCTOR - used for daughter cell, pass pointer to parent
    Cell(const Cell*);

    /* getters */
    CellPhase phase() const {return mPhase;}
    Rcpp::S4* cellType() const {return mCellType;}
    double cycleLength() const {return mCycleLength;}
    virtual double area() const = 0;
    
    /* setters */   
    void setPhase(CellPhase phase) {mPhase = phase;}
    void setCellTpye(Rcpp::S4* type) {mCellType = type;}
    void setCycleLength(double len) {mCycleLength = len;}
};

#endif
