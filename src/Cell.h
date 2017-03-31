#ifndef CIS_CELL_H
#define CIS_CELL_H

#include <Rcpp.h>

#include "CellType.h"

enum CellPhase {INTERPHASE, MITOSIS, G0_PHASE, G1_PHASE, S_PHASE, G2_PHASE};

class Cell
{
private:

    // R object with all info about the cell's type
    CellType* mCellType;    

    // must be R object so that R functions can manipulate it
    double mCycleLength;

    // phase in the cell cycle
    CellPhase mPhase;

    // whether or not a drug has been applied
    std::vector<bool> mDrugApplied;

    // calculate post-division properties
    virtual void divide(Cell*, const Cell*) = 0;

public:

    // constructor for initial cells
    Cell(CellType*);

    // NOT COPY CONSTRUCTOR - used for daughter cell, pass pointer to parent
    Cell(const Cell*);

    /* getters */
    CellPhase phase() const {return mPhase;}
    CellType* cellType() const {return mCellType;}
    double cycleLength() const {return mCycleLength;}
    virtual double area() const = 0;
    bool drugApplied(unsigned i) {return mDrugApplied[i];}
    
    /* setters */   
    void setPhase(CellPhase phase) {mPhase = phase;}
    void setCellTpye(Rcpp::S4* type) {mCellType = type;}
    void setCycleLength(double len) {mCycleLength = len;}
};

#endif
